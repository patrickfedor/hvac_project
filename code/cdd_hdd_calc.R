#### Set up ####

library(tidyverse)
library(stars)
library(ncdf4)
library(furrr)
library(tigris)
library(sf)
library(mgcv)
library(RColorBrewer)

plan(multisession) #enable parallel processing

#### Process Data ####

# Download from bucket
system("gsutil ls gs://clim_data_reg_useast1/era5/daily_aggregates/2m_temperature/", intern = TRUE) %>%
  str_subset("200[1-9]|201[0-9]|202[0-3]") %>%  # full range 2001–2023
  walk(~ system(paste("gsutil cp", shQuote(.x), "/mnt/pers_disk/hvac")))


# system("gsutil ls gs://clim_data_reg_useast1/era5/daily_aggregates/2m_temperature/ | grep '202[0-3]' > /mnt/pers_disk/filtered_files.txt") # alternative way - edit grep for correct years
# system("cat /mnt/pers_disk/filtered_files.txt | xargs -I {} gsutil cp {} /mnt/pers_disk/hvac")


# Load NetCDF files
nc_files <- list.files("/mnt/pers_disk/hvac", pattern = "\\.nc$", full.names = TRUE)

# Define CDD & HDD functions
calculate_cdd <- function(temp) {
  cdd <- ifelse(temp >= 65, temp - 65, 0)
  r <- c(CDD = cdd)
  return(r) 
}

calculate_hdd <- function(temp) {
  hdd <- ifelse(temp < 65, 65 - temp, 0)
  r <- c(HDD = hdd)
  return(r) 
}


# Set up directories
output_dir <- "/mnt/pers_disk/hvac_processed"
# fs::dir_create(output_dir)

bucket_path <- "gs://clim_data_reg_useast1/misc_data/temporary/"


# Loop to process each NetCDF file
for (file in nc_files) {

  # Read NetCDF file
  nc_data <- read_ncdf(file)
  
  # Convert units
  nc_data <- nc_data %>%
    mutate(t2m = t2m %>% units::set_units(degF))

  # Apply CDD/HDD calculations across space & time
  cdd_result <- st_apply(nc_data, 
                         MARGIN = c("longitude", "latitude"),
                         FUN = calculate_cdd,
                         .fname = 'CDD') %>%
    aperm(c(3,2,1))
  
  hdd_result <- st_apply(nc_data, 
                         MARGIN = c("longitude", "latitude"),
                         FUN = calculate_hdd,
                         .fname = 'CDD') %>%
    aperm(c(3,2,1))

  # Define output 
  file_basename <- tools::file_path_sans_ext(basename(file))
  cdd_output_file <- file.path(output_dir, paste0(file_basename, "_CDD.nc"))
  hdd_output_file <- file.path(output_dir, paste0(file_basename, "_HDD.nc"))

  # Save results
  write_stars(cdd_result, cdd_output_file)
  write_stars(hdd_result, hdd_output_file)
  
}


# List all processed NetCDF files 
processed_files <- list.files(output_dir, pattern = "\\.nc$", 
                              full.names = TRUE)

# Upload to bucket
for (file in processed_files) {
  system(paste("gsutil mv", file, file.path(bucket_path, 
                                            basename(file))))
}

#### CSV Data ####

# List energy data in Cloud bucket
csv_files <- system("gsutil ls gs://clim_data_reg_useast1/misc_data/temporary//HVAC", intern = T)

# Copy into local disk (pers_disk)
destination_dir <- "/mnt/pers_disk/hvac"
# fs::dir_create(destination_dir) #uncomment if folder does not exist
walk(csv_files, ~ system(paste("gsutil cp", .x, destination_dir)))

# List downloaded CSV files
csv_files <- list.files("/mnt/pers_disk/hvac", pattern = "_data.csv$",
                        full.names = TRUE)

# Read and covert to df
csv_data <- map_dfr(csv_files, ~ read_csv(.x))

# Make date column into date object and filter to be only years for analysis
csv_data <- csv_data %>% mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
  filter(date <= as.Date("2023-12-31"))



#### CDD Analysis ####


# Download U.S. states shapefile
us_states <- states(cb= F) %>%
  st_transform(4326)  

# Filter for your states of interest
target_states <- c("TX", "NJ", "DE", "WV", "VA", "OH", "PA", "MD")
state_shapes <- us_states %>% 
  filter(STUSPS %in% target_states)

# List and filter relevant CDD from cloud bucket
nc_files <- system("gsutil ls gs://clim_data_reg_useast1/misc_data/temporary//*", intern = T)

cdd_files <- grep("CDD.nc$", nc_files, value = T)  
# cdd_files <- grep("-(06|07|08)-", cdd_files, value = T)  #specific months 

# Set up directories
destination_dir <- "/mnt/pers_disk/hvac"
# fs::dir_create(destination_dir)

# Download the files to local disk 
walk(cdd_files, ~ system(paste("gsutil cp", .x, destination_dir)))


# List downloaded CDD NetCDF files
nc_files <- list.files("/mnt/pers_disk/hvac", pattern = "CDD.nc$",
                       full.names = TRUE)

# Read into R
cdd_stars <- future_map(nc_files, ~ read_stars(.x, proxy = TRUE), .progress = TRUE)

# Function to extract and aggregate CDD for each state
extract_state_cdd <- function(cdd, states) {
  cdd %>%
    st_crop(states) %>%  
    aggregate(by = states, FUN = mean, na.rm = TRUE) 
}

# Apply function across all months
state_cdd <- future_map(cdd_stars, extract_state_cdd, states = state_shapes, .progress = TRUE)

# Create a time index for bands in stars object
date_vector <- seq(as.Date('2001-01-01'), as.Date('2023-12-31'), by = "1 day")

# date_vector <- date_vector[month(date_vector) %in% c(6, 7, 8)] # uncomment for specific months

# Convert to CDD data to df
cdd_df <- map_dfr(state_cdd, function(s) {
  s %>% st_as_sf() %>% 
  mutate(state = state_shapes$NAME) %>% 
  st_drop_geometry() %>% 
  pivot_longer(-state, names_to = 'day', values_to = "CDD")
})

cdd_df <- cdd_df %>% 
  arrange(state) %>% 
  mutate(day = rep(date_vector, 8)) #add date vector into df

# Aggregate CDD df to monthly
cdd_monthly <- cdd_df %>%
  mutate(year = year(day), month = month(day)) %>%
  group_by(state, year, month) %>%
  summarize(cdd_mean = mean(CDD, na.rm = TRUE), 
            cdd_sum = sum(CDD, na.rm = T), .groups = "drop")

# Save data

write_csv(cdd_df, "/mnt/pers_disk/hvac_processed/cdd_df.csv")

write_csv(cdd_monthly, "/mnt/pers_disk/hvac_processed/cdd_monthly.csv")

system('gsutil mv ')
##### Plots ####

# Base df used for plotting
cdd_df <- read_csv('/mnt/pers_disk/hvac_processed/cdd_df.csv')

# Define a color palette
colors <- brewer.pal(12, "Set3") 

# Overall CDD by month across years
monthly_cdd <- cdd_df %>%
  mutate(year = lubridate::year(day), month = lubridate::month(day)) %>%
  group_by(year, month) %>%
  summarise(monthly_cdd = sum(CDD, na.rm = TRUE))

# Plot the time series for overall CDD trend with month names
ggplot(monthly_cdd, aes(x = year, y = monthly_cdd, group = month, color = factor(month, labels = month.name))) +
  geom_line() +
  scale_color_manual(values = colors) +
  labs(title = "Monthly CDD Trends (2001-2023)", x = "Year", y = "Cooling Degree Days (CDD)", color = "Month") +
  theme_minimal()


# Daily CDD Trends by Year
cdd_trend_yearly <- cdd_df %>%
  mutate(year = lubridate::year(day), doy = lubridate::yday(day)) %>%
  group_by(year, doy) %>%
  summarise(mean_cdd = mean(CDD, na.rm = TRUE), .groups = "drop")

ggplot(cdd_trend_yearly, aes(x = doy, y = mean_cdd, color = factor(year))) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Trend lines
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  labs(title = "Trends in Cooling Degree Days (CDD) Over the Year",
       subtitle = "How has CDD changed over different years?",
       x = "Day of Year",
       y = "Mean Cooling Degree Days (CDD)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Mean Monthly CDD Trends by Month
cdd_trend <- cdd_df %>%
  mutate(year = lubridate::year(day), month = lubridate::month(day, label = TRUE)) %>%
  group_by(year, month) %>%
  summarise(mean_cdd = mean(CDD, na.rm = TRUE), .groups = "drop")

ggplot(cdd_trend, aes(x = year, y = mean_cdd, group = month, color = month)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "lm") + # Trend line
  facet_wrap(~ month, scales = "fixed", ncol = 4) +  
  labs(title = "Trend of Cooling Degree Days (CDD) Over Time",
       subtitle = "Do certain months show an increasing trend?",
       x = "Year",
       y = "Mean CDD",
       color = "Month") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Heatmap of Daily Mean CDD by Year
ggplot(cdd_trend_yearly, aes(x = doy, y = factor(year), fill = mean_cdd)) +
  geom_tile() +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  scale_fill_viridis_c(option = "magma", name = "Mean CDD") +
  labs(title = "Heatmap of Cooling Degree Days (CDD) Over the Year",
       subtitle = "Do certain periods show increasing CDD over time?",
       x = "Day of Year",
       y = "Year") +
  theme_minimal()

# Boxplots of Monthly CDD by Year (not a good visual)
ggplot(cdd_trend, aes(x = factor(year), y = mean_cdd, fill = factor(month))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ month, scales = "free_y", ncol = 4) +
  scale_fill_viridis_d() +
  labs(title = "Monthly Cooling Degree Days (CDD) Over Time",
       subtitle = "Are CDD values becoming more extreme?",
       x = "Year",
       y = "Mean CDD",
       fill = "Month") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Daily CDD Anomaly Over Time
cdd_trend_yearly <- cdd_trend_yearly %>%
  group_by(doy) %>%
  mutate(cdd_anomaly = mean_cdd - mean(mean_cdd, na.rm = TRUE))

ggplot(cdd_trend_yearly, aes(x = doy, y = cdd_anomaly, color = factor(year))) +
  geom_line(size = 1) +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  labs(title = "Cooling Degree Days (CDD) Anomaly Over the Year",
       subtitle = "How much do recent years deviate from the historical average?",
       x = "Day of Year",
       y = "CDD Anomaly",
       color = "Year") +
  theme_minimal()

#CDD Anomaly Heatmap
ggplot(cdd_trend_yearly, aes(x = doy, y = factor(year), fill = cdd_anomaly)) +
  geom_tile() +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "CDD Anomaly") +
  labs(title = "Heatmap of CDD Anomalies Over Time",
       subtitle = "Red = Higher than average | Blue = Lower than average",
       x = "Day of Year",
       y = "Year") +
  theme_minimal()


# Smoothed CDD Anomaly Trends by Period
cdd_trend_yearly <- cdd_trend_yearly %>%
  group_by(doy) %>%
  mutate(cdd_anomaly = mean_cdd - mean(mean_cdd, na.rm = TRUE)) %>%
  mutate(period = case_when(
    year <= 2010 ~ "2001-2010",
    year > 2010 & year <= 2020 ~ "2011-2020",
    year > 2020 ~ "2021-2023"
  ))

ggplot(cdd_trend_yearly, aes(x = doy, y = cdd_anomaly, color = period)) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, size = 1.5) +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  labs(title = "CDD Anomaly Trends: Early vs. Recent Years",
       subtitle = "Are recent years warmer compared to earlier years?",
       x = "Day of Year",
       y = "CDD Anomaly",
       color = "Period") +
  theme_minimal()


##### Correlation #####

# Set up CSV data
energy_monthly <- csv_data %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(state, year, month) %>% 
  summarize(total_generation_gwh = sum(generation_gwh, na.rm = T), .groups = 'drop')

###### Statewide Correlation (Total Power vs. CDD) #######

# Merge monthly CDD data with energy generation data
merged_data <- left_join(cdd_monthly, energy_monthly, by = c("state", "year", "month"))

# Correlation across all states and months
cor_result <- cor(merged_data$total_generation_gwh, merged_data$cdd_mean, use = "complete.obs")
print(cor_result)

# Correlation computed per state
statewise_correlation <- merged_data %>%
  group_by(state) %>%
  summarize(correlation = cor(total_generation_gwh, cdd_mean, use = "complete.obs"))
print(statewise_correlation)

###### Correlation by Energy source ######

# Monthly totals
energy_monthly_by_source <- csv_data %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(state, year, month, variable) %>%
  summarize(generation = sum(generation_gwh, na.rm = TRUE), .groups = "drop")

# Merge with CDD data
merged_data2 <- left_join(cdd_monthly, energy_monthly_by_source, by = c("state", "year", "month"))

# Compute correlation by energy type
cor_by_source <- merged_data2 %>%
  group_by(variable) %>%
  summarize(correlation = cor(generation, cdd_mean, use = "complete.obs"))
print(cor_by_source)

# Compute correlation by state and energy type
cor_by_source_state <- merged_data2 %>%
  group_by(state, variable) %>%
  summarize(correlation = cor(generation, cdd_mean, use = "complete.obs"))
print(cor_by_source_state)

# Rsquared Values #

# Loess (default geom_smooth method)
r2_values <- merged_data2 %>%
  group_by(variable) %>%
  summarise(
    r_squared = {
      model <- loess(generation ~ cdd_mean, data = pick(everything()))
      1 - sum(model$residuals^2) / sum((generation - mean(generation))^2)
    }
  )

# # Linear Rsquared values if using 'lm' function
# r2_values <- merged_data2 %>%
#   group_by(variable) %>%
#   summarise(
#     r_squared = summary(lm(generation ~ hdd_mean, 
#                            data = pick(everything())))$r.squared
#   )


# Attache R values to main df
merged_data2 <- merged_data2 %>%
  left_join(r2_values, by = "variable")

## Plots

# CDD vs Generation by Energy Source (with R labels)
ggplot(merged_data2, aes(x = cdd_mean, y = generation, color = variable)) +
  geom_point(alpha = 0.6) +  
  geom_smooth(color = 'black') + 
  facet_wrap(~ variable, scales = "free") +  
  geom_text(data = r2_values, aes(x = Inf, y = Inf, label = paste0("R² = ", round(r_squared, 2))),
            hjust = 1.1, vjust = 1.2, inherit.aes = FALSE, size = 3) +  
  labs(title = "Correlation between CDD and Power Generation by Energy Source",
       x = "CDD",
       y = "Power Generation (GWh)") +
  theme_minimal() +
  theme(legend.position = "none")


# # Linear regression for each energy source
# ggplot(merged_data2, aes(x = hdd_mean, y = generation, color = variable)) +
#   geom_point(alpha = 0.6) +  
#   geom_smooth(color = 'black') +  
#   facet_wrap(~ variable, scales = "free") +  
#   labs(title = "Correlation between HDD and Power Generation by Energy Source",
#        x = "Heating Degree Days (HDD)",
#        y = "Power Generation (GWh)") +
#   theme_minimal() +
#   theme(legend.position = "none")

# Linear regression for statewide analysis
ggplot(merged_data, aes(x = cdd_mean, y = total_generation_gwh)) +
  geom_point(alpha = 0.5, aes(color = state)) +  # Scatter plot with transparency
  geom_smooth(method = "lm") +  # Linear regression line
  facet_wrap(~ state, scales = "free") +  # Facet by state to show individual regressions
  labs(title = "Statewide Linear Regression: CDD vs. Power Generation",
       x = "CDD",
       y = "Power Generation (GWh)") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide the legend for color, as each plot is by state


# Bar Chart correlation by State
ggplot(statewise_correlation, aes(x = reorder(state, correlation), y = correlation, fill = correlation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  
  labs(title = "Statewide Correlation between CDD and Power Generation",
       x = "State", y = "Correlation with CDD") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +  # Color gradient based on correlation
  theme_minimal()

# Correlation by Energy Source
ggplot(cor_by_source, aes(x = reorder(variable, correlation), y = correlation, fill = correlation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  
  labs(title = "Statewide Correlation between CDD and Energy Source",
       x = "Energy Source", y = "Correlation with CDD") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +  # Color gradient based on correlation
  theme_minimal()


















#### HDD Analysis ####

# Same as CDD just was easier to break the analysis into two parts

# Download U.S. states shapefile
us_states <- states(cb= F) %>%
  st_transform(4326)  

# Filter for states of interest
target_states <- c("TX", "NJ", "DE", "WV", "VA", "OH", "PA", "MD")
state_shapes <- us_states %>% 
  filter(STUSPS %in% target_states)


# # List and filter relevant HDD from cloud bucket
nc_files <- system("gsutil ls gs://clim_data_reg_useast1/misc_data/temporary//*", intern = T)

hdd_files <- grep("HDD.nc$", nc_files, value = T)

# hdd_files <- grep("-(01|02|12)-", hdd_files, value = T) #specific months

# Set up directories
destination_dir <- "/mnt/pers_disk/hvac"
# fs::dir_create(destination_dir)

# Download the files to local disk 
walk(hdd_files, ~ system(paste("gsutil cp", .x, destination_dir)))


# List downloaded HDD NetCDF files in directory
nc_files <- list.files("/mnt/pers_disk/hvac", pattern = "HDD.nc$",
                       full.names = TRUE)

# Read into R
hdd_stars <- future_map(nc_files, ~ read_stars(.x, proxy = TRUE), .progress = TRUE)

# Function to extract and aggregate HDD for each state
extract_state_hdd <- function(hdd, states) {
  hdd %>%
    st_crop(states) %>%  
    aggregate(by = states, FUN = mean, na.rm = TRUE) 
}

# Apply function across all months
state_hdd <- future_map(hdd_stars, extract_state_hdd, states = state_shapes, .progress = TRUE)

# Create a time vector for bands in stars object
date_vector <- seq(as.Date('2001-01-01'), as.Date('2023-12-31'), by = "1 day")

# date_vector2 <- date_vector[month(date_vector) %in% c(1, 2, 12)] # for specific months

# Create HDD data to df
hdd_df <- map_dfr(state_hdd, function(s) {
  s %>% st_as_sf() %>% 
    mutate(state = state_shapes$NAME) %>% 
    st_drop_geometry() %>% 
    pivot_longer(-state, names_to = 'day', values_to = "HDD")
})

hdd_df <- hdd_df %>% 
  arrange(state) %>% 
  mutate(day = rep(date_vector, 8))

# Aggregate HDD df to monthly
hdd_monthly <- hdd_df %>%
  mutate(year = year(day), month = month(day)) %>%
  group_by(state, year, month) %>%
  summarize(hdd_mean = mean(HDD, na.rm = TRUE), .groups = "drop")


# Save data
write_csv(hdd_df, "/mnt/pers_disk/hvac_processed/hdd_df.csv")

write_csv(hdd_monthly, "/mnt/pers_disk/hvac_processed/hdd_monthly.csv")

##### Plots #####

# Base df used for plotting
hdd_df <- read_csv('/mnt/pers_disk/hvac_processed/hdd_df.csv')

# Define a color palette 
colors <- brewer.pal(12, "Set3") 

# Overall HDD by month across years
monthly_hdd <- hdd_df %>%
  mutate(year = lubridate::year(day), month = lubridate::month(day)) %>%
  group_by(year, month) %>%
  summarise(monthly_hdd = sum(HDD, na.rm = TRUE))

# Plot the time series for overall HDD trend with month names
ggplot(monthly_hdd, aes(x = year, y = monthly_hdd, group = month, color = factor(month, labels = month.name))) +
  geom_line() +
  scale_color_manual(values = colors) +
  labs(title = "Monthly HDD Trends (2001-2023)", x = "Year", y = "Heating Degree Days (HDD)", color = "Month") +
  theme_minimal()

# Daily HDD Trends by Yea
hdd_trend_yearly <- hdd_df %>%
  mutate(year = lubridate::year(day), doy = lubridate::yday(day)) %>%
  group_by(year, doy) %>%
  summarise(mean_hdd = mean(HDD, na.rm = TRUE), .groups = "drop")

ggplot(hdd_trend_yearly, aes(x = doy, y = mean_hdd, color = factor(year))) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Trend lines
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  labs(title = "Trends in Heating Degree Days (HDD) Over the Year",
       subtitle = "How has HDD changed over different years?",
       x = "Day of Year",
       y = "Mean Heating Degree Days (HDD)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Mean Monthly HDD Trends by Month
hdd_trend <- hdd_df %>%
  mutate(year = lubridate::year(day), month = lubridate::month(day, label = TRUE)) %>%
  group_by(year, month) %>%
  summarise(mean_hdd = mean(HDD, na.rm = TRUE), .groups = "drop")

ggplot(hdd_trend, aes(x = year, y = mean_hdd, group = month, color = month)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "lm") + # Trend line
  facet_wrap(~ month, scales = "fixed", ncol = 4) +  
  labs(title = "Trend of Heating Degree Days (HDD) Over Time",
       subtitle = "Do certain months show an increasing trend?",
       x = "Year",
       y = "Mean HDD",
       color = "Month") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Heatmap of Daily Mean HDD by Year
ggplot(hdd_trend_yearly, aes(x = doy, y = factor(year), fill = mean_hdd)) +
  geom_tile() +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  scale_fill_viridis_c(option = "magma", name = "Mean HDD") +
  labs(title = "Heatmap of Heating Degree Days (HDD) Over the Year",
       subtitle = "Do certain periods show increasing HDD over time?",
       x = "Day of Year",
       y = "Year") +
  theme_minimal()

# Boxplot of Monthly HDD by Year
ggplot(hdd_trend, aes(x = factor(year), y = mean_hdd, fill = factor(month))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ month, scales = "free_y", ncol = 4) +
  scale_fill_viridis_d() +
  labs(title = "Monthly Heating Degree Days (HDD) Over Time",
       subtitle = "Are HDD values becoming more extreme?",
       x = "Year",
       y = "Mean CDD",
       fill = "Month") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Daily CDD Anomaly Over Time
hdd_trend_yearly <- hdd_trend_yearly %>%
  group_by(doy) %>%
  mutate(hdd_anomaly = mean_hdd - mean(mean_hdd, na.rm = TRUE))

ggplot(hdd_trend_yearly, aes(x = doy, y = hdd_anomaly, color = factor(year))) +
  geom_line(linewidth = 1) +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  labs(title = "Heating Degree Days (HDD) Anomaly Over the Year",
       subtitle = "How much do recent years deviate from the historical average?",
       x = "Day of Year",
       y = "HDD Anomaly",
       color = "Year") +
  theme_minimal()

# Hdd Anomaly Heatmap
ggplot(hdd_trend_yearly, aes(x = doy, y = factor(year), fill = hdd_anomaly)) +
  geom_tile() +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "HDD Anomaly") +
  labs(title = "Heatmap of HDD Anomalies Over Time",
       subtitle = "Red = Higher than average | Blue = Lower than average",
       x = "Day of Year",
       y = "Year") +
  theme_minimal()

# Smoothed HDD Anomaly Trends by Period
hdd_trend_yearly <- hdd_trend_yearly %>%
  group_by(doy) %>%
  mutate(hdd_anomaly = mean_hdd - mean(mean_hdd, na.rm = TRUE)) %>%
  mutate(period = case_when(
    year <= 2010 ~ "2001-2010",
    year > 2010 & year <= 2020 ~ "2011-2020",
    year > 2020 ~ "2021-2023"
  ))

ggplot(hdd_trend_yearly, aes(x = doy, y = hdd_anomaly, color = period)) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, size = 1.5) +  
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  
  labs(title = "HDD Anomaly Trends: Early vs. Recent Years",
       subtitle = "Are recent years warmer compared to earlier years?",
       x = "Day of Year",
       y = "HDD Anomaly",
       color = "Period") +
  theme_minimal()







##### Correlation #####

# Set up CSV data
energy_monthly <- csv_data %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(state, year, month) %>% 
  summarize(total_generation_gwh = sum(generation_gwh, na.rm = T), .groups = 'drop')

###### Statewide Correlation (Total Power vs. CDD) #######

# Merge monthly HDD data with energy generation data
merged_data <- merged_data <- left_join(hdd_monthly, energy_monthly, by = c("state", "year", "month"))

# Correlation across all states and months
cor_result <- cor(merged_data$total_generation_gwh, merged_data$hdd_mean, use = "complete.obs")
print(cor_result)

# Correlation computed per state
statewise_correlation <- merged_data %>%
  group_by(state) %>%
  summarize(correlation = cor(total_generation_gwh, hdd_mean, use = "complete.obs"))
print(statewise_correlation)

###### By Energy source ######

# Monthly Totals
energy_monthly_by_source <- csv_data %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(state, year, month, variable) %>%
  summarize(generation = sum(generation_gwh, na.rm = TRUE), .groups = "drop")

# Merge with HDD data
merged_data2 <- left_join(hdd_monthly, energy_monthly_by_source, by = c("state", "year", "month"))

# Compute correlation by energy type
cor_by_source <- merged_data2 %>%
  group_by(variable) %>%
  summarize(correlation = cor(generation, hdd_mean, use = "complete.obs"))

# Compute correlation by state and energy type
cor_by_source_state <- merged_data2 %>%
  group_by(state, variable) %>%
  summarize(correlation = cor(generation, hdd_mean, use = "complete.obs"))
print(cor_by_source)

# Rsquared Values #

# Loess (default geom_smooth method)
r2_values <- merged_data2 %>%
  group_by(variable) %>%
  summarise(
    r_squared = {
      model <- loess(generation ~ hdd_mean, data = pick(everything()))
      1 - sum(model$residuals^2) / sum((generation - mean(generation))^2)
    }
  )

# # Linear Rsquared
# r2_values <- merged_data2 %>%
#   group_by(variable) %>%
#   summarise(
#     r_squared = summary(lm(generation ~ hdd_mean, 
#                            data = pick(everything())))$r.squared
#   )


# Attach R values to main df
merged_data2 <- merged_data2 %>%
  left_join(r2_values, by = "variable")

## Plots

# HDD vs Generation by Energy Source (with R labels)
ggplot(merged_data2, aes(x = hdd_mean, y = generation, color = variable)) +
  geom_point(alpha = 0.6) +  
  geom_smooth(color = 'black') + 
  facet_wrap(~ variable, scales = "free") +  
  geom_text(data = r2_values, aes(x = Inf, y = Inf, label = paste0("R² = ", round(r_squared, 2))),
            hjust = 1.1, vjust = 1.2, inherit.aes = FALSE, size = 3) +  
  labs(title = "Correlation between HDD and Power Generation by Energy Source",
       x = "HDD",
       y = "Power Generation (GWh)") +
  theme_minimal() +
  theme(legend.position = "none")


# Linear regression for each energy source
ggplot(merged_data2, aes(x = hdd_mean, y = generation, color = variable)) +
  geom_point(alpha = 0.6) +  
  geom_smooth(color = 'black') +  
  facet_wrap(~ variable, scales = "free") +  
  labs(title = "Correlation between HDD and Power Generation by Energy Source",
       x = "Heating Degree Days (HDD)",
       y = "Power Generation (GWh)") +
  theme_minimal() +
  theme(legend.position = "none")

# Linear regression for statewise analysis
ggplot(merged_data, aes(x = hdd_mean, y = total_generation_gwh)) +
  geom_point(alpha = 0.5, aes(color = state)) +  # Scatter plot with transparency
  geom_smooth(method = "lm") +  # Linear regression line
  facet_wrap(~ state, scales = "free") +  # Facet by state to show individual regressions
  labs(title = "Statewide Linear Regression: HDD vs. Power Generation",
       x = "Heating Degree Days (HDD)",
       y = "Power Generation (GWh)") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide the legend for color, as each plot is by state


# Bar Chart correlation by State
ggplot(statewise_correlation, aes(x = reorder(state, correlation), y = correlation, fill = correlation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  
  labs(title = "Statewide Correlation between HDD and Power Generation",
       x = "State", y = "Correlation with HDD") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +  # Color gradient based on correlation
  theme_minimal()

# Correlation plot Energy Source
ggplot(cor_by_source, aes(x = reorder(variable, correlation), y = correlation, fill = correlation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  
  labs(title = "Statewide Correlation between HDD and Energy Source",
       x = "Energy Source", y = "Correlation with HDD") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +  # Color gradient based on correlation
  theme_minimal()