#### Set up ####

library(tidyverse)
library(stars)
library(ncdf4)
library(furrr)
library(tigris)
library(sf)
library(mgcv)
library(RColorBrewer)
library(rnaturalearth)
library(tidycensus)
library(units)

#### Load in Data ####

# List SSP data from cloud bucket
nc_files <- system("gsutil ls gs://clim_data_reg_useast1/misc_data/temporary//HVAC", intern = T)

# Identify files: AC penetration (ends in 2050.nc) and population (contains "pop")
ssp_files <- grep("2050.nc$", nc_files, value = T)

pop_file <- grep('pop', nc_files, value = T)

# Set up directories
destination_dir <- "/mnt/pers_disk/hvac"
# fs::dir_create(destination_dir) # Uncomment if directory doesn't exist

# Download the files to local disk 
walk(ssp_files, ~ system(paste("gsutil cp", .x, destination_dir)))
walk(pop_file, ~ system(paste("gsutil cp", .x, destination_dir)))

# List downloaded CDD NetCDF files
ac_files <- list.files("/mnt/pers_disk/hvac", pattern = "ac",
                       full.names = TRUE)

pop_files <- list.files("/mnt/pers_disk/hvac", pattern = "pop",
                         full.names = TRUE)

# AC penetration Data Handling

# Read AC data
ac_data <- read_ncdf(ac_files)
ac_df <- as.data.frame(ac_data, xy = TRUE)

ac_df$SSP <- as.vector(ac_df$SSP)  # for saving purposes

# Save cleaned AC penetration data
write_csv(ac_df, "/mnt/pers_disk/hvac_processed/ac_penetration.csv")

# Population Data Handling

# Read population data
pop_data <- read_ncdf(pop_files)
pop_df <- as.data.frame(pop_data, xy = TRUE)

# pop_df$SSP <- as.vector(pop_df$SSP)  # for saving purposes
# pop_df$population <- as.numeric(pop_df$population) # for saving purposes

# Save population data as CSV
write_csv(pop_df, "/mnt/pers_disk/hvac_processed/population.csv")





#### AC Penetration Analysis (SSP5-8.5) ####

# Read In AC penetration data
ac_df <- read_csv("/mnt/pers_disk/hvac_processed/ac_penetration.csv")

# Convert AC data to sf and filter for SSP5-8.5
ssp_585 <- ac_df %>%
  filter(SSP == 5) %>%
  mutate(decade = paste0(floor(Time / 10) * 10, "s")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Load country polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Spatial join with countries (for global summary)
ssp_585_joined <- st_join(ssp_585, world["admin"])

# Compute mean AC penetration by country and decade
ranked_585 <- ssp_585_joined %>%
  st_drop_geometry() %>%
  group_by(admin, decade) %>%
  summarise(mean_ac = mean(ac_penetration, na.rm = TRUE)) %>%
  ungroup()

# Rank countries by mean AC penetration within each decade
ranked_585 <- ranked_585 %>%
  group_by(decade) %>%
  arrange(desc(mean_ac)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Calculate total change in AC penetration between 2000s and 2040s
total_change <- ranked_585 %>%
  filter(decade %in% c("2000s", "2040s")) %>%
  select(admin, decade, mean_ac) %>%
  pivot_wider(names_from = decade, values_from = mean_ac) %>%
  rename(ac_2005 = `2000s`, ac_2045 = `2040s`) %>%
  mutate(total_increase = ac_2045 - ac_2005) %>%
  drop_na(total_increase) %>%
  arrange(desc(total_increase))


# Join ac change data to country geometries
change_map <- left_join(world, total_change, by = c("admin"))


# Plot: Increases in AC penetration by country (2000s–2040s)
ggplot(change_map) +
  geom_sf(aes(fill = total_increase), color = NA) +
  scale_fill_viridis_c(option = "inferno", na.value = "grey80", name = "AC Penetration\nIncrease (2000s–2040s)") +
  labs(
    title = "Countries with Fastest Growth in AC Penetration (SSP5-8.5)",
    caption = "Data: Nature Communications (2024)"
  ) +
  theme_minimal()

# Join ranked AC penetration data to country geometries (for faceted map)
map_data <- left_join(world, ranked_585, by = c("admin")) %>%
  filter(!is.na(mean_ac))

# Plot: AC penetration by country by decade
ggplot(map_data) +
  geom_sf(aes(fill = mean_ac), color = NA) +
  scale_fill_viridis_c(option = "inferno", na.value = "grey80", name = "AC Penetration") +
  labs(
    title = "AC Penetration by Country by Decade",
    subtitle = "SSP5-8.5 Scenario",
    caption = "Data: Nature Communications (2024)"
  ) +
  facet_wrap(~decade) +
  theme_minimal()






#### State-Level AC Penetration & Electricity Demand (SSP5-8.5) ####

# Load and Filter AC Penetration Data from above (ssp_585) 

# Load U.S. state boundaries
states <- states(cb = TRUE) %>%
  filter(STUSPS %in% state.abb) %>%
  st_transform(4326)  

# Align CRS and spatially join AC data with states
ssp_us_states <- ssp_585 %>%
  st_transform(st_crs(states)) %>%
  st_join(states["NAME"]) %>%
  drop_na(NAME)

# Aggregate mean AC penetration by state and decade
ac_by_state <- ssp_us_states %>%
  st_drop_geometry() %>%
  group_by(NAME, decade) %>%
  summarise(ac_pen = mean(ac_penetration, na.rm = TRUE), .groups = "drop") %>%
  rename(state = NAME) # Rename 'NAME' column to 'state' for consistency

# Join back with state polygons for plotting
ac_by_state <- inner_join(states, ac_by_state, by = c("NAME" = "state")) %>%
  st_as_sf()

# Plot: AC Penetration by State and Decade
ggplot(ac_by_state) +
  geom_sf(aes(fill = ac_pen), color = "white", size = 0.1) +
  scale_fill_viridis_c(name = "AC Penetration (%)", option = "magma", direction = -1) +
  facet_wrap(~decade) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +  # Zoom
  labs(
    title = "State-Level AC Penetration by Decade (SSP5-8.5)",
    caption = "Data: Nature Communications (2024)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


###### Load and Prepare Population Data #####

# Read in population data
pop_df <- read_csv("/mnt/pers_disk/hvac_processed/population.csv")

# Convert pop data to sf and filter for SSP5-8.5
pop_ssp5 <- pop_df %>%
  filter(SSP == 5) %>%  
  mutate(decade = paste0(floor(Time / 10) * 10, "s")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Join population data with state boundaries
pop_us_states <- pop_ssp5 %>%
  st_transform(st_crs(states)) %>%
  st_join(states["NAME"]) %>%
  drop_na(NAME)

# Aggregate total population by state and decade
pop_by_state <- pop_us_states %>%
  group_by(NAME, decade) %>%
  summarise(pop = sum(population, na.rm = TRUE), .groups = "drop") %>%
  rename(state = NAME) # Rename 'NAME' column to 'state' for consistency

# Join population and AC data
hvac_state_table <- left_join(ac_by_state %>% 
                                st_drop_geometry(), 
                              pop_by_state, 
                              by = c("NAME" = "state", 'decade')) %>%
  mutate(
    ac_users = (ac_pen / 100) * pop, # Calculate number of AC users
    demand_kwh = ac_users * 800, # Estimate demand in kWh (assuming 800 kWh per user)
    demand_gwh = demand_kwh / 1e6 # Convert to GWh
  )

# # Save data as CSV file (I saved with baseline and additional demand below)
# write_csv(hvac_state_table, "/mnt/pers_disk/hvac_processed/hvac_state_table.csv")

###### Calculate Additional Demand Relative to 2020s as baseline #####

# Extract baseline demand (2020s)
baseline <- hvac_state_table %>%
  filter(decade == "2020s") %>%
  select(NAME, base_demand = demand_gwh, base_pop = pop)

# Join and compute additional demand
hvac_state_additional <- hvac_state_table %>%
  left_join(baseline, by = 'NAME') %>%
  mutate(addl_demand_gwh = demand_gwh - base_demand,
         addl_pop = pop - base_pop)

# Subset to future decades
hvac_future <- hvac_state_additional %>%
  filter(decade %in% c("2030s", "2040s"))

# Plot: Bar Chart of Additional Demand by State and Decade
ggplot(hvac_future, aes(x = decade, y = addl_demand_gwh, fill = decade)) +
  geom_col() +
  facet_wrap(~ NAME) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Future Increase in Residential Electricity Demand for AC",
    subtitle = "Relative to 2020s baseline (SSP5-8.5)",
    x = "Decade",
    y = "Additional Demand (GWh)",
    caption = "Data: Nature Communications (2024)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# Plot: Heatmap of Additional Demand by State and Decade
ggplot(hvac_state_additional %>% 
         filter(decade %in% c("2030s", "2040s"))) +
  aes(x = decade, y = NAME, fill = addl_demand_gwh) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Additional Demand (GWh)", option = "plasma") +
  labs(
    title = "Projected Additional Electricity Demand for Residential AC",
    subtitle = "Relative to 2020s baseline (SSP5-8.5)",
    x = "Decade",
    y = "State"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

# Join additional demand data back to state polygons and filter for future decades
map_data <- left_join(states, hvac_state_additional, by = "NAME") %>%
  st_as_sf() %>%
  filter(decade %in% c("2030s", "2040s"))  # only future decades

# Plot: Map of Additional Demand by State
ggplot(map_data) +
  geom_sf(aes(fill = addl_demand_gwh), color = "white", size = 0.1) +
  scale_fill_viridis_c(
    name = "Additional Demand (GWh)",
    option = "inferno",  # More distinct than plasma/magma
    direction = -1
  ) +
  facet_wrap(~decade) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  labs(
    title = "Additional Electricity Demand for Residential AC (vs 2020s Baseline)",
    subtitle = "SSP5-8.5 Scenario",
    caption = "Data: Nature Communications (2024)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Rescale addl_pop by decade
map_data_scaled <- map_data %>%
  group_by(decade) %>%
  mutate(
    addl_pop_scaled = scales::rescale(addl_pop, to = c(-1, 1))  # Standardize within each decade
  ) %>%
  ungroup()

# Plot: Relative Population Change by State and Decade
ggplot(map_data_scaled) +
  geom_sf(aes(fill = addl_pop_scaled), color = "white", size = 0.1) +
  scale_fill_gradient2(
    name = "Relative Change",
    low = "blue", mid = "white", high = "red",
    midpoint = 0
  ) +
  facet_wrap(~decade) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  labs(
    title = "Relative Change in Population by State",
    subtitle = "Rescaled within each decade (SSP5-8.5)",
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Save data as CSV file
write_csv(hvac_state_additional, "/mnt/pers_disk/hvac_processed/hvac_state_table.csv")

##### Analysis for Texas ######

# Preview key metrics for Texas
hvac_state_additional %>%
  filter(NAME == "Texas") %>%
  select(decade, ac_pen, pop, ac_users, demand_gwh, addl_demand_gwh)

# Bar plot Texas AC Penetration by Decade
hvac_state_table %>%
  filter(NAME == "Texas") %>%
  ggplot(aes(x = decade, y = ac_pen)) +
  geom_col(fill = "#21918c") +
  labs(title = "Texas: AC Penetration by Decade", y = "AC Penetration (%)", x = NULL) +
  theme_minimal()

# Line plot Texas Population by Decade
hvac_state_table %>%
  filter(NAME == "Texas") %>%
  ggplot(aes(x = decade, y = pop / 1e6)) +
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  labs(title = "Texas: Population by Decade", y = "Population (millions)", x = NULL) +
  theme_minimal()

# Line plot Texas Total Electricity Demand by Decade
hvac_state_table %>%
  filter(NAME == "Texas") %>%
  ggplot(aes(x = decade, y = demand_gwh)) +
  geom_line(group = 1, color = "black") +
  geom_point(size = 2) +
  labs(title = "Texas: Total Electricity Demand (GWh)", y = "Demand (GWh)", x = NULL) +
  theme_minimal()



