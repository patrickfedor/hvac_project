#### Set up ####

library(tidyverse)
library(stars)
library(viridis)
library(rnaturalearth)  # Country boundary shapefiles
library(sf)
library(exactextractr)
library(terra)
library(future)         # For parallel processing

plan(multisession) #enable parallel processing


#### Process Data ####

# List CDD/HDD data from cloud bucket
nc_files <- system("gsutil ls gs://cmip6_data/hvac/", intern = T)

# Identify files
degreeday_files <- grep(".nc$", nc_files, value = T)

# Set up directories
destination_dir <- "/mnt/pers_disk/hvac"
# fs::dir_create(destination_dir) # Uncomment if directory doesn't exist

# Download the files to local disk 
walk(degreeday_files, ~ system(paste("gsutil cp", .x, destination_dir)))

# Get the list of .nc files in the directory
nc_files <- list.files("/mnt/pers_disk/hvac", pattern = "deg.nc$", full.names = TRUE)

# Read .nc file into stars objects
nc_list <- lapply(nc_files, function(f) {
  obj <- read_stars(f)
  attr(obj, "filename") <- f
  obj
})

# Load country boundary data as an sf object
countries <- ne_countries(scale = "medium", returnclass = "sf")

# Separate CDD and HDD datasets using filename pattern
cdd_list <- nc_list[grepl("cdd", sapply(nc_list, attr, "filename"))]
hdd_list <- nc_list[grepl("hdd", sapply(nc_list, attr, "filename"))]


# Function to extract warming level from filename 
# (e.g., "cdd_2.0deg.nc" → 2.0)
get_warming_level <- function(st) {
  filename <- attr(st, "filename")
  as.numeric(gsub(".*_(\\d+\\.\\d+)deg\\.nc", "\\1", filename))
}

# Identify and extract the baseline (1.0°C warming) files
baseline_cdd <- cdd_list[[which(sapply(cdd_list, get_warming_level) == 1.0)]]
baseline_hdd <- hdd_list[[which(sapply(hdd_list, get_warming_level) == 1.0)]]

#### Calculate Percent Change ####

# Function to compute percent change relative to baseline
percent_change <- function(future, baseline) {
  change <- (future - baseline) / baseline * 100
  return(change)
}

# Apply percent change function to each warming level file
cdd_change_list <- lapply(cdd_list, function(future) {
  change <- percent_change(future, baseline_cdd)
  attr(change, "filename") <- attr(future, "filename")  # Retain filename attribute
  change
  change
})

hdd_change_list <- lapply(hdd_list, function(future) {
  change <- percent_change(future, baseline_hdd)
  attr(change, "filename") <- attr(future, "filename")
  change
})


#### Fix Raster Metadata ####

# Function to correct raster dimensions and CRS 
fix_raster_extent <- function(raster) {
  x_coords <- seq(-180, 180, length.out = dim(raster)[1])
  y_coords <- seq(90, -90, length.out = dim(raster)[2])
  
  raster <- st_set_dimensions(raster, 1, values = x_coords, names = "x")
  raster <- st_set_dimensions(raster, 2, values = y_coords, names = "y")
  st_crs(raster) <- 4326
  return(raster)
}


# Apply fix to all percent change rasters
hdd_change_list <- lapply(hdd_change_list, fix_raster_extent)
cdd_change_list <- lapply(cdd_change_list, fix_raster_extent)

# Function to plot a stars object with country outlines
plot_stars <- function(stars_obj, title) {
  df <- as.data.frame(stars_obj, xy = TRUE)
  col_name <- names(df)[3]
  ggplot(df, aes(x = x, y = y, fill = !!sym(col_name))) +
    geom_raster() +
    scale_fill_viridis_c(option = "turbo", na.value = "white", limits = c(-100, 0),
                         oob = scales::squish) +
    labs(title = title, fill = "% Change") +
    geom_sf(data = countries, fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
    coord_sf(crs = st_crs(4326), expand = FALSE) +
    theme_minimal()
}


# Example: plot HDD % change at a specific warming level (e.g., 2°C)
plot_stars(cdd_change_list[[5]], paste("CDD % Change", get_warming_level(cdd_change_list[[5]]), "°C"))

#### Optional: Raw Difference Plot ####

# Compute and visualize raw difference in HDD values (not percent change)
raw_diff <- cdd_list[[3]] - baseline_hdd
raw_diff <- fix_raster_extent(raw_diff)

# Plot
plot_stars(raw_diff, "Difference CDD 2°C - 1.0°C")






#### Summarize Country-Level Mean Change at 3°C ####

# Extract the 3°C warming scenario for CDD
cdd_3deg <- cdd_list[[which(sapply(cdd_list, get_warming_level) == 3.0)]]
cdd_3deg <- fix_raster_extent(cdd_3deg)


# Convert stars object to terra SpatRaster for extraction
cdd_3deg <- st_as_stars(cdd_3deg) %>% 
  as("Raster") %>% 
  rast()

# Make sure countries layer is in same CRS
countries <- st_transform(countries, st_crs(cdd_3deg))

# Extract mean percent change per country
country_summary <- exact_extract(cdd_3deg, countries, 'mean', progress = TRUE)

# Combine with country names in df
summary_cdd_3deg <- data.frame(
  country = countries$name,
  mean_cdd_3deg = country_summary
) %>%
  arrange(country)


# Preview summary
print(head(summary_cdd_3deg))
print(tail(summary_cdd_3deg, n = 100))

# Save results to CSV
write_csv(summary_cdd_3deg, "country_cdd_3C.csv")

# Upload to cloud bucket 
system("gsutil cp /home/pfedor/country_cdd_3C.csv gs://clim_data_reg_useast1/misc_data/temporary/country_cdd_3C.csv")



#### Top countries ####

# Calculate mean and total CDD by country for the baseline
baseline_cdd_fixed <- fix_raster_extent(baseline_cdd)
baseline_cdd_rast <- st_as_stars(baseline_cdd_fixed) %>% as("Raster") %>% rast()

# Match countries to raster CRS
countries_proj <- st_transform(countries, st_crs(baseline_cdd_rast))


# Mean CDD by country
country_cdd_means <- exact_extract(baseline_cdd_rast, countries_proj, 'mean')

summary_cdd <- tibble(
  country = countries_proj$name,
  mean_cdd = country_cdd_means
)

# Sort and get top 15
top_cdd_countries <- summary_cdd %>%
  arrange(desc(mean_cdd)) %>%
  slice_head(n = 15)

# Print result
print(top_cdd_countries)



# Total CDD by country
summary_cdd_total <- tibble(
  country = countries_proj$name,
  total_cdd = total_cdd_by_country
)

# Sort and get top 15
top_total_cdd_countries <- summary_cdd_total %>%
  arrange(desc(total_cdd)) %>%
  slice_head(n = 15)

# View result
print(top_total_cdd_countries)

write_csv(top_cdd_countries, "top15_cdd_countries.csv")
