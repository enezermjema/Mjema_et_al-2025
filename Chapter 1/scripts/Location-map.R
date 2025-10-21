# Required libraries

library(tidyverse)
library(sf)
library(rnaturalearth)
library(leaflet)
library(ggmap)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# Loading dataset
df_all <- read.csv("Desktop/git/Chapter 1/data/winter_spring_all_SB.csv")

# Declaring factor column
df_all$season <- as.factor(df_all$season)
df_all$year <- as.factor(df_all$year)
df_all$location <- as.factor(df_all$location)

# Dropping NAs in coordinate columns
df_all <- df_all %>%
  drop_na(longitude_GMS) %>%
  drop_na(latitude_GMS)

# Creating map object

europe_sf <- ne_countries(scale = "medium", returnclass = "sf", 
                          continent = "Europe") # extracting Europe polygons info

germany_sf <- europe_sf %>%
  filter(sovereignt == "Germany") # Filtering for our country of interest

sample_sf <- st_as_sf(df_all, coords = c("longitude_GMS", "latitude_GMS"), crs = 4326)

# Plot Germany's boundary and sample collection points

sample_sf <- sample_sf %>%
  group_by(location) %>%
  mutate(sample_size = n())

# final map plot

final_map <- ggplot() +
  geom_sf(data = europe_sf, fill = "gray90", color = "black") +
  geom_sf(data = germany_sf, fill = "lightyellow", color = "black") +
  geom_sf(data = sample_sf, aes(colour = location)) +
  scale_size_continuous(range = c(2, 10)) +  # adjust point size range
  labs(x = "Longitude (°E)", y = "Latitude (°E)") +
  # Zooming into Germany by adjusting coordinate limits (bounding box)
  coord_sf(
    xlim = c(5, 16),  # Longitude range for Germany and a bit of its surroundings
    ylim = c(47, 55)  # Latitude range for Germany
  ) +
  theme_pub1()

# saving

ggsave(filename = "Desktop/git/Chapter 1/plots/collection_map.pdf", 
       plot = final_map, width = 6, height = 4, units = "in", dpi = 450)
