# libraries

library(tidyverse)
library(readxl)
library(parzer)

source("Desktop/git/functions/lat_long.R")

# Loading trait measurements
winter_2021 <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", 
                         sheet = 2) 
spring_2021 <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", 
                         sheet = 1) 
spring_2022 <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", 
                         sheet = 4)
spring_2023 <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", 
                         sheet = 5)
spring_2024 <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", 
                         sheet = 6)
spring_2025 <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", 
                         sheet = 7)
location <- read_xlsx("Desktop/git/Chapter 1/data/Phenotype Measurements - Updated.xlsx", sheet = 3)

# Winter trait measurements were separated from the location info 
# thus first we combine them

winter_2021_2 <- merge(winter_2021, location, by = "plant_ID")

# Transforming location coordinates

## Winter
winter_2021_2$longitude_GMS <- sapply(winter_2021_2$longitude_GMS, convert_latitude_format) # making a prepare dms formated columns
winter_2021_2$latitude_GMS <- sapply(winter_2021_2$latitude_GMS, convert_latitude_format) # making a prepare dms formated columns

winter_2021_2$longitude_GMS <- parse_lon(winter_2021_2$longitude_GMS)
winter_2021_2$latitude_GMS <- parse_lat(winter_2021_2$latitude_GMS)

## Spring

spring_2025 <- spring_2025 %>%
  mutate(
    latitude_GMS = str_split_i(Coordinates, ", ", 1),
    longitude_GMS = str_split_i(Coordinates, ", ", 2)
  )

# Joining collections

spring_cols <- Reduce(
  intersect, list(
    colnames(spring_2021), colnames(spring_2022), 
    colnames(spring_2023), colnames(spring_2024), 
    colnames(spring_2025)
  )
) # checking and extracting common columns in spring data 

spring_2021_col <- spring_2021 %>%
  select(any_of(spring_cols))
spring_2022_col <- spring_2022 %>%
  select(any_of(spring_cols))
spring_2023_col <- spring_2023 %>%
  select(any_of(spring_cols))
spring_2024_col <- spring_2024 %>%
  select(any_of(spring_cols))
spring_2025_col <- spring_2025 %>%
  select(any_of(spring_cols))

spring_all <- rbind(
  spring_2021_col, spring_2022_col, spring_2023_col, 
  spring_2024_col, spring_2025_col
)

# Re-setting columns datatypes

spring_all <- spring_all %>%
  mutate_at(c("leaves_number", "petiole_longest_leaf_cm", "length_longest_leaf_cm", 
              "width_longest_leaf_cm", "stem_mm", "Cauline_leaves", "flowers", "plant_length_cm", 
              "Sideshoots", "Sidebranch", "longitude_GMS", "latitude_GMS", "temperature_C"), as.numeric) %>%
  mutate_at(c("year", "season", "location"), as.factor) %>%
  select(-comments) # for spring

winter_all <- winter_2021_2 %>%
  mutate_at(c("petiole_longest_leaf_cm", "length_longest_leaf_cm", "width_longest_leaf_cm", 
              "temperature_C", "leaves_number", "longitude_GMS", "latitude_GMS"), as.numeric) %>%
  mutate_at(c("year", "season", "location"), as.factor)  # for winter

# renaming stem

spring_all <- spring_all %>%
  rename(stem_cm = stem_mm)

# Generating extra columns (Trait ratios)

winter_all <- winter_all %>%
  mutate(
    petiole_length_ratio = petiole_longest_leaf_cm / length_longest_leaf_cm,
    aspect_ratio = width_longest_leaf_cm / length_longest_leaf_cm
  )

spring_all <- spring_all %>%
  mutate(petiole_length_ratio = petiole_longest_leaf_cm / length_longest_leaf_cm,
         aspect_ratio = width_longest_leaf_cm / length_longest_leaf_cm,
         leaf_to_stem = length_longest_leaf_cm / stem_cm,
         leaf_plant_length = leaves_number / plant_length_cm,
         stem_plant_length = stem_cm / plant_length_cm,
         Sideshoot_to_plant_length = Sideshoots / plant_length_cm,
         flower_to_plnat_length = flowers / plant_length_cm,
         plant_length_latitude = plant_length_cm / latitude_GMS
  )


# Joining Winter and Spring collections

winter_spring_col <- Reduce(
  intersect, list(colnames(winter_all), colnames(spring_all))
)

winter_same <- winter_all %>%
  select(any_of(winter_spring_col))
spring_same <- spring_all %>%
  select(any_of(winter_spring_col))

winter_spring_all <- rbind(winter_same, spring_same)

# Extracting only Spiekeroog and Brachwitz samples

winter_spring_all_SB <- winter_spring_all %>%
  filter(location == "Spiekeroog" | location == "Brachwitz") %>%
  droplevels()

# Saving resulting files

write.csv(winter_spring_all_SB, "Desktop/git/Chapter 1/data/winter_spring_all_SB.csv", 
          quote = FALSE, row.names = FALSE)
write.csv(spring_all, "Desktop/git/Chapter 1/data/spring_all.csv", 
          quote = FALSE, row.names = FALSE)
write.csv(winter_all, "Desktop/git/Chapter 1/data/winter_all.csv", 
          quote = FALSE, row.names = FALSE)
