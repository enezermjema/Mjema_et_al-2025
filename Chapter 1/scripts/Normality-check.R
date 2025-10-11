# Required libraries

library(tidyverse)
library(ggpubr)
library(car)

source("Desktop/git/functions/norm_var_tests.R")

# Winter and spring data
data_winter <- read.csv("Desktop/git/Chapter 1/data/winter_all.csv")
data_spring <- read.csv("Desktop/git/Chapter 1/data/spring_all.csv")

# Data formatting

## Winter

df_winter_spieke <- data_winter %>%
  filter(location == "Spiekeroog") %>%
  drop_na(length_longest_leaf_cm) %>%
  mutate_at(c("season", "year", "location"), as.factor) %>%
  mutate_at(c("aspect_ratio", "petiole_length_ratio"), as.numeric) %>%
  select(c("length_longest_leaf_cm", "petiole_longest_leaf_cm", 
           "width_longest_leaf_cm", "leaves_number", "season", "year", 
           "location", "aspect_ratio", "petiole_length_ratio")) %>%
  droplevels()

## Spring

data_spring_2 <- data_spring %>%
  filter(location == "Spiekeroog" | location == "Brachwitz") %>%
  drop_na(length_longest_leaf_cm) %>%
  mutate_at(c("season", "year", "location"), as.factor) %>%
  mutate_at(c("aspect_ratio", "petiole_length_ratio"), as.numeric) %>%
  select(c("length_longest_leaf_cm", "petiole_longest_leaf_cm", "stem_cm", 
           "width_longest_leaf_cm", "leaves_number", "Cauline_leaves", "flowers",
           "plant_length_cm", "Sideshoots", "Sidebranch", "season", "year", "location", 
           "aspect_ratio", "petiole_length_ratio")) %>%
  droplevels() %>%
  droplevels()

df_spring_spieke <- data_spring_2 %>%
  filter(location == "Spiekeroog") %>%
  droplevels()

df_spring_brachwitz <- data_spring_2 %>%
  filter(location == "Brachwitz") %>%
  droplevels()

# Renaming some columns (Just make them presentable)

df_winter_spieke2 <- df_winter_spieke %>%
  rename(
    `Length longest leaf` = length_longest_leaf_cm,
    `Petiole longest leaf` = petiole_longest_leaf_cm,
    `Width longest leaf` = width_longest_leaf_cm, 
    `Leaf number` = leaves_number
  )

df_spring_spieke2 <- df_spring_spieke %>%
  rename(
    `Length longest leaf` = length_longest_leaf_cm,
    `Petiole longest leaf` = petiole_longest_leaf_cm,
    `Stem width` = stem_cm,
    `Width longest leaf` = width_longest_leaf_cm, 
    `Leaf number` = leaves_number,
    `Cauline leaf` = Cauline_leaves,
    Flowers = flowers,
    `Plant length` = plant_length_cm
  )

df_spring_brachwitz2 <- df_spring_brachwitz %>%
  rename(
    `Length longest leaf` = length_longest_leaf_cm,
    `Petiole longest leaf` = petiole_longest_leaf_cm,
    `Stem width` = stem_cm,
    `Width longest leaf` = width_longest_leaf_cm, 
    `Leaf number` = leaves_number,
    `Cauline leaf` = Cauline_leaves,
    Flowers = flowers,
    `Plant length` = plant_length_cm
  )

# Shapiro test

## Winter

columns_winter <- names(df_winter_spieke2[, c(1:4, 8, 9)]) # extracting names
winter_spieke_shapiro <- shapiro_tests(
  data = df_winter_spieke2, 
  factor_column = "year", 
  variables = columns_winter) # shapiro test

winter_shapiro <- shapiro_to_df(winter_spieke_shapiro) # creating a dataframe

## Spring

columns_spring <- names(df_spring_spieke2[, c(1:10, 14, 15)])

spring_spieke_shapiro <- shapiro_tests(
  data = df_spring_spieke2, 
  factor_column = "year", 
  variables = columns_spring) # spiekeroog spring samples

spring_brach_shapiro <- shapiro_tests(
  data = df_spring_brachwitz2, 
  factor_column = "year", 
  variables = columns_spring) # brachwitz spring samples

springSpieke_shapiro <- shapiro_to_df(spring_spieke_shapiro)
springBrach_shapiro <- shapiro_to_df(spring_brach_shapiro)


# Saving resulting df

write.csv(
  winter_shapiro, 
  "Desktop/git/Chapter 1/resulting data/winter_shapiro.csv", quote = FALSE)
write.csv(
  springSpieke_shapiro, 
  "Desktop/git/Chapter 1/resulting data/springSpieke_shapiro.csv", quote = FALSE)
write.csv(
  springBrach_shapiro, 
  "Desktop/git/Chapter 1/resulting data/springBrach_shapiro.csv", quote = FALSE)
