library(tidyverse)

data_winter <- read.csv("Desktop/git/Chapter 1/data/winter_all.csv")
data_spring <- read.csv("Desktop/git/Chapter 1/data/spring_all.csv")

# 750 winter plants

data_spring %>%
  group_by(location, year) %>%
  summarise(n())


df_winter_spieke <- data_winter %>%
  filter(location == "Spiekeroog") %>%
  drop_na(length_longest_leaf_cm) %>%
  mutate_at(c("season", "year", "location"), as.factor) %>%
  mutate_at(c("aspect_ratio", "petiole_length_ratio"), as.numeric) %>%
  select(c("length_longest_leaf_cm", "petiole_longest_leaf_cm", 
           "width_longest_leaf_cm", "leaves_number", "season", "year", 
           "location", "aspect_ratio", "petiole_length_ratio")) %>%
  droplevels()

# Spring

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

data_spring_3 <- data_spring_2 %>%
  mutate(
    group = case_when(
      year == "2021" ~ "Spring 2021",
      year == "2022" ~ "Spring 2022",
      year == "2023" ~ "Spring 2023",
      year == "2024" ~ "Spring 2024",
      year == "2025" ~ "Spring 2025"
    )
  )

df_spring_spieke <- data_spring_3 %>%
  filter(location == "Spiekeroog") %>%
  droplevels()

df_spring_brachwitz <- data_spring_3 %>%
  filter(location == "Brachwitz") %>%
  droplevels()


df_spring_spieke$group <- factor(
  df_spring_spieke$group,
  levels = c("Spring 2021", "Spring 2022", "Spring 2023", 
             "Spring 2024", "Spring 2025")
)
df_spring_brachwitz$group <- factor(
  df_spring_brachwitz$group,
  levels = c("Spring 2023", "Spring 2024", "Spring 2025")
)

# trait data points

df_winter_spieke %>%
  select(where(is.numeric)) %>%
  summarise(trait_data_points = sum(!is.na(across(everything()))))

df_spring_spieke %>%
  select(where(is.numeric), year) %>%
  group_by(year) %>%
  summarise(sum(!is.na(across(everything()))))

df_spring_brachwitz %>%
  select(where(is.numeric), year) %>%
  group_by(year) %>%
  summarise(sum(!is.na(across(everything()))))


# Transcriptome

trans_seq <- read.csv("Desktop/git/ena submissions/ena-samples list.csv")

trans_seq %>%
  group_by(location, year) %>%
  summarise(n())


