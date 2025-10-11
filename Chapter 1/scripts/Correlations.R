# Required libraries

library(tidyverse)
library(ggpubr)
library(ggstatsplot)
library(ggcorrplot)
library(patchwork)
library(showtext)

showtext_auto()   # automatically use showtext for all plots
# Register Arial
font_add("Arial", regular = "Arial.ttf")

source("Desktop/git/functions/custom_theme.R")

# loading raw data

data_winter <- read.csv("Desktop/git/Chapter 1/data/winter_all.csv")
data_spring <- read.csv("Desktop/git/Chapter 1/data/spring_all.csv")

# Data prep

## winter
df_winter_spieke <- data_winter %>%
  filter(location == "Spiekeroog") %>%
  drop_na(length_longest_leaf_cm) %>%
  mutate_at(c("season", "year", "location"), as.factor) %>%
  select(c("length_longest_leaf_cm", "petiole_longest_leaf_cm", 
           "width_longest_leaf_cm", "leaves_number", "season", 
           "year", "location", "longitude_GMS", "latitude_GMS")) %>%
  droplevels()

## Spring
data_spring_2 <- data_spring %>%
  filter(location == "Spiekeroog" | location == "Brachwitz") %>%
  drop_na(length_longest_leaf_cm) %>%
  filter(plant_ID != "H047") %>%
  mutate_at(c("season", "year", "location"), as.factor) %>%
  select(c("length_longest_leaf_cm", "petiole_longest_leaf_cm", "stem_cm", 
           "width_longest_leaf_cm", "leaves_number", "Cauline_leaves", 
           "flowers", "plant_length_cm", "Sideshoots", "Sidebranch", "season", "year", 
           "location", "longitude_GMS", "latitude_GMS")) %>%
  droplevels()

df_spring_spieke <- data_spring_2 %>%
  filter(location == "Spiekeroog") %>%
  droplevels()

df_spring_brachwitz <- data_spring_2 %>%
  filter(location == "Brachwitz") %>%
  droplevels()

# Renaming columns

df_winter_spieke2 <- df_winter_spieke %>%
  rename(
    `Length longest leaf` = length_longest_leaf_cm,
    `Petiole longest leaf` = petiole_longest_leaf_cm,
    `Width longest leaf` = width_longest_leaf_cm, 
    `Leaf number` = leaves_number,
    Latitude = latitude_GMS
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
    `Plant length` = plant_length_cm,
    Latitude = latitude_GMS    
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
    `Plant length` = plant_length_cm,
    Latitude = latitude_GMS    
  )

# trait correlations for each location and across years, respectively

## winter spiekeroog

cor_winter <- cor(df_winter_spieke2[, c(1:4, 9)], 
                  method = "spearman", use = "pairwise.complete.obs")
## spring spiekeroog
cor_spieke_2021 <- cor(df_spring_spieke2[df_spring_spieke2$year == 2021, c(1:10)], 
                       method = "spearman", use = "pairwise.complete.obs")
cor_spieke_2022 <- cor(df_spring_spieke2[df_spring_spieke2$year == 2022, c(1:10, 15)], 
                       method = "spearman", use = "pairwise.complete.obs")
cor_spieke_2023 <- cor(df_spring_spieke2[df_spring_spieke2$year == 2023, c(1:10, 15)], 
                       method = "spearman", use = "pairwise.complete.obs")
cor_spieke_2024 <- cor(df_spring_spieke2[df_spring_spieke2$year == 2024, c(1:10, 15)], 
                       method = "spearman", use = "pairwise.complete.obs")
cor_spieke_2025 <- cor(df_spring_spieke2[df_spring_spieke2$year == 2025, c(1:10, 15)], 
                       method = "spearman", use = "pairwise.complete.obs")

## spring brachwitz
cor_brachwitz_2023 <- cor(df_spring_brachwitz2[df_spring_brachwitz2$year == 2023, c(1:10, 15)], 
                          method = "spearman", use = "pairwise.complete.obs")
cor_brachwitz_2024 <- cor(df_spring_brachwitz2[df_spring_brachwitz2$year == 2024, c(1:10, 15)], 
                          method = "spearman", use = "pairwise.complete.obs")
cor_brachwitz_2025 <- cor(df_spring_brachwitz2[df_spring_brachwitz2$year == 2025, c(1:10, 15)], 
                          method = "spearman", use = "pairwise.complete.obs")

# Computing pvalues

# spiekeroog samples
p_winter <- cor_pmat(df_winter_spieke2[, c(1:4, 9)], 
                     method = "spearman", )

p_spieke_2021 <- cor_pmat(df_spring_spieke2[df_spring_spieke2$year == 2021, c(1:10)], 
                          method = "spearman")
p_spieke_2022 <- cor_pmat(df_spring_spieke2[df_spring_spieke2$year == 2022, c(1:10, 15)], 
                          method = "spearman")
p_spieke_2023 <- cor_pmat(df_spring_spieke2[df_spring_spieke2$year == 2023, c(1:10, 15)], 
                          method = "spearman")
p_spieke_2024 <- cor_pmat(df_spring_spieke2[df_spring_spieke2$year == 2024, c(1:10, 15)], 
                          method = "spearman")
p_spieke_2025 <- cor_pmat(df_spring_spieke2[df_spring_spieke2$year == 2025, c(1:10, 15)], 
                          method = "spearman")

# spring brachwitz
p_brachwitz_2023 <- cor_pmat(df_spring_brachwitz2[df_spring_brachwitz2$year == 2023, c(1:10, 15)], 
                             method = "spearman")
p_brachwitz_2024 <- cor_pmat(df_spring_brachwitz2[df_spring_brachwitz2$year == 2024, c(1:10, 15)], 
                             method = "spearman")
p_brachwitz_2025 <- cor_pmat(df_spring_brachwitz2[df_spring_brachwitz2$year == 2025, c(1:10, 15)], 
                             method = "spearman")

# Correlation heatmaps

p1 <- ggcorrplot(
  cor_winter, 
  type = "lower", 
  p.mat = p_winter, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Winter 2021"
)

p2 <- ggcorrplot(
  cor_spieke_2021, 
  type = "lower", 
  p.mat = p_spieke_2021, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2021"
)

p3 <- ggcorrplot(
  cor_spieke_2022, 
  type = "lower", 
  p.mat = p_spieke_2022, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2022"
)

p4 <- ggcorrplot(
  cor_spieke_2023, 
  type = "lower", 
  p.mat = p_spieke_2023, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2023"
)

p5 <- ggcorrplot(
  cor_spieke_2024, 
  type = "lower", 
  p.mat = p_spieke_2024, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2024"
)

p6 <- ggcorrplot(
  cor_spieke_2025, 
  type = "lower", 
  p.mat = p_spieke_2025, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2025"
)

p7 <- ggcorrplot(
  cor_brachwitz_2023, 
  type = "lower", 
  p.mat = p_brachwitz_2023, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2023"
)

p8 <- ggcorrplot(
  cor_brachwitz_2024, 
  type = "lower", 
  p.mat = p_brachwitz_2024, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  show.legend = FALSE,
  title = "Spring 2024"
)
p9 <- ggcorrplot(
  cor_brachwitz_2025, 
  type = "lower", 
  p.mat = p_brachwitz_2025, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#6D9EC1", "white", "#E46726"),
  ggtheme = theme_pub1(), hc.order = FALSE,
  title = "Spring 2025"
)

# Combining some plots
spieke_all <- p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(nrow = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "a", title = "Spiekeroog")

brach_all <- p7 + p8 + 
  plot_layout(nrow = 1, byrow = TRUE) +
  plot_annotation(tag_levels = "a", title = "Brachwitz")

# Saving 

ggsave(filename = "Desktop/git/Chapter 1/plots/correlation-main(brachwitz2025).pdf", 
       plot = p9, width = 7, height = 7, units = "in", dpi = 450)
ggsave(filename = "Desktop/git/Chapter 1/plots/correlation-spiekeroog.pdf", 
       plot = spieke_all, width = 18, height = 11, units = "in", dpi = 450)
ggsave(filename = "Desktop/git/Chapter 1/plots/correlation-brachwitz.pdf", 
       plot = brach_all, width = 15, height = 11, units = "in", dpi = 450)
