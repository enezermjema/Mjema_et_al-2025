# Loading libraries

library(tidyverse)
library(ggpubr)
library(ggsignif)
library(readxl)
library(showtext)
library(patchwork)

showtext_auto()   # automatically use showtext for all plots
# Register Arial
font_add("Arial", regular = "Arial.ttf")

source("Desktop/git/functions/custom_theme.R")


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

# Data prep (setting y/n for 2024 and 2025 etc)

spring_2024 <- spring_2024 %>%
  mutate(Mainshoot = if_else(Mainshoot, "y", "n"))

spring_2025 <- spring_2025 %>%
  mutate(Mainshoot = case_when(Mainshoot == 1 ~ "y", 
                               Mainshoot == 0 ~ "n"))

spring_2021 <- spring_2021 %>%
  rename("Mainshoot" = `Mainshoot_y/n`)

# Spring collections

spring_col <- Reduce(
  intersect, list(colnames(spring_2021), 
                  colnames(spring_2022), colnames(spring_2023), 
                  colnames(spring_2024), colnames(spring_2025)))

spring_2345 <- Reduce(intersect, 
                      list(colnames(spring_2022), colnames(spring_2023), 
                           colnames(spring_2024), colnames(spring_2025)))

spring_2021_col <- spring_2021 %>%
  select(any_of(spring_col)) %>%
  select(c(plant_ID, Mainshoot, Sideshoots, 
           Sidebranch, Mainshoot, flowers, location, year))
spring_2022_col <- spring_2022 %>%
  select(any_of(spring_2345)) %>%
  select(c(plant_ID, Mainshoot, siliques_number, Sideshoots, 
           Sidebranch, Mainshoot, flowers, location, year))
spring_2023_col <- spring_2023 %>%
  select(any_of(spring_2345))  %>%
  select(c(plant_ID, Mainshoot, siliques_number, Sideshoots, 
           Sidebranch, Mainshoot, flowers, location, year))
spring_2024_col <- spring_2024 %>%
  select(any_of(spring_2345))  %>%
  select(c(plant_ID, Mainshoot, siliques_number, Sideshoots, 
           Sidebranch, Mainshoot, flowers, location, year))
spring_2025_col <- spring_2025 %>%
  select(any_of(spring_2345))  %>%
  select(c(plant_ID, Mainshoot, siliques_number, Sideshoots, 
           Sidebranch, Mainshoot, flowers, location, year))

spring_all <- rbind(spring_2022_col, spring_2023_col, spring_2024_col, spring_2025_col)

# re-setting columns datatypes

spring_all <- spring_all %>%
  mutate_at(c("flowers", "Sidebranch", 
              "Sideshoots", "siliques_number"), as.numeric) %>%
  mutate_at(c("Mainshoot", "year", "location"), as.factor) # for spring

spring_all_NA <- spring_all %>%
  filter(location %in% c("Spiekeroog", "Brachwitz")) %>%
  dplyr::mutate(reproductive_succ = flowers + siliques_number) %>%
  drop_na() %>%
  droplevels()

spring_2021_col <- spring_2021_col %>%
  drop_na()

# Renaming columns

spring_all_NA <- spring_all_NA %>%
  rename(
    Siliques = siliques_number,
    Flowers = flowers,
    `Reproductive success` = reproductive_succ
  )
spring_2021_col <- spring_2021_col %>%
  rename(
    Flowers = flowers
  )

# preping for boxplot

long_data <- spring_all_NA %>%
  pivot_longer(cols = c(Flowers, Siliques, `Reproductive success`, 
                        Sideshoots, Sidebranch), 
               names_to = "Variable", 
               values_to = "Value")

# individual boxes

sideshoots <- ggplot(spring_all_NA) +
  aes(x = Mainshoot, y = Sideshoots) +
  geom_boxplot(aes(fill = Mainshoot), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.05, colour = "black") +
  coord_cartesian(ylim = c(0, 12)) +
  scale_fill_manual(values = c("y" = "#4367B5", "n" = "#737373"), name = "Main shoot") +
  labs(x = "Mainshoot (y/n)", y = "Side shoots") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09, y_position = 10,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()

sidebranch <- ggplot(spring_all_NA) +
  aes(x = Mainshoot, y = Sidebranch) +
  geom_boxplot(aes(fill = Mainshoot), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.05, colour = "black") +
  coord_cartesian(ylim = c(0, 10)) +
  scale_fill_manual(values = c("y" = "#4367B5", "n" = "#737373"), name = "Main shoot") +
  labs(x = "Mainshoot (y/n)", y = "Side branch") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09, y_position = c(8),
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()

flowers <- ggplot(spring_all_NA) +
  aes(x = Mainshoot, y = Flowers) +
  geom_boxplot(aes(fill = Mainshoot), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.05, colour = "black") +
  coord_cartesian(ylim = c(0, 52)) +
  scale_fill_manual(values = c("y" = "#4367B5", "n" = "#737373"), name = "Main shoot") +
  labs(x = "Mainshoot (y/n)", y = "Flowers") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09, y_position = c(45),
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()

silique <- ggplot(spring_all_NA) +
  aes(x = Mainshoot, y = Siliques) +
  geom_boxplot(aes(fill = Mainshoot), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.05, colour = "black") +
  coord_cartesian(ylim = c(0, 60)) +
  scale_fill_manual(values = c("y" = "#4367B5", "n" = "#737373"), name = "Main shoot") +
  labs(x = "Mainshoot (y/n)", y = "Siliques") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09, y_position = c(50),
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()

reproductive <- ggplot(spring_all_NA) +
  aes(x = Mainshoot, y = `Reproductive success`) +
  geom_boxplot(aes(fill = Mainshoot), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.05, colour = "black") +
  coord_cartesian(ylim = c(0, 70)) +
  scale_fill_manual(values = c("y" = "#4367B5", "n" = "#737373"), name = "Main shoot") +
  labs(x = "Mainshoot (y/n)", y = "Reproductive success",
       caption = "Sample sizes: y = 1880, n = 128") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09, y_position = c(60),
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub1()

# saving plots

mainshoot_main <- sideshoots + flowers +silique
mainshoot_supple <- sidebranch + reproductive

ggsave(filename = "Desktop/git/Chapter 1/plots/mainshoot_main.pdf", 
       plot = mainshoot_main, width = 8, height = 4, units = "in", dpi = 450)
ggsave(filename = "Desktop/git/Chapter 1/plots/mainshoot_supple.pdf", 
       plot = mainshoot_supple, width = 7, height = 4, units = "in", dpi = 450)
