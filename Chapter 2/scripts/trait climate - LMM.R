# Required libraries

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(extrafont)

font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# weather info

spieke_spring <- read.csv("Desktop/git/Chapter 2/resulting data/spieke-spring RDA.csv", 
                          check.names = FALSE)
brach_spring <- read.csv("Desktop/git/Chapter 2/resulting data/brach-spring RDA.csv", 
                         check.names = FALSE)

# Traits dataset

spring_coll <- read.csv("Desktop/git/Chapter 2/data/spring_all.csv")

# Extracting the four months prior to collection

spieke_spring2 <- spieke_spring %>%
  mutate(grouping = case_when(
    year == 2020 | year == 2021 & month %in% c(1, 2, 3) ~ "A",
    year == 2021 & month %in% c(12) | year == 2022 & month %in% c(1, 2, 3) ~ "B", 
    year == 2022 & month %in% c(12) | year == 2023 & month %in% c(1, 2, 3) ~ "C", 
    year == 2023 & month %in% c(12) | year == 2024 & month %in% c(1, 2, 3) ~ "D",
    year == 2024 & month %in% c(12) | year == 2025 & month %in% c(1, 2, 3) ~ "E"
  )) %>%
  drop_na(grouping) %>%
  droplevels()

brach_spring2 <- brach_spring %>%
  mutate(grouping = case_when(
    year == 2022 & month %in% c(12) | year == 2023 & month %in% c(1, 2, 3) ~ "A",
    year == 2023 & month %in% c(12) | year == 2024 & month %in% c(1, 2, 3) ~ "B",
    year == 2024 & month %in% c(12) | year == 2025 & month %in% c(1, 2, 3) ~ "C"
  )) %>%
  drop_na(grouping) %>%
  droplevels()

# Weather means by groups

spring_spieke <- spieke_spring2 %>%
  select(c(3:12, 15)) %>%
  group_by(grouping) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

spring_brach <- brach_spring2 %>%
  select(c(3:12, 15)) %>%
  group_by(grouping) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

# Adding year

spring_spieke <- spring_spieke %>%
  mutate(year = c(2021, 2022, 2023, 2024, 2025))
spring_brach <- spring_brach %>%
  mutate(year = c(2023, 2024, 2025))

# Preparing trait df

spring_df <- spring_coll %>%
  filter(location == "Spiekeroog" | location == "Brachwitz") %>%
  droplevels()

spring_df <- spring_df %>%
  rename(
    `Length longest leaf` = length_longest_leaf_cm,
    `Petiole longest leaf` = petiole_longest_leaf_cm,
    `Stem width` = stem_cm,
    `Width longest leaf` = width_longest_leaf_cm, 
    `Leaf number` = leaves_number,
    `Cauline leaf` = Cauline_leaves,
    Flowers = flowers,
    `Plant length` = plant_length_cm,
    `Aspect ratio` = aspect_ratio,
    `Petiole length ratio` = petiole_length_ratio
  ) # Renaming some columns

# Subsetting traits by location and merging weather data

traits <- spring_df[, c(1:14)] %>%
  drop_na()

traits_brach <- traits %>%
  filter(location == "Brachwitz") %>%
  merge(., spring_brach, by = "year") %>%
  drop_na()
traits_brach$year <- as.factor(traits_brach$year)

traits_spieke <- traits %>%
  filter(location == "Spiekeroog") %>%
  merge(., spring_spieke, by = "year") %>%
  drop_na()
traits_spieke$year <- as.factor(traits_spieke$year)

# merging locations

all <- bind_rows(traits_spieke, traits_brach) %>%
  column_to_rownames("plant_ID") %>%
  droplevels() %>%
  drop_na()

# scaling and transforming 
all_long <- all %>%
  mutate_at(c(2:11, 15:24), ~scale(., center = TRUE, scale = TRUE)) %>%   #, 15:24
  pivot_longer(cols = 2:11, 
              names_to = "trait", 
              values_to = "trait_value")%>%
  mutate_at(c(1), as.factor) %>%
  droplevels()

all_long_spieke <- all_long %>%
  filter(location == "Spiekeroog") %>%
  droplevels()

all_long_brachwitz <- all_long %>%
  filter(location == "Brachwitz") %>%
  droplevels()

# lmm

corr <- cor(select(all_long, precipitation_sum, `temperature_2m_mean (°C)`, `wind_speed_10m_mean (km/h)`, 
                   `dew_point_2m_mean (°C)`, `relative_humidity_2m_mean (%)`, `sunshine_duration (s)`, 
                   `soil_temperature_0_to_7cm_mean (°C)`), use = "pairwise.complete.obs")

#all_lmm <- lmer(trait_value ~ (
#  precipitation_sum + `temperature_2m_mean (°C)` + `soil_temperature_0_to_7cm_mean (°C)` + 
#    `wind_speed_10m_mean (km/h)` + `relative_humidity_2m_mean (%)` + `dew_point_2m_mean (°C)` + 
#    `sunshine_duration (s)`) * trait + (1 | year), data = all_long)

all_lmm <- lmer(trait_value ~ (
  precipitation_sum + `temperature_2m_mean (°C)` + 
    `wind_speed_10m_mean (km/h)` + `sunshine_duration (s)`) * 
    trait + (1 | year), data = all_long)

summary(all_lmm)
anova(all_lmm)

sig_effects <- tidy(all_lmm, effects = "fixed") %>%
  filter(p.value < 0.05)

ranef_df <- as.data.frame(VarCorr(all_lmm))
ranef_df

write.csv(sig_effects, "Desktop/git/Chapter 2/resulting data/sig_effects.csv", 
          row.names = FALSE)

# Spieke
corr2 <- cor(select(all_long_spieke, precipitation_sum, `temperature_2m_mean (°C)`, `wind_speed_10m_mean (km/h)`, 
                   `dew_point_2m_mean (°C)`, `relative_humidity_2m_mean (%)`, `sunshine_duration (s)`, 
                   `soil_temperature_0_to_7cm_mean (°C)`), use = "pairwise.complete.obs")


all_lmm_spieke <- lmer(trait_value ~ (
  precipitation_sum + `temperature_2m_mean (°C)` + 
    `wind_speed_10m_mean (km/h)` + `sunshine_duration (s)`) * 
    trait + (1 | year), data = all_long_spieke)

summary(all_lmm_spieke)
anova(all_lmm_spieke)

sig_effects_spieke <- tidy(all_lmm_spieke, effects = "fixed") %>%
  filter(p.value < 0.05)

ranef_df1 <- as.data.frame(VarCorr(all_lmm_spieke))
ranef_df1

write.csv(sig_effects_spieke, "Desktop/git/Chapter 2/resulting data/sig_effects_spieke.csv", 
          row.names = FALSE)

# Brachwitz

all_lmm_brach <- lmer(trait_value ~ (
  precipitation_sum + `temperature_2m_mean (°C)` + 
    `wind_speed_10m_mean (km/h)` + `sunshine_duration (s)`) * 
    trait + (1 | year), data = all_long_brachwitz)

summary(all_lmm_brach)
anova(all_lmm_brach)

sig_effects_brach <- tidy(all_lmm_brach, effects = "fixed") %>%
  filter(p.value < 0.05)

ranef_df2 <- as.data.frame(VarCorr(all_lmm_brach))
ranef_df2

write.csv(sig_effects_brach, "Desktop/git/Chapter 2/resulting data/sig_effects_brach.csv", 
          row.names = FALSE)


# plotting 
cols = c("traitFlowers" = "#1b9e77", "traitLeaf number" = "#d95f02", 
         "traitLength longest leaf" = "#7570b3", "traitPetiole longest leaf" = "#e7298a", 
         "traitPlant length" = "#66a61e", "traitSidebranch" = "#e6ab02",
         "traitSideshoots" = "#a6761d", "traitStem width" = "#666666", 
         "traitWidth longest leaf" = "#8dd3c7")
label = c("traitFlowers" = "Flowers", "traitLeaf number" = "Leaf number", 
          "traitLength longest leaf" = "Length longest leaf", 
          "traitPetiole longest leaf" = "Petiole longest leaf", 
          "traitPlant length" = "Plant length", "traitSidebranch" = "Side branches",
          "traitSideshoots" = "Side shoots", "traitStem width" = "Stem width", 
          "traitWidth longest leaf" = "Width longest leaf")

lmm_plot <- tidy(all_lmm, effects = "fixed") %>%
  filter(str_detect(term, ":trait")) %>%
  filter(p.value < 0.05) %>%
  separate(term, into = c("climate", "trait"), sep = ":") %>%
  ggplot(aes(x = climate, y = estimate, fill = trait)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = cols, label = label, name = "Trait measurements") +
  labs(y = "Effect on trait value", title = "Weather-trait interaction effects") +
  theme_pub4()
lmm_plot  

ggsave(filename = "Desktop/git/Chapter 2/plots/Figure-2c.pdf", 
       plot = lmm_plot, width = 7, height = 5, units = "in", dpi = 500)




