# Required libraries

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)

# weather data

spieke <- read.csv("Desktop/git/Chapter 2/resulting data/Anomalies spiekeroog.csv")
brach <- read.csv("Desktop/git/Chapter 2/resulting data/Anomalies brachwitz.csv")

# Traits

spring_coll <- read.csv("Desktop/git/Chapter 2/data/spring_all.csv")

# location wise

spieke_traits <- spring_coll %>%
  filter(location == "Spiekeroog") %>%
  droplevels()
brach_traits <- spring_coll %>%
  filter(location == "Brachwitz") %>%
  droplevels()

# extracting petiole length

spieke_pl <- spieke_traits %>%
  select(c(3, 14))
brach_pl <- brach_traits %>%
  select(c(3, 14))

# Annual mean

spieke_clim <- spieke  %>%
  mutate(grouping = case_when(
    year == 2020 & month %in% c(10, 11, 12) | year == 2021 & month %in% c(1, 2, 3, 4, 5) ~ "A",
    year == 2021 & month %in% c(10, 11, 12) | year == 2022 & month %in% c(1, 2, 3, 4, 5) ~ "B", 
    year == 2022 & month %in% c(10, 11, 12) | year == 2023 & month %in% c(1, 2, 3, 4, 5) ~ "C", 
    year == 2023 & month %in% c(10, 11, 12) | year == 2024 & month %in% c(1, 2, 3, 4, 5) ~ "D",
    year == 2024 & month %in% c(10, 11, 12) | year == 2025 & month %in% c(1, 2, 3, 4, 5) ~ "E"
  )) %>%
  drop_na(grouping) %>%
  droplevels()
spieke_clim <- spieke  %>%
  mutate(grouping = case_when(
    year == 2021 & month %in% c(1, 2, 3, 4, 5) ~ "A",
    year == 2022 & month %in% c(1, 2, 3, 4, 5) ~ "B", 
    year == 2023 & month %in% c(1, 2, 3, 4, 5) ~ "C", 
    year == 2024 & month %in% c(1, 2, 3, 4, 5) ~ "D",
    year == 2025 & month %in% c(1, 2, 3, 4, 5) ~ "E"
  )) %>%
  drop_na(grouping) %>%
  droplevels()

brach_clim <- brach  %>%
  mutate(grouping = case_when(
    year == 2022 & month %in% c(10, 11, 12) | year == 2023 & month %in% c(1, 2, 3, 4, 5) ~ "A",
    year == 2023 & month %in% c(10, 11, 12) | year == 2024 & month %in% c(1, 2, 3, 4, 5) ~ "B",
    year == 2024 & month %in% c(10, 11, 12) | year == 2025 & month %in% c(1, 2, 3, 4, 5) ~ "C"
  )) %>%
  drop_na(grouping) %>%
  droplevels()

# Weather means by groups

spring_spieke <- spieke_clim %>%
  select(c(3:9, 24)) %>% #18
  group_by(grouping) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

spring_brach <- brach_clim %>%
  select(c(3:9, 24)) %>%
  group_by(grouping) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

# Adding year

spring_spieke <- spring_spieke %>%
  mutate(year = c(2021, 2022, 2023, 2024, 2025))
spring_brach <- spring_brach %>%
  mutate(year = c(2023, 2024, 2025))

# Merging trait and weather

traits_brach <- brach_pl %>%
  merge(., spring_brach, by = "year") %>% #spring_brach
  drop_na()
traits_brach$year <- as.factor(traits_brach$year)

traits_spieke <- spieke_pl %>%
  merge(., spring_spieke, by = "year") %>% # spring_spieke  spieke_clim
  drop_na()
traits_spieke$year <- as.factor(traits_spieke$year)


traits_spieke <- traits_spieke %>%
  mutate_at(c(2, 4:10), ~scale(.)) %>%
  mutate_at(c(4:10), as.factor) %>%
  droplevels()

# model for spiekeroog collection only

m1 <- lmer(petiole_longest_leaf_cm ~ 1 + (1 | year), 
           data = traits_spieke)

m2 <- lmer(petiole_longest_leaf_cm ~ temperature_2m_mean...C. + (1 | year), 
     data = traits_spieke)

summary(m1)
summary(m2)

model_comparison <- anova(m1, m2) %>%
  data.frame() %>%
  mutate(Conditional_R2 = c(MuMIn::r.squaredGLMM(m1)[2], MuMIn::r.squaredGLMM(m2)[2]))

write.csv(model_comparison, "Desktop/git/Chapter 2/resulting data/model_comparison.csv", 
          row.names = FALSE)
