# Required libraries

library(tidyverse)
library(vegan)
library(ggrepel)
library(adespatial)
library(grid)
library(extrafont)

font_import(pattern = "Arial")
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

x_all <- scale(all[, c(15, 16, 19, 21:24)]) %>%
  data.frame() %>%
  round(., digits = 3)
y_all <- scale(all[, c(2:11)]) %>%
  data.frame(check.names = FALSE) %>%
  round(., digits = 3)

# forward selection
# Selecting significant variables

all_sig = forward.sel(Y = y_all, X = x_all, alpha = 0.01)

x_all_sig <- subset(x_all, select = all_sig$variables)

# RDA 

rda_all <- rda(y_all ~ ., data = x_all_sig, scale = FALSE)
rda_all

RsquareAdj(rda_all)
rda_res <- anova(rda_all, by = "terms", permutations = 1000) %>%
  data.frame() %>%
  drop_na() %>%
  mutate(`contrib(%)` = (Variance / rda_all$CCA$tot.chi) * 100) %>%
  select(-1) # calculating percentage contributions

rda_res <- rda_res %>%
  mutate(`Weather variable`  = c(
    "Dew point (°C)", "Precipitation (mm)", "Soil temperature (°C)",
    "Wind speed (km/h)", "Sunshine duration (s)", "Air temperature (°C)",
    "Relative humidity (%)")) %>%
  data.frame(check.names = FALSE, row.names = NULL) # renaming weather parameters

write.csv(rda_res, "Desktop/git/Chapter 2/resulting data/rda-all.csv")

# Plotting a rda biplot

# Extract site scores
sites <- as.data.frame(scores(rda_all, display = "sites", scaling = 2))
sites$SampleID <- rownames(sites)

# Extract environmental variable scores
env_vectors <- as.data.frame(scores(rda_all, display = "bp", scaling = 2))
env_vectors <- env_vectors %>%
  mutate(Variables  = c(
    "Dew point (°C)", "Precipitation (mm)", "Soil temperature (°C)",
    "Wind speed (km/h)", "Sunshine duration (s)", "Air temperature (°C)",
    "Relative humidity (%)"),
    Variable_id = c(1, 2, 3, 4, 5, 6, 7)) %>%
  data.frame(row.names = NULL)

# Traits
species <- as.data.frame(scores(rda_all, display = "species", scaling = 2))
species <- species %>%
  rownames_to_column(var = "Traits")

# Biplots
key_text <- paste(env_vectors$Variable_id, env_vectors$Variables, sep = " = ")

biplot <- ggplot() +
  # Species (traits as arrows + labels)
  geom_segment(data = species, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "grey") +
  geom_text(data = species, aes(x = RDA1, y = RDA2, label = Traits),
            vjust = -0.5, size = 3, color = "navyblue") +
  
  # Predictors (bp arrows + labels)
  geom_segment(data = env_vectors, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text(data = env_vectors, aes(x = RDA1, y = RDA2, label = Variable_id),
            vjust = 2, size = 3, color = "black") +
  xlim(c(min(species$RDA1, env_vectors$RDA1) - 0.5,
         max(species$RDA1, env_vectors$RDA1) + 1)) +  # push right border  
  labs(
    x = "RDA1", y = "RDA2", subtitle = paste0(
      "Explained variation: ", round(RsquareAdj(rda_all)$adj.r.squared * 100, 2), "%")) +
  theme_pub1() +
  annotation_custom(
    grob = textGrob(paste(key_text, collapse = "\n"),
                    x = unit(0.75, "npc"),   # push to right
                    y = unit(0.26, "npc"), 
                    just = "left", gp = gpar(col = "black", fontsize = 10))
  ) +
  theme_pub1()

ggsave(filename = "Desktop/git/Chapter 2/plots/Figure-2c.pdf", 
       plot = biplot, width = 10, height = 7, units = "in", dpi = 500)
