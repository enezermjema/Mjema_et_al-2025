# Loading libraries

library(tidyverse)
library(vegan)
library(factoextra)
library(FactoMineR)
library(ggpubr)
library(patchwork)
library(showtext)

showtext_auto()   # automatically use showtext for all plots
# Register Arial
font_add("Arial", regular = "Arial.ttf")

source("Desktop/git/functions/custom_theme.R")

# Datasets

data_winter <- read.csv("Desktop/git/Chapter 1/data/winter_all.csv")
data_spring <- read.csv("Desktop/git/Chapter 1/data/spring_all.csv")

# data prep

## winter
df_winter_spieke <- data_winter %>%
  filter(location == "Spiekeroog") %>%
  drop_na(length_longest_leaf_cm) %>%
  mutate_at(c("season", "year", "location"), as.factor) %>%
  mutate_at(c("aspect_ratio", "petiole_length_ratio"), as.numeric) %>%
  select(c("length_longest_leaf_cm", "petiole_longest_leaf_cm", 
           "width_longest_leaf_cm", "leaves_number", "season", "year", 
           "location", "aspect_ratio", "petiole_length_ratio")) %>%
  droplevels() # 660 plants

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
  droplevels()   # 2371 plants

# Renaming 

df_winter_spieke2 <- df_winter_spieke %>%
  rename(
    `Length longest leaf` = length_longest_leaf_cm,
    `Petiole longest leaf` = petiole_longest_leaf_cm,
    `Width longest leaf` = width_longest_leaf_cm, 
    `Leaf number` = leaves_number,
    `Aspect ratio` = aspect_ratio,
    `Petiole length ratio` = petiole_length_ratio
  )

data_spring_3 <- data_spring_2 %>%
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
  )

# Joining seasons

winter_spring_col <- Reduce(
  intersect, list(colnames(df_winter_spieke2), 
                  colnames(data_spring_3)))

winter_same <- df_winter_spieke2 %>%
  select(any_of(winter_spring_col))
spring_same <- data_spring_3 %>%
  select(any_of(winter_spring_col))

winter_spring_all <- rbind(winter_same, spring_same) # 3031 plants

# FAMD analysis
## Exploring the effect of Season and year of collection

winterSpring_drop <- winter_spring_all %>%
  drop_na()

all_famd <- FAMD(winterSpring_drop[, -c(8, 9)], graph = FALSE)

#fviz_screeplot(all_famd)
#fviz_famd_var(all_famd, repel = TRUE)

## contribution plot
all_contrib_1 <- fviz_contrib(all_famd, "var", axes = 1) + 
  theme_pub3() + 
  theme(axis.text.x = element_text(hjust = 1), 
        axis.title.x = element_blank())
all_contrib_2 <- fviz_contrib(all_famd, "var", axes = 2) + 
  theme_pub3() + 
  theme(axis.text.x = element_text(hjust = 1), 
        axis.title.x = element_blank())
all_contrib <- all_contrib_1 + all_contrib_2

## plotting FAMD results

famd_coord <- all_famd$ind$coord %>%
  as.data.frame()   # extracting famd object as a df

winterSpring_drop2 <- winterSpring_drop %>%
  mutate(group = case_when(
    year == 2021 & season == "Winter" ~ "Winter 2021",
    year == 2021 & season == "Spring" ~ "Spring 2021",
    year == 2022 & season == "Spring" ~ "Spring 2022",
    year == 2023 & season == "Spring" ~ "Spring 2023",
    year == 2024 & season == "Spring" ~ "Spring 2024",
    year == 2025 & season == "Spring" ~ "Spring 2025"
  ))   # creating a new column for easy mapping on a plot
winterSpring_drop2$group <- as.factor(winterSpring_drop2$group)

famd_coord2 <- bind_cols(
  famd_coord, winterSpring_drop2[, c(5:7, 10)]
)

famd_main <- ggplot(famd_coord2, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group, shape = location),
             size = 1.5, alpha = 1) +
  scale_color_manual(values = c("Winter 2021" = "skyblue", "Spring 2021" = "#4D4D4D", "Spring 2022" = "#737373",
                                "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", "Spring 2025" = "#CCCCCC"), 
                     name = "Collection") +
  scale_shape_manual(values = c("Spiekeroog" = 16, "Brachwitz" = 17), name = "Location") +
  labs(x = "Dimension 1 (29.2%)",
       y = "Dimension 2 (17%)") +
  theme_pub1()


# PCA
## Analyzing continuous columns for each season respectively

## Winter
df_winter_spieke2 <- df_winter_spieke2 %>%
  drop_na()

pca_winter <- PCA(df_winter_spieke2[, c(1:4)], 
                  graph = FALSE, scale.unit = TRUE)

winter_plot <- fviz_pca_biplot(
  pca_winter, col.ind = "lightgrey", col.var = "black", 
  repel = TRUE, label = c("ind", "var"), 
  title = "PCA on 659 Winter plants", geom.ind = "point") +
  theme_pub1()

## Spring
data_spring_3 <- data_spring_3 %>%
  drop_na()

pca_spring <- PCA(data_spring_3[, c(1:10)], 
                  graph = FALSE, scale.unit = TRUE)

spring_loc <- fviz_pca_biplot(
  pca_spring, col.ind = as.factor(data_spring_3$location), 
  col.var = "black", repel = TRUE, 
  label = c("ind", "var"), title = "PCA on 2193 Spring plants", 
  geom.ind = "point") +
  theme_pub1()   # location 

spring_year <- fviz_pca_biplot(
  pca_spring, col.ind = as.factor(data_spring_3$year), 
  col.var = "black", repel = TRUE, 
  label = c("ind", "var"), title = "PCA on 2193 Spring plants", 
  geom.ind = "point") +
  theme_pub1()   # year


# saving plots

ggsave(filename = "Desktop/git/Chapter 1/plots/famd_biplot.pdf", 
       plot = famd_main, width = 6, height = 4, units = "in", dpi = 450)

ggsave(filename = "Desktop/git/Chapter 1/plots/famd_contribution.pdf", 
       plot = all_contrib, width = 8, height = 4, units = "in", dpi = 450)

ggsave(filename = "Desktop/git/Chapter 1/plots/spring-PCA.pdf", 
       plot = spring_loc / spring_year, width = 7, 
       height = 10, units = "in", dpi = 450)

ggsave(filename = "Desktop/git/Chapter 1/plots/winter-PCA.pdf", 
       plot = winter_plot, width = 6, height = 4, units = "in", dpi = 450)

