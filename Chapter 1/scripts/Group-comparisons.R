# Loading libraries

library(tidyverse)
library(ggpubr)
library(ggstatsplot)
library(ggcorrplot)
library(patchwork)
library(ggsignif)
library(showtext)

showtext_auto()   # automatically use showtext for all plots
# Register Arial
font_add("Arial", regular = "Arial.ttf")

source("Desktop/git/functions/plot_comparison.R")
source("Desktop/git/functions/custom_theme.R")

# loading raw data

data_winter <- read.csv("Desktop/git/Chapter 1/data/winter_all.csv")
data_spring <- read.csv("Desktop/git/Chapter 1/data/spring_all.csv")

# winter

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

# combining winter and spring

winter_spring <- rbind(df_winter_spieke[, c(1:9)], 
                       data_spring_2[, c(1:2, 4:5, 11:15)]) %>%
  data.frame() %>%
  drop_na() %>%
  droplevels()

winter_spring <- winter_spring %>%
  filter(petiole_longest_leaf_cm > 0)  # 2877 plants

winter_spring <- winter_spring %>%
  mutate(
    group = case_when(
      season == "Winter" ~ "Winter 2021",
      season == "Spring" & year == "2021" ~ "Spring 2021",
      year == "2022" ~ "Spring 2022",
      year == "2023" ~ "Spring 2023",
      year == "2024" ~ "Spring 2024",
      year == "2025" ~ "Spring 2025"
    ))

# Seasonally comparison (Winter (2021) vs Spring(2022))
# for leaf related traits

winter_spring2122 <- winter_spring %>%
  filter(group %in% c("Winter 2021", "Spring 2022")) %>%
  droplevels()
winter_spring2122$group <- factor(
  winter_spring2122$group,
  levels = c("Winter 2021", "Spring 2022")
)

ws1 <- ggplot(winter_spring2122) +
  aes(x = group, y = length_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Winter 2021" = "#4367B5", "Spring 2022" = "#737373"), 
                    name = "Collection") +
  labs(x = NULL, y = "Length longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()
ws2 <- ggplot(winter_spring2122) +
  aes(x = group, y = petiole_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Winter 2021" = "#4367B5", "Spring 2022" = "#737373"), 
                    name = "Collection") +
  labs(x = NULL, y = "Petiole longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()
ws3 <- ggplot(winter_spring2122) +
  aes(x = group, y = width_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Winter 2021" = "#4367B5", "Spring 2022" = "#737373"), 
                    name = "Collection") +
  labs(x = NULL, y = "Width longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub2()
ws4 <- ggplot(winter_spring2122) +
  aes(x = group, y = leaves_number) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Winter 2021" = "#4367B5", "Spring 2022" = "#737373"), 
                    name = "Collection") +
  labs(x = NULL, y = "Leaf number") +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  ) +
  theme_pub1()
ws <- ws1 + ws2 + ws3 + ws4 +
  plot_layout(nrow = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "a", 
                  title = "Winter vs Spring comparisons (2021 - 2022)")

ggsave(filename = "Desktop/git/Chapter 1/plots/leafTraits-winterVSspring.pdf", 
       plot = ws, width = 9, height = 9, units = "in", dpi = 450)

# Spring (Spiekeroog) across years (2021 - 2025)
ss1 <- ggplot(df_spring_spieke) +
  aes(x = group, y = length_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Length longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss2 <- ggplot(df_spring_spieke) +
  aes(x = group, y = petiole_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Petiole longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss3 <- ggplot(df_spring_spieke) +
  aes(x = group, y = width_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Width longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss4 <- ggplot(df_spring_spieke) +
  aes(x = group, y = leaves_number) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Leaf number") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss5 <- ggplot(df_spring_spieke) +
  aes(x = group, y = Cauline_leaves) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Caulie leaf number") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss6 <- ggplot(df_spring_spieke) +
  aes(x = group, y = stem_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Stem width") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss7 <- ggplot(df_spring_spieke) +
  aes(x = group, y = flowers) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Flower number") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss8 <- ggplot(df_spring_spieke) +
  aes(x = group, y = plant_length_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "#333333", "Spring 2022" = "#737373",
                               "Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Plant length") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
ss <- ss1 + ss2 + ss3 + ss4 + ss5 + ss6 + ss7 + ss8 +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "a", title = "Spiekeroog spring comparisons (2021 - 2025)")

ggsave(filename = "Desktop/git/Chapter 1/plots/allTraits - springSpiekeroog.pdf", 
       plot = ss, width = 19.69291, height = 12, units = "in", dpi = 450)

# Petiole length ratio (For weather results)
petioleRatio <- ggplot(df_spring_spieke) +
  aes(x = group, y = petiole_length_ratio) +
  geom_boxplot(aes(fill = group, color = group), outlier.shape = NA, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2.2)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.03, colour = "black") +
  scale_fill_manual(values = c("Spring 2021" = "skyblue", "Spring 2022" = "blue",
                               "Spring 2023" = "green", "Spring 2024" = "orange", 
                               "Spring 2025" = "red"), name = "Collection") +
  scale_color_manual(values = c("Spring 2021" = "skyblue", "Spring 2022" = "blue",
                                "Spring 2023" = "green", "Spring 2024" = "orange", 
                                "Spring 2025" = "red"), name = "Collection") +
  labs(x = NULL, y = "Petiole length ratio") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), 
                       c(2, 5), c(3, 4), c(3, 5), c(4, 5)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.04, y_position = c(0.8),
    textsize = 4, tip_length = 0.005, margin_top = 0.01
  )  +
  theme_pub2()
ggsave(filename = "Desktop/git/Chapter 1/plots/petioleRatio - spiekeroog.pdf", 
       plot = petioleRatio, width = 7, height = 6, units = "in", dpi = 450)

# Spring (Brachwitz) across years (2023 - 2025)
sb1 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = length_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Length longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb2 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = petiole_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Petiole longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb3 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = width_longest_leaf_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Width longest leaf") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb4 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = leaves_number) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Leaf number") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb5 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = Cauline_leaves) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Caulie leaf number") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb6 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = stem_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Stem width") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb7 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = flowers) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Flower number") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb8 <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = plant_length_cm) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "#989898", "Spring 2024" = "#B4B4B4", 
                               "Spring 2025" = "#CCCCCC"), name = "Collection") +
  labs(x = NULL, y = "Plant length") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.09,
    textsize = 5, tip_length = 0.01, margin_top = 0.01
  )  +
  theme_pub2()
sb <- sb1 + sb2 + sb3 + sb4 + sb5 + sb6 + sb7 + sb8 +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "a", title = "Brachwitz spring comparisons (2023 - 2025)")

ggsave(filename = "Desktop/git/Chapter 1/plots/allTraits - springBrachwitz.pdf", 
       plot = sb, width = 19.69291, height = 12, units = "in", dpi = 450)

# Petiole length ratio (For weather results)

petioleRatio_brach <- ggplot(df_spring_brachwitz) +
  aes(x = group, y = petiole_length_ratio) +
  geom_boxplot(aes(fill = group, color = group), outlier.shape = NA, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.03, colour = "black") +
  scale_fill_manual(values = c("Spring 2023" = "green", "Spring 2024" = "orange", 
                               "Spring 2025" = "red"), name = "Collection") +
  scale_color_manual(values = c("Spring 2023" = "green", "Spring 2024" = "orange", 
                                "Spring 2025" = "red"), name = "Collection") +
  labs(x = NULL, y = "Petiole length ratio") +
  geom_signif(
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    map_signif_level = TRUE, test = "wilcox.test",
    size = 0.3, step_increase = 0.04, y_position = c(0.8),
    textsize = 4, tip_length = 0.005, margin_top = 0.01
  )  +
  theme_pub2()

ggsave(filename = "Desktop/git/Chapter 1/plots/petioleRatio - brachwitz.pdf", 
       plot = petioleRatio_brach, width = 7, height = 6, units = "in", dpi = 450)
