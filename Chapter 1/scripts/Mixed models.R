# Libraries

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggcorrplot)
library(Cairo)
#library(showtext)

#showtext_auto()   # automatically use showtext for all plots
# Register Arial
#font_add("Arial", regular = "Arial.ttf")

source("Desktop/git/functions/custom_theme.R")

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

# Working on spring data

traits <- colnames(data_spring_2[, c(1:10)])

spring_scaled <- data_spring_2 %>%
  mutate(across(all_of(traits), ~ scale(.), .names = "{.col}_s")) %>%
  data.frame(check.names = TRUE)  # scaling measurements to be comparable

traits_scaled <- paste0(traits, "_s")


# mixed model and extracting results 
# Function written with chatgpt assistance

fit_pair <- function(response, predictor, data) {
  form <- as.formula(paste(response, "~", predictor, "+ (1 | year:location)"))
  m <- try(lmer(form, data = data, REML = FALSE), silent = TRUE)
  
  if (inherits(m, "try-error")) {
    return(tibble(estimate = NA_real_, se = NA_real_, p.value = NA_real_))
  }
  
  td <- broom.mixed::tidy(m, effects = "fixed") %>%
    filter(term == predictor)
  
  tibble(
    estimate = td$estimate,
    se = td$std.error,
    p.value = td$p.value
  )
}


results <- expand.grid(response = traits_scaled, predictor = traits_scaled,
                       stringsAsFactors = FALSE) %>%
  filter(response != predictor) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(out = list(fit_pair(response, predictor, spring_scaled))) %>%
  unnest(out) %>%
  ungroup()

results <- results %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         sig = case_when(
           is.na(p.adj) ~ NA_character_,
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ ""
         ))


# Visualizing with ggcorrplot

## estimate pairwise matrix

corr_mat <- results %>%
  select(response, predictor, estimate) %>%
  pivot_wider(names_from = predictor, values_from = estimate) %>%
  column_to_rownames("response") %>%
  as.matrix()

all_traits <- union(rownames(corr_mat), colnames(corr_mat)) #extracting unique trait names
corr_mat <- corr_mat[all_traits, all_traits] # reordering the matrix

corr_mat[!is.finite(corr_mat)] <- 1 # Setting NAs to 1 as they are present when predictor = response

corr_mat_sym <- (corr_mat + t(corr_mat)) / 2 # Symmetrize estimates

## Significance matrix
## Same steps as for etsimates above

p_mat <- results %>%
  select(response, predictor, p.adj) %>%
  pivot_wider(names_from = predictor, values_from = p.adj) %>%
  column_to_rownames("response") %>%
  as.matrix()

p_mat <- p_mat[all_traits, all_traits]
p_mat[!is.finite(p_mat)] <- 0

p_mat_sym <- pmax(p_mat, t(p_mat), na.rm = TRUE)

# Renaming
corr_mat_sym2 <- corr_mat_sym
p_mat_sym2 <- p_mat_sym

colnames(corr_mat_sym2) <- c("Petiole longest leaf", "Stem width", "Width longest leaf", "Leaf number", 
                             "Cauline leaf", "Flowers", "Plant length", "Sideshoots", "Sidebranch", "Length longest leaf")
rownames(corr_mat_sym2) <- c("Petiole longest leaf", "Stem width", "Width longest leaf", "Leaf number", 
                             "Cauline leaf", "Flowers", "Plant length", "Sideshoots", "Sidebranch", "Length longest leaf")

colnames(p_mat_sym2) <- c("Petiole longest leaf", "Stem width", "Width longest leaf", "Leaf number", 
                             "Cauline leaf", "Flowers", "Plant length", "Sideshoots", "Sidebranch", "Length longest leaf")
rownames(p_mat_sym2) <- c("Petiole longest leaf", "Stem width", "Width longest leaf", "Leaf number", 
                             "Cauline leaf", "Flowers", "Plant length", "Sideshoots", "Sidebranch", "Length longest leaf")

# Correlation heatmap

spring <- ggcorrplot(
  corr_mat_sym2, 
  type = "lower", 
  p.mat = p_mat_sym2, 
  sig.level = 0.05, 
  #insig = "blank", 
  lab = TRUE, 
  colors = c("#E46726", "white", "#6D9EC1"),
  ggtheme = theme_pub1(), hc.order = TRUE
)

ggsave(filename = "Desktop/git/Chapter 1/plots/Figure-1d.pdf", 
       plot = spring, width = 7, height = 7, units = "in", dpi = 600, 
       device = cairo_pdf)
