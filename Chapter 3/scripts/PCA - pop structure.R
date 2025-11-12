library(tidyverse)
library(patchwork)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# PCA results from TASSEL

winter_pca <- read.delim("Desktop/git/Chapter 3/resulting data/diversity/winter_PCA.txt")
spring_pca <- read.delim("Desktop/git/Chapter 3/resulting data/diversity/spring _PCA.txt")

proport_winter <- read.delim("Desktop/git/Chapter 3/resulting data/diversity/Eigenvalues_winter.txt")
proport_spring <- read.delim("Desktop/git/Chapter 3/resulting data/diversity/Eigenvalues_spring.txt")

# Trait meta information

trait_winter <- read.delim("Desktop/git/Chapter 3/data/diversity/winter.txt")
trait_spring <- read.delim("Desktop/git/Chapter 3/data/diversity/spring.txt")


# Visualizing PCA 

## winter

winter_pca12 <- ggplot(winter_pca) +
  aes(x = PC1, y = PC2) +
  geom_point(size = 2, color = "skyblue") +
  labs(title = "PCA on winter SNPs") +
  xlab("PC1 (28.7%)") +
  ylab("PC2 (7.5%)") +
  theme_pub1()


## Spring

spring_pca2 <- spring_pca %>%
  select(c(1:3)) %>%
  merge(., trait_spring[, c(1, 3, 5)], by = "Taxa") %>%
  mutate_at(c(4, 5), as.factor)

cols <- c("2022" = "blue", "2023" = "grey40", "2024" = "orange")

spring_pca12 <- ggplot(spring_pca2) +
  aes(x = PC1, y = PC2, color = year, shape = Location) +
  geom_point(size = 2) +
  scale_color_manual(values = cols, name = "Year") +
  scale_shape_manual(values = c("Spiekeroog" = 16, "Brachwitz" = 17), name = "Location") +
  labs(title = "PCA on spring SNPs") +
  xlab("PC1 (14.8%)") +
  ylab("PC2 (9.2%)") +
  theme_pub1()

pca_all <- winter_pca12 + spring_pca12

ggsave(filename = "Desktop/git/Chapter 3/plots/SNPs - PCA.pdf", 
       plot = pca_all, width = 13, height = 5, units = "in", dpi = 550)
