library(tidyverse)
library(patchwork)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# raw counts
load("Desktop/git/Chapter 3/resulting data/DGE/spring_countsMeta.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/winter_countsMeta.Rdata")

# batch corrected counts
load("Desktop/git/Chapter 3/resulting data/DGE/winter_batch_corrected.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_batch_corrected.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_brach_batch_corrected.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_batch_corrected.Rdata")

# normalized counts
load("Desktop/git/Chapter 3/resulting data/DGE/winter_norm_cpm.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_brach_norm_cpm.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_norm_cpm.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_norm_cpm.Rdata")

# merging location and seasons

## raw counts
common_raw <- Reduce(
  intersect, list(
    colnames(winter_countsMeta), 
    colnames(spring_countsMeta)
  )
)
winters_raw <- winter_countsMeta %>%
  select(any_of(common_raw))
spring_raw <- spring_countsMeta %>%
  select(any_of(common_raw))

rawCounts_all <- rbind(
  winters_raw, spring_raw) %>%
  droplevels()

## batch corrected
common_batch <- Reduce(
  intersect, list(
    colnames(winter_batch_corrected), 
    colnames(spring_brach_batch_corrected),
    colnames(spring_spieke_batch_corrected)
  )
)
winters_batch <- winter_batch_corrected %>%
  select(any_of(common_batch))
spring_brach_batch <- spring_brach_batch_corrected %>%
  select(any_of(common_batch))
spring_spieke_batch <- spring_spieke_batch_corrected %>%
  select(any_of(common_batch))

batchCounts_all <- rbind(
  winters_batch, spring_brach_batch, spring_spieke_batch) %>%
  droplevels()

## norm counts
common_norm <- Reduce(
  intersect, list(
    colnames(winter_norm_cpm), 
    colnames(spring_brach_norm_cpm),
    colnames(spring_spieke_norm_cpm)
  )
)
winters_norm <- winter_norm_cpm %>%
  select(any_of(common_norm))
spring_brach_norm <- spring_brach_norm_cpm %>%
  select(any_of(common_norm))
spring_spieke_norm <- spring_spieke_norm_cpm %>%
  select(any_of(common_norm))

normCounts_all <- rbind(
  winters_norm, spring_brach_norm, spring_spieke_norm) %>%
  droplevels()

# all springs
common_norm2 <- Reduce(
  intersect, list(
    colnames(winter_norm_cpm), 
    colnames(spring_norm_cpm)
  )
)
winters_norm2 <- winter_norm_cpm %>%
  select(any_of(common_norm2))
spring_norm <- spring_norm_cpm %>%
  select(any_of(common_norm2))

normCounts_all2 <- rbind(
  winters_norm2, spring_norm) %>%
  droplevels()

# PCAs
pca_raw <- prcomp(log2(rawCounts_all[, 6:ncol(rawCounts_all)] + 1))
pca_batch <- prcomp(log2(batchCounts_all[, 6:ncol(batchCounts_all)] + 1))
pca_norm <- prcomp(log2(normCounts_all[, 6:ncol(normCounts_all)] + 1))


# Extracting PCA results (explained variance)
pca12_raw <- as.data.frame(pca_raw$x[, c(1, 2)])
pcaRaw_summary <- summary(pca_raw)$importance

pca12_batch <- as.data.frame(pca_batch$x[, c(1, 2)])
pcaBatch_summary <- summary(pca_batch)$importance

pca12_norm <- as.data.frame(pca_norm$x[, c(1, 2)])
pcaNorm_summary <- summary(pca_norm)$importance


# Visualizing (PC biplot)

## batch colors
cols <- c(
  "1" = "red", "2" = "blue", "3" = "green", "4" = "orange", "5" = "purple",
  "6" = "cyan", "7" = "magenta", "8" = "brown", "9" = "pink", "10" = "darkgreen",
  "11" = "skyblue", "12" = "gold", "13" = "darkred", "14" = "black", "15" = "lightblue",
  "16" = "darkorange", "17" = "violet", "18" = "gray", "19" = "tan", "20" = "darkblue"
)


## raw
a <- ggplot(pca12_raw) +
  aes(x = PC1, y = PC2, color = rawCounts_all$Batch) + 
  geom_point() +
  scale_color_manual(values = cols, name = "Batch number") +
  labs(title = "Raw expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaRaw_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaRaw_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub2()

b <- ggplot(pca12_raw) +
  aes(x = PC1, y = PC2, color = rawCounts_all$Season) + 
  geom_point() +
  scale_color_manual(values = c("Winter" = "red", "Spring" = "blue"), name = "Season") +
  labs(title = "Raw expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaRaw_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaRaw_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub2()

c <- ggplot(pca12_raw) +
  aes(x = PC1, y = PC2, color = rawCounts_all$Location) + 
  geom_point() +
  scale_color_manual(values = c("Spiekeroog" = "skyblue", "Brachwitz" = "brown"), name = "Location") +
  labs(title = "Raw expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaRaw_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaRaw_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub2()

# batch
d <- ggplot(pca12_batch) +
  aes(x = PC1, y = PC2, color = batchCounts_all$Batch) + 
  geom_point() +
  scale_color_manual(values = cols, name = "Batch number") +
  labs(title = "Batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaBatch_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaBatch_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub2()

e <- ggplot(pca12_batch) +
  aes(x = PC1, y = PC2, color = batchCounts_all$Season) + 
  geom_point() +
  scale_color_manual(values = c("Winter" = "red", "Spring" = "blue"), name = "Season") +
  labs(title = "Batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaBatch_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaBatch_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub2()

f <- ggplot(pca12_batch) +
  aes(x = PC1, y = PC2, color = batchCounts_all$Location) + 
  geom_point() +
  scale_color_manual(values = c("Spiekeroog" = "skyblue", "Brachwitz" = "brown"), name = "Location") +
  labs(title = "Batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaBatch_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaBatch_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub2()

# norm

g <- ggplot(pca12_norm) +
  aes(x = PC1, y = PC2, color = normCounts_all$Batch) + 
  geom_point() +
  scale_color_manual(values = cols, name = "Batch number", guide = guide_legend(ncol = 2)) +
  labs(title = "Normalized batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaNorm_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaNorm_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub1()

h <- ggplot(pca12_norm) +
  aes(x = PC1, y = PC2, color = normCounts_all$Season) + 
  geom_point() +
  scale_color_manual(values = c("Winter" = "red", "Spring" = "blue"), name = "Season") +
  labs(title = "Normalized batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaNorm_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaNorm_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub1()

i <- ggplot(pca12_norm) +
  aes(x = PC1, y = PC2, color = normCounts_all$Location) + 
  geom_point() +
  scale_color_manual(values = c("Spiekeroog" = "skyblue", "Brachwitz" = "brown"), name = "Location") +
  labs(title = "Normalized batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaNorm_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaNorm_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub1()

dge_counts <- a + d + g + b + e + h + c + f + i

ggsave(filename = "Desktop/git/Chapter 3/plots/Counts - DGE.pdf", 
       plot = dge_counts, width = 18, height = 10, units = "in", dpi = 550)

