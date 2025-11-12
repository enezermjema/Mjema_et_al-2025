library(tidyverse)
library(patchwork)
library(ggrepel)
library(gprofiler2)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# normalized counts
# Winter vs Spring using only the spiekeroog collection

load("Desktop/git/Chapter 3/resulting data/DGE/winter_norm_cpm.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_norm_cpm.Rdata")

# merging counts

common_norm <- Reduce(
  intersect, list(
    colnames(winter_norm_cpm), 
    colnames(spring_spieke_norm_cpm)
  )
)
winters_norm <- winter_norm_cpm %>%
  select(any_of(common_norm)) %>%
  filter(!(plant_ID %in% c(488, 505, 560, 514, 512, 808, 498, 561, 558))) # these samples were flagged as outliers from mahalanobis

spring_spieke_norm <- spring_spieke_norm_cpm %>%
  select(any_of(common_norm)) %>%
  filter(!(plant_ID %in% c("4_1257", "C745", "4_1209", "4_1179")))

normCounts_all <- rbind(
  winters_norm, spring_spieke_norm) %>%
  droplevels()

# filtering low variable genes

counts_all <- normCounts_all[, c(1, 6:ncol(normCounts_all))] %>%
  data.frame(row.names = NULL) %>%
  column_to_rownames("plant_ID") %>%
  mutate_all(~ log(. + 1))

gene_var <- apply(counts_all, 2, var) # computing column variance
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:1000] # top 3000 most variable genes

expr_top <- counts_all[, top_genes]

# PCA

pca_norm <- prcomp(expr_top, scale. = TRUE)

pca12_norm <- as.data.frame(pca_norm$x[, c(1, 2)]) %>%
  data.frame(Season = normCounts_all$Season, 
             year = normCounts_all$year)
pcaNorm_summary <- summary(pca_norm)$importance

###-------------------------------------------------------------------------------------
# extreme outliers
# A PCA with the samples was done prior to computing m distances
# The outliers were later removed (see above)
# You don't have to run this. And if needed run it before sample filtering

mah <- pca12_norm %>%
  rownames_to_column("plant_ID") %>%
  group_by(Season) %>%
  mutate(
    m_dist = mahalanobis(cbind(PC1, PC2),
                         colMeans(cbind(PC1, PC2)),
                         cov(cbind(PC1, PC2)))
  ) %>%
  ungroup()

threshold <- qchisq(0.975, df = 2)  # 97.5% cutoff
outliers <- mah %>% filter(m_dist > threshold)
###-------------------------------------------------------------------------------------

# Biplots

a <- ggplot(pca12_norm) +
  aes(x = PC1, y = PC2, color = normCounts_all$Season) + 
  geom_point(size = 3, alpha = .8) +
  scale_color_manual(values = c("Winter" = "skyblue", "Spring" = "red3"), name = "Season") +
  labs(title = "Winter vs Spring (Spiekeroog, 2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaNorm_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaNorm_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub1()

b <- ggplot(pca12_norm) +
  aes(x = PC1, y = PC2, color = normCounts_all$year) + 
  geom_point() +
  #scale_color_manual(values = c("Spiekeroog" = "skyblue", "Brachwitz" = "brown"), name = "Location") +
  #labs(title = "Normalized batch-corrected expression counts (2021 - 2024)") +
  xlab(paste0("PC", 1, " (", 
              sprintf(pcaNorm_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", 2," (",
              sprintf(pcaNorm_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme_pub1()

# saving

ggsave(filename = "Desktop/git/Chapter 3/plots/Figure 3b.pdf", 
       plot = a, width = 7, height = 4, units = "in", dpi = 550)

save(normCounts_all, file = "Desktop/git/Chapter 3/resulting data/DGE/normCounts_all - Seasonal DGE.Rdata")
