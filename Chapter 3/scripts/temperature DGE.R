library(tidyverse)
library(readxl)
library(patchwork)
library(gprofiler2)
library(extrafont)

font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# Loading all counts
load("Desktop/git/Chapter 3/resulting data/DGE/normCounts_all - Location DGE.Rdata")
normCounts_spring <- normCounts_all
load("Desktop/git/Chapter 3/resulting data/DGE/normCounts_all - Seasonal DGE.Rdata")
normCounts <- normCounts_all

#load("Desktop/git/Chapter 3/resulting data/DGE/winter_norm_cpm.Rdata")
#load("Desktop/git/Chapter 3/resulting data/DGE/spring_brach_norm_cpm.Rdata")
#load("Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_norm_cpm.Rdata")

# phenotypes

pheno_2021 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 1)
pheno_2022_mk <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 4)
pheno_2023 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 5)
pheno_2024 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 6)

pheno_2021a <- pheno_2021 %>%
  select(c(1, 9))
pheno_2022a <- pheno_2022_mk %>%
  select(c(1, 2))
pheno_2023a <- pheno_2023 %>%
  select(c(1, 2))
pheno_2024a <- pheno_2024 %>%
  select(c(1, 3))

pheno_all <- rbind(pheno_2021a, pheno_2022a, pheno_2023a, pheno_2024a)
pheno_all$temperature_C <- as.numeric(pheno_all$temperature_C)

# merging all counts

common_norm <- Reduce(
  intersect, list(
    colnames(normCounts),
    colnames(normCounts_spring)
  )
)
winters_norm <- normCounts %>%
  filter(Season == "Winter") %>%
  select(any_of(common_norm))
spring_counts <- normCounts_spring %>%
  select(any_of(common_norm))

normCounts_all <- rbind(
  winters_norm, spring_counts) %>%
  droplevels()

df_all <- merge(pheno_all, normCounts_all, by = "plant_ID") %>%
  distinct() %>%
  drop_na(temperature_C)

# 10 and 90th percentile

prob1 <- quantile(df_all$temperature_C, probs = .1)
prob9 <- quantile(df_all$temperature_C, probs = .9)

df_lower <- df_all %>%
  filter(temperature_C <= prob1) %>%
  droplevels()
df_upper <- df_all %>%
  filter(temperature_C >= prob9) %>%
  droplevels()

# extracting sample counts

gene_lower <- df_lower[, c(7:ncol(df_lower))] %>%
  droplevels() %>%
  t()

gene_upper <- df_upper[, c(7:ncol(df_upper))] %>%
  droplevels() %>%
  t()

# Gene-wise DGE comparison

wilcox_output <- sapply(1:nrow(gene_lower), function(i) {
  wilcox.test(as.numeric(gene_lower[i, ]), 
              as.numeric(gene_upper[i, ]))
})

# apply correction
adjusted_p_vals <- p.adjust(wilcox_output["p.value", ], 
                            method = "BH")

# find the log fold changes
log_fold_changes <- log2(rowMeans(as.matrix(gene_lower)) /
                           rowMeans(as.matrix(gene_upper)))

# extracting and tidying resulting df

deg_res <- data.frame(id = row.names(gene_lower),
                      log2FoldChange = log_fold_changes,
                      p_adj = adjusted_p_vals) %>%
  filter(log2FoldChange != Inf & log2FoldChange != -Inf) %>%
  drop_na(log2FoldChange, p_adj)

deg_res$state <- "NO"
deg_res$state[deg_res$log2FoldChange > 1 & 
                deg_res$p_adj < 0.01] <- "UP"
deg_res$state[deg_res$log2FoldChange < -1 & 
                deg_res$p_adj < 0.01] <- "DOWN"

deg_res$state <- as.factor(deg_res$state)

# Extracting significant DGEs

dge_up <- deg_res %>%
  filter(state == "UP") %>%
  rename(Gene = id) %>%
  arrange(desc(log2FoldChange)) %>%
  droplevels() %>%
  data.frame(row.names = NULL)

dge_down <- deg_res %>%
  filter(state == "DOWN") %>%
  rename(Gene = id) %>%
  arrange(desc(log2FoldChange)) %>%
  droplevels() %>%
  data.frame(row.names = NULL)


# Visualization

dge_plot <- ggplot(data = deg_res, aes(x = log2FoldChange, y = -log10(p_adj), col = state)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed') + 
  geom_point(size = 1.5) + 
  scale_color_manual(values = c("red", "lightgrey", "blue"), 
                     labels = c("Downregulated (Up Upper LST)", "Not significant", "Upregulated (Up Upper LST)")) + 
  coord_cartesian(ylim = c(0, 45), xlim = c(-8, 9)) + 
  labs(color = 'State', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-7, 8, 1)) + 
  labs(title = 'Leaf surface temperature') +
  theme_pub1()

ggsave(filename = "Desktop/git/Chapter 3/plots/Leaf surface temperature DGE - volcano.pdf", 
       plot = dge_plot, width = 8, height = 4, units = "in", dpi = 550)

write.csv(dge_up, "Desktop/git/Chapter 3/resulting data/DGE - results/Leaf surface temperature - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(dge_down, "Desktop/git/Chapter 3/resulting data/DGE - results/Leaf surface temperature - dge down.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(deg_res[, -1], "Desktop/git/Chapter 3/resulting data/DGE - results/Leaf surface temperature - dge all.csv", 
          row.names = FALSE, quote = FALSE)

# GO term analysis of significant DGEs

organism = "athaliana"

up <- as.vector(dge_up$Gene)
down <- as.vector(dge_down$Gene)

# enrichments

gprof_up <- gost(query = up,
                 organism = "athaliana",
                 sources = c("GO:BP"),
                 correction_method = "fdr",
                 user_threshold = 1e-4)

gprof_down <- gost(query = down,
                   organism = "athaliana",
                   sources = c("GO:BP"),
                   correction_method = "fdr",
                   user_threshold = 1e-4)

# tidying GO results

down_term <- gprof_down$result
up_term <- gprof_up$result

up_term$Upregulated <- "Lower LST"
down_term$Upregulated <- "Upper LST"

combined <- rbind(up_term, down_term)

# Unique to each

unique_lower <- setdiff(up_term$term_name, down_term$term_name)
unique_upper <- setdiff(down_term$term_name, up_term$term_name)

# Shared

shared_terms <- intersect(up_term$term_name, down_term$term_name)

# extracting terms for visualization
lower <- up_term %>%
  filter(term_name %in% c(unique_lower)) %>%
  arrange(p_value)
upper <- down_term %>%
  filter(term_name %in% c(unique_upper)) %>%
  arrange(p_value)
shared <- down_term %>%
  filter(term_name %in% c(shared_terms)) %>%
  arrange(p_value)

# visualization (Bar plot)

l <- lower %>%
  mutate(fdr = -log10(p_value)) %>%
  slice_head(n = 25) %>%
  ggpubr::ggbarplot(
    x = "term_name", 
    y = "fdr", 
    fill = "blue", 
    color = "white", palette = "jco", 
    sort.val = "asc", sort.by.groups = TRUE, 
    xlab = "Biological process", ylab = "-log(FDR)", 
    title = "GO-Term Biological Process", rotate = TRUE, 
    legend.title = "DE status", ggtheme = theme_pub1()
  )
u <- upper %>%
  mutate(fdr = -log10(p_value)) %>%
  slice_head(n = 25) %>%
  ggpubr::ggbarplot(
    x = "term_name", 
    y = "fdr", 
    fill = "red", 
    color = "white", palette = "jco", 
    sort.val = "asc", sort.by.groups = TRUE, 
    xlab = "Biological process", ylab = "-log(FDR)", 
    title = "GO-Term Biological Process", rotate = TRUE, 
    legend.title = "DE status", ggtheme = theme_pub1()
  )
c <- shared %>%
  mutate(fdr = -log10(p_value)) %>%
  #slice_head(n = 24) %>%
  ggpubr::ggbarplot(
    x = "term_name", 
    y = "fdr", 
    fill = "black", 
    color = "white", palette = "jco", 
    sort.val = "asc", sort.by.groups = TRUE, 
    xlab = "Biological process", ylab = "-log(FDR)", 
    title = "GO-Term Biological Process", rotate = TRUE, 
    legend.title = "DE status", ggtheme = theme_pub1()
  )

ggsave(filename = "Desktop/git/Chapter 3/plots/Leaf surface temperature DGE - GO lower.pdf", 
       plot = l, width = 8, height = 5, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/Leaf surface temperature - GO upper.pdf", 
       plot = u, width = 8, height = 5, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/Leaf surface temperature - GO shared.pdf", 
       plot = c, width = 8, height = 4, units = "in", dpi = 550)

write.csv(lower[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO lower LST - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(upper[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO upper LST - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(shared[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO shared (LST) - dge up.csv", 
          row.names = FALSE, quote = FALSE)
