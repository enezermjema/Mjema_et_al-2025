library(tidyverse)
library(patchwork)
library(gprofiler2)
library(extrafont)

font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# Loading spring counts from Spiekeroog and Brachwitz

load("Desktop/git/Chapter 3/resulting data/DGE/normCounts_all - Location DGE.Rdata")

# extracting locations

counts_spieke <- normCounts_all[normCounts_all$Location == "Spiekeroog", 
                                6:ncol(normCounts_all)] %>%
  droplevels() %>%
  t()

counts_brach <- normCounts_all[normCounts_all$Location == "Brachwitz", 
                                       6:ncol(normCounts_all)] %>%
  droplevels() %>%
  t()

# Gene-wise DGE comparison

## adopted from https://github.com/stressedplants/SinglePlantOmics/blob/main/desync_gene_expr/R/2_DGE.R

wilcox_output <- sapply(1:nrow(counts_spieke), function(i) {
  wilcox.test(as.numeric(counts_spieke[i, ]), 
              as.numeric(counts_brach[i, ]))
})

# applying correction for multiple correction

adjusted_p_vals <- p.adjust(wilcox_output["p.value", ], 
                            method = "BH")

# log fold changes

log_fold_changes <- log2(rowMeans(as.matrix(counts_spieke)) /
                           rowMeans(as.matrix(counts_brach)))

# extracting and tidying resulting df

deg_res <- data.frame(id = row.names(counts_spieke),
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
                     labels = c("Downregulated (Up Brachwitz)", "Not significant", "Upregulated (Up Spiekeroog)")) + 
  labs(color = 'State', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-7, 9, 1)) + 
  labs(title = 'Spiekeroog vs Brachwitz') +
  theme_pub1()

ggsave(filename = "Desktop/git/Chapter 3/plots/Locational DGE - volcano.pdf", 
       plot = dge_plot, width = 8, height = 4, units = "in", dpi = 550)

write.csv(dge_up, "Desktop/git/Chapter 3/resulting data/DGE - results/Locational - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(dge_down, "Desktop/git/Chapter 3/resulting data/DGE - results/Locational - dge down.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(deg_res[, -1], "Desktop/git/Chapter 3/resulting data/DGE - results/Locational - dge all.csv", 
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

up_term$Upregulated <- "Spiekeroog"
down_term$Upregulated <- "Brachwitz"

combined <- rbind(up_term, down_term)

# Unique to each

unique_spieke <- setdiff(up_term$term_name, down_term$term_name)
unique_brach <- setdiff(down_term$term_name, up_term$term_name)

# Shared

shared_terms <- intersect(up_term$term_name, down_term$term_name)

# extracting terms for visualization
spieke <- up_term %>%
  filter(term_name %in% c(unique_spieke)) %>%
  arrange(p_value)
brach <- down_term %>%
  filter(term_name %in% c(unique_brach)) %>%
  arrange(p_value)
shared <- down_term %>%
  filter(term_name %in% c(shared_terms)) %>%
  arrange(p_value)

# visualization (Bar plot)

s <- spieke %>%
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
b <- brach %>%
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

ggsave(filename = "Desktop/git/Chapter 3/plots/Locational DGE - GO spiekeroog.pdf", 
       plot = s, width = 8, height = 5, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/Locational DGE - GO brachwitz.pdf", 
       plot = b, width = 8, height = 5, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/Locational DGE - GO shared.pdf", 
       plot = c, width = 8, height = 4, units = "in", dpi = 550)

write.csv(spieke[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO spiekeroog - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(brach[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO brachwitz - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(shared[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO shared (Locational) - dge up.csv", 
          row.names = FALSE, quote = FALSE)
