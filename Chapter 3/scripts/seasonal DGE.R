library(tidyverse)
library(patchwork)
library(gprofiler2)
library(extrafont)

font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# normalized counts
# Winter vs Spring using only the spiekeroog collection

load("Desktop/git/Chapter 3/resulting data/DGE/normCounts_all - Seasonal DGE.Rdata")

# extracting seasons

counts_winter <- normCounts_all[normCounts_all$Season == "Winter", 
                                6:ncol(normCounts_all)] %>%
  droplevels() %>%
  t()

counts_spring_spieke <- normCounts_all[normCounts_all$Season == "Spring", 
                                       6:ncol(normCounts_all)] %>%
  droplevels() %>%
  t()

# Gene-wise DGE comparison

## adopted from https://github.com/stressedplants/SinglePlantOmics/blob/main/desync_gene_expr/R/2_DGE.R

wilcox_output <- sapply(1:nrow(counts_winter), function(i) {
  wilcox.test(as.numeric(counts_winter[i, ]), 
              as.numeric(counts_spring_spieke[i, ]))
})

# applying correction for multiple correction

adjusted_p_vals <- p.adjust(wilcox_output["p.value", ], 
                            method = "BH")

# log fold changes

log_fold_changes <- log2(rowMeans(as.matrix(counts_winter)) /
                           rowMeans(as.matrix(counts_spring_spieke)))

# extracting and tidying resulting df

deg_res <- data.frame(id = row.names(counts_winter),
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
                     labels = c("Downregulated (Up spring)", "Not significant", "Upregulated (Up winter)")) + 
  coord_cartesian(ylim = c(0, 215), xlim = c(-7, 10)) + 
  labs(color = 'State', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-7, 9, 1)) + 
  labs(title = 'Winter vs Spring') +
  theme_pub1()

ggsave(filename = "Desktop/git/Chapter 3/plots/Seasonal DGE - volcano.pdf", 
       plot = dge_plot, width = 8, height = 4, units = "in", dpi = 550)

write.csv(dge_up, "Desktop/git/Chapter 3/resulting data/DGE - results/seasonal - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(dge_down, "Desktop/git/Chapter 3/resulting data/DGE - results/seasonal - dge down.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(deg_res[, -1], "Desktop/git/Chapter 3/resulting data/DGE - results/seasonal - dge all.csv", 
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

up_term$Upregulated <- "Winter"
down_term$Upregulated <- "Spring"

combined <- rbind(up_term, down_term)

# Unique to each

unique_winter <- setdiff(up_term$term_name, down_term$term_name)
unique_spring <- setdiff(down_term$term_name, up_term$term_name)

# Shared

shared_terms <- intersect(up_term$term_name, down_term$term_name)

# extracting terms for visualization
winter <- up_term %>%
  filter(term_name %in% c(unique_winter)) %>%
  arrange(p_value)
spring <- down_term %>%
  filter(term_name %in% c(unique_spring)) %>%
  arrange(p_value)
shared <- down_term %>%
  filter(term_name %in% c(shared_terms)) %>%
  arrange(p_value)

# visualization (Bar plot)

w <- winter %>%
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
s <- spring %>%
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

ggsave(filename = "Desktop/git/Chapter 3/plots/Seasonal DGE - GO winter.pdf", 
       plot = w, width = 8, height = 5, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/Seasonal DGE - GO spring.pdf", 
       plot = s, width = 8, height = 5, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/Seasonal DGE - GO shared.pdf", 
       plot = c, width = 8, height = 4, units = "in", dpi = 550)

write.csv(winter[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO winter - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(spring[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO spring - dge up.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(shared[, -14], "Desktop/git/Chapter 3/resulting data/GO terms/GO shared - dge up.csv", 
          row.names = FALSE, quote = FALSE)
