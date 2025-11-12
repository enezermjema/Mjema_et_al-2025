library(tidyverse)
library(patchwork)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# Loading alignment stats from STAR

log_dir <- list.files("Desktop/git/Chapter 3/data/alignment_stats/", 
                      pattern = "Log.final.out", full.names = TRUE)

# Parsing stat files, extracting relevant info
# chatgpt assisted script

parse_star_log_percentages <- function(file_path) {
  lines <- readLines(file_path)
  metrics <- c(
    "Uniquely mapped reads %" = "unique_pct",
    "% of reads mapped to multiple loci" = "multi_pct",
    "% of reads unmapped: too many mismatches" = "unmap_mismatch_pct",
    "% of reads unmapped: too short" = "unmap_short_pct",
    "% of reads unmapped: other" = "unmap_other_pct"
  )
  
  results <- sapply(names(metrics), function(pattern) {
    line <- lines[grep(pattern, lines, fixed = TRUE)]
    if (length(line) == 0) return(NA)
    val <- sub(".*\\|", "", line)
    val <- trimws(gsub("%", "", val))
    as.numeric(val)
  })
  
  # Extract everything before "Log.final.out"
  sample <- basename(file_path) %>%
    sub("Log\\.final\\.out$", "", .)
  
  tibble(plant_ID = sample, !!!setNames(results, metrics))
}

# Parsed dataframes

all_stats <- map_dfr(log_dir, parse_star_log_percentages)

all_stats <- all_stats %>%
  mutate(unmapped = unmap_mismatch_pct + unmap_short_pct + unmap_other_pct)

# Mapping IDs to the ones used in the DGE analysis

id_counts <- read.csv("Desktop/git/Chapter 3/resulting data/DGE/normCounts_all.csv")

allstats_df <- merge(all_stats[, c(1:3, 7)], 
                     id_counts[, c(1:5)], by = "plant_ID")

# Prep for plotting

long_stats <- allstats_df %>%
  pivot_longer(
    cols = c(unique_pct, multi_pct, unmapped),
    names_to = "category",
    values_to = "percent"
  )

long_stats <- long_stats %>%
  mutate(category = recode(category,
                           unique_pct = "Uniquely mapped",
                           multi_pct = "Multi-mapped",
                           unmapped = "Unmapped"
  ))

long_stats <- long_stats %>%
  mutate(category = factor(category, levels = c(
    "Uniquely mapped",
    "Multi-mapped",
    "Unmapped"
  )))

# visualization

## individual plants

season <- ggplot(long_stats, aes(x = plant_ID, y = percent, fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Season, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL, y = "% of reads mapped", fill = "Mapping category",
    title = "Alignment summary"
  ) +
  theme_pub1() +
  theme(axis.text.x = element_blank())

location <- ggplot(long_stats, aes(x = plant_ID, y = percent, fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Location, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL, y = "Percentage of reads", fill = "Mapping category",
    title = "STAR Mapping Stats by Location"
  ) +
  theme_pub1() +
  theme(axis.text.x = element_blank())

season
location


# Averages
summary_stats <- long_stats %>%
  group_by(Season, Location, category) %>%
  summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop")

desired_order <- c("Winter.Spiekeroog", "Spring.Spiekeroog", "Spring.Brachwitz")

summary_stats <- summary_stats %>%
  mutate(
    Season_Location = interaction(Season, Location),
    Season_Location = factor(Season_Location, levels = desired_order)
  )


average_align <- ggplot(summary_stats, aes(x = Season_Location, 
                                           y = percent, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL, 
    y = "% of average mapped reads",
    fill = "Category", 
    title = "Average mapped reads"
  ) +
  scale_x_discrete(labels = c("Winter-Spiekeroog", 
                            "Spring-Spiekeroog", "Spring-Brachwitz")) +
  theme_pub1()

average_align

ggsave(filename = "Desktop/git/Chapter 3/plots/Figure 3a.pdf", 
       plot = season, width = 7, height = 4, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/mapping stats - location.pdf", 
       plot = location, width = 7, height = 4, units = "in", dpi = 550)
ggsave(filename = "Desktop/git/Chapter 3/plots/mapping stats - means.pdf", 
       plot = average_align, width = 7, height = 4, units = "in", dpi = 550)
