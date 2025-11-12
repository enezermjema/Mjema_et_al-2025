library(tidyverse)
library(sva)

# Loading data

load("Desktop/git/Chapter 3/resulting data/DGE/spring_countsMeta.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/winter_countsMeta.Rdata")

# winter dataset

winters_df <- winter_countsMeta %>%
  dplyr::filter(Tissue == "SHOOT" & Location == "Spiekeroog") %>%
  droplevels()

winter_counts <- winters_df[, c(1, 8:ncol(winters_df))] %>%
  column_to_rownames("plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE)

# spring dataset

spring_spieke <- spring_countsMeta[spring_countsMeta$Location == "Spiekeroog", 
                                   c(1, 6:ncol(spring_countsMeta))] %>%
  data.frame(row.names = NULL) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  droplevels()
spring_spieke_pheno <- spring_countsMeta[spring_countsMeta$Location == "Spiekeroog", c(1:5)] %>%
  data.frame(row.names = NULL) %>%
  column_to_rownames(var = "plant_ID") %>%
  droplevels()

spring_brach <- spring_countsMeta[spring_countsMeta$Location == "Brachwitz", 
                                  c(1, 6:ncol(spring_countsMeta))] %>%
  data.frame(row.names = NULL) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  droplevels()
spring_brach_pheno <- spring_countsMeta[spring_countsMeta$Location == "Brachwitz", c(1:5)] %>%
  data.frame(row.names = NULL) %>%
  column_to_rownames(var = "plant_ID") %>%
  droplevels()

# SVA batch correction

# winter batches
batch_winter <- ComBat_seq(
  counts = as.matrix(winter_counts), 
  batch = winters_df$Batch
)
batch_winter_df <- batch_winter %>%
  t() %>%
  data.frame(check.names = FALSE)

winter_batch_corrected <- batch_winter_df %>%
  rownames_to_column("plant_ID") %>%
  merge(winters_df[, c(1:7)], ., by = "plant_ID")

# Spring for each location
batch_spring_spieke <- ComBat_seq(
  counts = as.matrix(spring_spieke), 
  batch = spring_spieke_pheno$Batch
)
batch_spring_spieke_df <- batch_spring_spieke %>%
  t() %>%
  data.frame(check.names = FALSE)

spring_spieke_batch_corrected <- batch_spring_spieke_df %>%
  rownames_to_column("plant_ID") %>%
  merge(spring_countsMeta[, c(1:5)], ., by = "plant_ID") %>%
  droplevels()

#-----

batch_spring_brach <- ComBat_seq(
  counts = as.matrix(spring_brach), 
  batch = spring_brach_pheno$Batch
)
batch_spring_brach_df <- batch_spring_brach %>%
  t() %>%
  data.frame(check.names = FALSE)

spring_brach_batch_corrected <- batch_spring_brach_df %>%
  rownames_to_column("plant_ID") %>%
  merge(spring_countsMeta[, c(1:5)], ., by = "plant_ID") %>%
  droplevels()

# All spring samples

spring_countsMeta_BS <- spring_countsMeta %>%
  filter(Location %in% c("Spiekeroog", "Brachwitz")) %>%
  droplevels()

spring_counts <- spring_countsMeta_BS[, c(1, 6:ncol(spring_countsMeta_BS))] %>%
  data.frame(row.names = NULL) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE)

covar_mat <- cbind(spring_countsMeta_BS$Location, spring_countsMeta_BS$year)

batch_spring <- ComBat_seq(
  counts = as.matrix(spring_counts), 
  batch = spring_countsMeta_BS$Batch,
  group = NULL,
  covar_mod = covar_mat
)

batch_spring_df <- batch_spring %>%
  t() %>%
  data.frame(check.names = FALSE)

spring_batch_corrected <- batch_spring_df %>%
  rownames_to_column("plant_ID") %>%
  merge(spring_countsMeta_BS[, c(1:5)], ., by = "plant_ID") %>%
  droplevels()

# saving

save(winter_batch_corrected, file = "Desktop/git/Chapter 3/resulting data/DGE/winter_batch_corrected.Rdata")
save(spring_spieke_batch_corrected, file = "Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_batch_corrected.Rdata")
save(spring_brach_batch_corrected, file = "Desktop/git/Chapter 3/resulting data/DGE/spring_brach_batch_corrected.Rdata")
save(spring_batch_corrected, file = "Desktop/git/Chapter 3/resulting data/DGE/spring_batch_corrected.Rdata")
