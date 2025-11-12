library(tidyverse)
library(DGEobj.utils)
library(edgeR)

# batch corrected counts
load("Desktop/git/Chapter 3/resulting data/DGE/winter_batch_corrected.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_batch_corrected.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_brach_batch_corrected.Rdata")
load("Desktop/git/Chapter 3/resulting data/DGE/spring_batch_corrected.Rdata")

# removing irrelevant genes

winter_batch_corrected_filt <- winter_batch_corrected %>%
  select(-starts_with("G")) %>%
  select(-contains("."))
spring_brach_batch_corrected_filt <- spring_brach_batch_corrected %>%
  select(-starts_with("G")) %>%
  select(-contains("."))
spring_spieke_batch_corrected_filt <- spring_spieke_batch_corrected %>%
  select(-starts_with("G")) %>%
  select(-contains("."))
spring_batch_corrected_filt <- spring_batch_corrected %>%
  select(-starts_with("G")) %>%
  select(-contains("."))

# Normalizing counts
# CPM normalization

## winter

winter_filtered <- winter_batch_corrected_filt %>%
  select(c(1, 8:ncol(winter_batch_corrected_filt))) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  filter(rowMeans(.) > 10)
winter_norm_cpm <- convertCounts(countsMatrix = as.matrix(winter_filtered), 
                                 unit = "CPM")

winter_norm_cpm <- winter_norm_cpm %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "plant_ID") %>%
  merge(winter_batch_corrected_filt[, c(1, 2, 5:7)], ., by = "plant_ID")
  

## spring

spring_brach_filtered <- spring_brach_batch_corrected_filt %>%
  select(c(1, 6:ncol(spring_brach_batch_corrected_filt))) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  filter(rowMeans(.) > 10)
spring_brach_norm_cpm <- convertCounts(countsMatrix = as.matrix(spring_brach_filtered), 
                                 unit = "CPM")   

spring_brach_norm_cpm <- spring_brach_norm_cpm %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "plant_ID") %>%
  merge(spring_brach_batch_corrected_filt[, c(1:5)], ., by = "plant_ID")   # brachwitz

spring_spieke_filtered <- spring_spieke_batch_corrected_filt %>%
  select(c(1, 6:ncol(spring_spieke_batch_corrected_filt))) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  filter(rowMeans(.) > 10)
spring_spieke_norm_cpm <- convertCounts(countsMatrix = as.matrix(spring_spieke_filtered), 
                                       unit = "CPM")    

spring_spieke_norm_cpm <- spring_spieke_norm_cpm %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "plant_ID") %>%
  merge(spring_spieke_batch_corrected_filt[, c(1:5)], ., by = "plant_ID")   # spiekeroog

spring_filtered <- spring_batch_corrected_filt %>%
  select(c(1, 6:ncol(spring_batch_corrected_filt))) %>%
  column_to_rownames(var = "plant_ID") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  filter(rowMeans(.) > 10)
spring_norm_cpm <- convertCounts(countsMatrix = as.matrix(spring_filtered), 
                                        unit = "CPM")  

spring_norm_cpm <- spring_norm_cpm %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "plant_ID") %>%
  merge(spring_batch_corrected_filt[, c(1:5)], ., by = "plant_ID")   # spring 

# saving

save(winter_norm_cpm, file = "Desktop/git/Chapter 3/resulting data/DGE/winter_norm_cpm.Rdata")
save(spring_brach_norm_cpm, file = "Desktop/git/Chapter 3/resulting data/DGE/spring_brach_norm_cpm.Rdata")
save(spring_spieke_norm_cpm, file = "Desktop/git/Chapter 3/resulting data/DGE/spring_spieke_norm_cpm.Rdata")
save(spring_norm_cpm, file = "Desktop/git/Chapter 3/resulting data/DGE/spring_norm_cpm.Rdata")
