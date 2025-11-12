library(readxl)
library(tidyverse)

# counts

counts_2025_01 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Jan2025_01_updated3.txt", 
                             header = TRUE, check.names = FALSE)
counts_2024_11 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Dec2024_11_updated3.txt", 
                             header = TRUE, check.names = FALSE)
counts_2024_06 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Jul2024_3.txt", 
                             header = TRUE, check.names = FALSE)
counts_2025_03 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Apr2025_03_edited.txt", 
                             header = TRUE, check.names = FALSE)
counts_2025_05 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.May2025_03.txt", 
                             header = TRUE, check.names = FALSE)

# PCR plates

## MLU extraction

MJ2 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 2)
MJ3 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 3)
MJ4 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 4)
MJ5 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 5)
MJ6 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 6)
MJ14 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 14)
MJ15 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 15)
MJ16 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 16)
MJ17 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 17)
MJ19 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 19)
MJ20 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 20)
MJ22 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                  skip = 1, sheet = 22)

# Phenotype data

pheno_2022_mk <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 4)
pheno_2023 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 5)
pheno_2024 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 6)

# Tidying the Phenotype data column type 

pheno_2022_mk <- pheno_2022_mk %>%
  mutate_at(c("root_length_cm", "length_longest_leaf_cm", "petiole_longest_leaf_cm", 
              "width_longest_leaf_cm", "leaves_number", "temperature_C", "longitude_GMS", 
              "latitude_GMS", "altitude", "Internodium_cm", "stem_mm", "Cauline_leaves", 
              "flowers", "siliques_number", "plant_length_cm", "Mainshoot", "Sideshoots", 
              "Sidebranch"), as.numeric) %>%
  select(-time)

pheno_2023 <- pheno_2023 %>%
  mutate_at(c("root_length_cm", "length_longest_leaf_cm", "petiole_longest_leaf_cm", 
              "width_longest_leaf_cm", "leaves_number", "temperature_C", "longitude_GMS", 
              "latitude_GMS", "altitude", "Internodium_cm", "stem_mm", "Cauline_leaves", 
              "flowers", "siliques_number", "plant_length_cm", "Mainshoot", "Sideshoots", 
              "Sidebranch"), as.numeric) %>%
  select(-time)

pheno_2024 <- pheno_2024 %>%
  mutate_at(c("root_length_cm", "length_longest_leaf_cm", "petiole_longest_leaf_cm", 
              "width_longest_leaf_cm", "leaves_number", "temperature_C", "longitude_GMS", 
              "latitude_GMS", "altitude", "Internodium_cm", "stem_mm", "Cauline_leaves", 
              "flowers", "siliques_number", "plant_length_cm", "Mainshoot", "Sideshoots", 
              "Sidebranch"), as.numeric) %>%
  select(-time)

# Merging phenotype colnames

same_names <- Reduce(intersect, list(colnames(pheno_2022_mk), 
                                     colnames(pheno_2023), colnames(pheno_2024)))

traits_spring_2022 <- pheno_2022_mk %>%
  select(any_of(same_names))

traits_spring_2023 <- pheno_2023 %>%
  select(any_of(same_names))

traits_spring_2024 <- pheno_2024 %>%
  select(any_of(same_names))

pheno_df <- rbind(traits_spring_2022[, -c(21, 22)], 
                  traits_spring_2023[, -c(21, 22)], traits_spring_2024[, -c(21, 22)])

# Merging counts

all_counts <- cbind(counts_2024_06[, -2], counts_2024_11[, -2], 
                    counts_2025_01[, -2], counts_2025_03[, -2], counts_2025_05[, -2])

allCounts_df <- all_counts %>%
  column_to_rownames(var = "Geneid")

# IDs from the plates

id_mj2 <- MJ2[, 3]
colnames(id_mj2)[1] <- "plant_ID"
id_mj2$Batch <- 7

id_mj3 <- MJ3[, 3]
colnames(id_mj3)[1] <- "plant_ID"
id_mj3$Batch <- 8

id_mj4 <- MJ4[, 3]
colnames(id_mj4)[1] <- "plant_ID"
id_mj4$Batch <- 9

id_mj5 <- MJ5[, 3]
colnames(id_mj5)[1] <- "plant_ID"
id_mj5$Batch <- 10

id_mj6 <- MJ6[, 3]
colnames(id_mj6)[1] <- "plant_ID"
id_mj6$Batch <- 11

id_mj14 <- MJ14[, 3]
colnames(id_mj14)[1] <- "plant_ID"
id_mj14$Batch <- 14

id_mj15 <- MJ15[, 3]
colnames(id_mj15)[1] <- "plant_ID"
id_mj15$Batch <- 15

id_mj16 <- MJ16[, 3]
colnames(id_mj16)[1] <- "plant_ID"
id_mj16$Batch <- 16

id_mj17 <- MJ17[, 3]
colnames(id_mj17)[1] <- "plant_ID"
id_mj17$Batch <- 17

id_mj19 <- MJ19[, 3]
colnames(id_mj19)[1] <- "plant_ID"
id_mj19$Batch <- 18

id_mj20 <- MJ20[, 3]
colnames(id_mj20)[1] <- "plant_ID"
id_mj20$Batch <- 19

id_mj22 <- MJ22[, 3]
colnames(id_mj22)[1] <- "plant_ID"
id_mj22$Batch <- 20

# Combining all plates

id_plates <- rbind(id_mj2, id_mj3, id_mj4, id_mj5, id_mj6, 
                   id_mj14, id_mj15, id_mj16, id_mj17, 
                   id_mj19, id_mj20, id_mj22)

id_plates$Batch <- as.factor(id_plates$Batch)

pheno_plates <- merge(pheno_df, id_plates, by = "plant_ID")

# Merging plates, phenotype and counts
# Final dataframes

spring_counts <- allCounts_df %>%
  select(any_of(as.character(pheno_plates$plant_ID)))

spring_pheno <- pheno_plates %>%
  filter(plant_ID %in% unique(colnames(spring_counts))) %>%
  select(c("plant_ID", "season", "year", "Batch", "location")) %>%
  droplevels()

# Transposing and merging df

spring_counts_t <- spring_counts %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "plant_ID")

spring_countsMeta <- merge(spring_pheno, spring_counts_t, by = "plant_ID") %>%
  distinct() %>%
  mutate_at(c(2:5), as.factor) %>%
  rename(
    Location = location,
    Season = season
  )

# saving resulting files
save(spring_countsMeta, 
     file = "Desktop/git/Chapter 3/resulting data/DGE/spring_countsMeta.Rdata")

write.csv(spring_countsMeta, "Desktop/git/Chapter 3/resulting data/DGE/spring_countsMeta.csv", 
          row.names = FALSE, quote = FALSE)
