library(readxl)
library(tidyverse)

# counts

counts_2022 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/2022samples.Aug2024_3.txt", 
                          header = TRUE, check.names = FALSE)
counts_2025_01 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Jan2025_01_updated3.txt", 
                             header = TRUE, check.names = FALSE)
counts_2024_11 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Dec2024_11_updated3.txt", 
                             header = TRUE, check.names = FALSE)
counts_2024_06 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/allsamples.Jul2024_3.txt", 
                             header = TRUE, check.names = FALSE)
counts_2024_08 <- read.delim("Desktop/git/Chapter 3/data/raw counts/readyCounts/plate4.Aug2024_3.txt", 
                             header = TRUE, check.names = FALSE)

# PCR plates

## Oldenburg extraction

plate1 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/OverviewPlates_RNA-DNA-Protien-Isolation.xlsx", 
                    skip = 7, sheet = 1)
plate2 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/OverviewPlates_RNA-DNA-Protien-Isolation.xlsx", 
                    skip = 7, sheet = 2)
plate3 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/OverviewPlates_RNA-DNA-Protien-Isolation.xlsx", 
                    skip = 7, sheet = 3)
plate4 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/OverviewPlates_RNA-DNA-Protien-Isolation.xlsx", 
                    skip = 7, sheet = 4)
plate5 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/OverviewPlates_RNA-DNA-Protien-Isolation.xlsx", 
                    skip = 7, sheet = 5)
plate6 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/OverviewPlates_RNA-DNA-Protien-Isolation.xlsx", 
                    skip = 7, sheet = 6)

## MLU extraction

MJ2 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 2)
MJ6 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 6)
MJ7 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 7)
MJ8 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/MJ DNA_RNA_Protein isolated Spiekeroog_Brachwitz - 1.xlsx", 
                 skip = 1, sheet = 8)

# Phenotype data

pheno_2021 <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", sheet = 1)

# location and time

site <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Phenotype Measurements - Updated.xlsx", 
                  sheet = 2)
time <- read_xlsx("Desktop/git/Chapter 3/data/metadata/Collection Trip 1222_RNA-DNA-Protien-Isolation.xlsx")

# Tidying the Phenotype data column type 

site_time <- merge(site, time[, c(4, 6, 9)], by = "plant_ID")
site_time$Season <- "Winter"

site_time <- site_time%>%
  mutate_at(c("Location", "Tissue", "Season"), as.factor)

# Merging counts

all_counts <- cbind(counts_2022[, -2], counts_2024_06[, -2], counts_2024_08[, -2], 
                    counts_2024_11[, -2], counts_2025_01[, -2])

allCounts_df <- all_counts %>%
  column_to_rownames(var = "Geneid")

# Understanding the size of shoot and root samples

allCounts_plates <- allCounts_df %>%
  select(any_of(as.character(unique(site_time$plant_ID))))

shoot_root <- site_time %>%
  filter(plant_ID %in% unique(colnames(allCounts_plates))) %>%
  group_by(Tissue) %>%
  summarise(size = n())

# For plate 1, 2, 3, 5, 6 and the MJ plates

## Here all the plates are isolated and batches are created/named. The batches sequence 
## follows the ordering of plates from the ones that were prepared in Oldenburg to Halle.

# IDs from the plates

id_1 <- plate1[, 5]
colnames(id_1)[1] <- "plant_ID"
id_1$Batch <- 1

id_2 <- plate2[, 5]
colnames(id_2)[1] <- "plant_ID"
id_2$Batch <- 2

id_3 <- plate3[, 5]
colnames(id_3)[1] <- "plant_ID"
id_3$Batch <- 3

id_4 <- plate4[, 5]
colnames(id_4)[1] <- "plant_ID"
id_4$Batch <- 4

id_5 <- plate5[, 5]
colnames(id_5)[1] <- "plant_ID"
id_5$Batch <- 5

id_6 <- plate6[, 5]
colnames(id_6)[1] <- "plant_ID"
id_6$Batch <- 6

id_mj2 <- MJ2[, 3]
colnames(id_mj2)[1] <- "plant_ID"
id_mj2$Batch <- 7

id_mj6 <- MJ6[, 3]
colnames(id_mj6)[1] <- "plant_ID"
id_mj6$Batch <- 11

id_mj7 <- MJ7[, 3]
colnames(id_mj7)[1] <- "plant_ID"
id_mj7$Batch <- 12

id_mj8 <- MJ8[, 3]
colnames(id_mj8)[1] <- "plant_ID"
id_mj8$Batch <- 13

# Combining all plates

id_plates <- rbind(id_1, id_2, id_3, id_4, id_5, id_6, id_mj2, 
                   id_mj6, id_mj7, id_mj8)
id_plates$Batch <- as.factor(id_plates$Batch)

# adding batch info to the pheno

siteTime_df <- merge(site_time, id_plates, by = "plant_ID")
siteTime_df$year <- as.factor(2021)

# removing duplicates

siteTime_df2 <- siteTime_df %>%
  distinct()
siteTime_df3 <- siteTime_df2 %>%
  filter(plant_ID != 968)

# Merging plates, phenotype and counts
# Final dataframes

winter_counts <- allCounts_plates %>%
  select(any_of(as.character(siteTime_df3$plant_ID)))

siteTime <- siteTime_df3 %>%
  filter(plant_ID %in% colnames(winter_counts)) %>%
  distinct(plant_ID, .keep_all = TRUE)

# Transposing and merging df

winter_counts_t <- winter_counts %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "plant_ID")

winter_countsMeta <- merge(siteTime, winter_counts_t, by = "plant_ID")

# saving resulting files
save(winter_countsMeta, 
     file = "Desktop/git/Chapter 3/resulting data/DGE/winter_countsMeta.Rdata")

write.csv(winter_countsMeta, "Desktop/git/Chapter 3/resulting data/DGE/winter_countsMeta.csv", 
          row.names = FALSE, quote = FALSE)
