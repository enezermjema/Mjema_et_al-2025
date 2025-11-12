library(tidyverse)
library(GAPIT)

myG <- read.delim("Desktop/git/Chapter 3/resulting data/diversity/winter.hmp.txt", 
                  header = FALSE)
myY <- read.delim("Desktop/git/Chapter 3/data/diversity/winter.txt", 
                  header = TRUE, check.names = FALSE)

setwd("Desktop/git/Chapter 3/resulting data/diversity/GWAS/")

myGAPIT <- GAPIT(
  Y = myY[, c(1, 2)],
  G = myG,
  PCA.total = 2,
  kinship.algorithm = "EMMAX",
  FDRcut = 0.05,
  Model.selection = TRUE,
  model = c("BLINK") # "MLM", , "BLINK", FarmCPU
)

