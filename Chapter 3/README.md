## Directory structure

- A lot of outputs to add the picture.

## Resulting datasets

- Some of the outputs are huge (over Github limits) especially the gene matrices. 
- These dataset can be accessed in the [git.nfdi4plants.org/eneza-yoeli.mjema](https://git.nfdi4plants.org/eneza-yoeli.mjema/molecular-and-phenotypic-footprints-of-climate-in-native-arabidopsis-thaliana).

## Scripts description

**winter counts - prep.R** and **spring counts - prep.R** - Preparing the raw gene counts from _featureCounts_ to their respective metadata information.

**batch correction - SVA.R** and **normalization CPM.R** - batch coorection and normalization scripts for both winter and spring samples.

**counting checking.R** - PCA plots of raw, batch corrected and normalized counts.

**mapping stats.R** - Visualization script for alignment statistics.

**seasonal DGE.R**, **temperature DGE.R**, **locational DGE.R** and **yearly DGE.R** - Nonparametric differential gene expression analyses with their respctive GO for biological processes analyses of; 
  (1) Winter vs Spring (Spiekeroog samples only)
  (2) Leaf surface temperature measuremnets (All plants divided into two groups as described in the manuscript)
  (3) Spiekeroog vs Brachwitz (Spring collection of the same years i.e 2023 and 2024 plants)
  (4) 2023 vs 2024 spiekeroog spring plants

**PCA - pop structure.R** - PCA plot of SNPs from winter and spring plants.
