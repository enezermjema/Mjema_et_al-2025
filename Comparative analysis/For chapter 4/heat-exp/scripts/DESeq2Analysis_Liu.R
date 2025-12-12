#Differential Gene Expression Analysis of Liu et al. 2020 with DESeq2
#TASK: Find the genes that are differential expressed

#Packages
library(DESeq2)
library(readr)
library(tximport)
library(vsn)
library(pheatmap)
library(ggplot2)

            #######################
#Read and prepare file with sample information, quant.sf files and tx2gene file
            #######################
setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\Liu\\DESeq2")

#Read Information about the samples (group)
#Group P: Primed (37°C 1.5 h each day for 1 week), Group PT: Primed and Triggered 
#Group N: Naive (21°C the whole time), Group T: Triggered (45°C 5 min once)
samples <- read.csv("sample_for_DESeq2.csv",header=TRUE)

#factorize group
samples$group <- factor(samples$group)

#Create vector of paths to the quant.sf files
files <- file.path("star_salmon",samples$sample,"quant.sf")
#Name files (we have to do that so that the columns will be named correctly in the DESeqDataSet)
names(files) <- samples$sample

#Read tx2gene
tx2gene <- read_tsv("tx2gene.tsv")

#Read quantification results and transfer them on Gene level
txi <- tximport(files=files,tx2gene=tx2gene,type="salmon")

            #######################
    #Create DESeqDataSets and pre-filter them
            #######################
#We do the analysis separately for the primed plants and the plants that were not primed to find all differential expressed genes

#Primed
keep <- samples$group == "P" | samples$group == "PT"
dds_P <- DESeqDataSetFromTximport(txi = list(abundance = txi$abundance[,keep],counts = txi$counts[,keep], length = txi$length[,keep], countsFromAbundance = txi$countsFromAbundance),
                                   colData = samples[keep, ],
                                   design = ~ group)
dim(dds_P)
#Pre-filter
#for each group we have 2 repetitions
smallestGroupSize <- 2
#at least 10 reads in at least 2 samples
keep_rows <- rowSums(counts(dds_P) >= 10) >= smallestGroupSize
dds_P <- dds_P[keep_rows,]
dim(dds_P)

#Not primed
keep_N <- samples$group == "N" | samples$group == "T"
dds_N <- DESeqDataSetFromTximport(txi = list(abundance = txi$abundance[,keep_N],counts = txi$counts[,keep_N], length = txi$length[,keep_N], countsFromAbundance = txi$countsFromAbundance),
                                   colData = samples[keep_N, ],
                                   design = ~ group)
dim(dds_N)
#Pre-filter
#for each group we have 2 repetitions
smallestGroupSize <- 2
#at least 10 reads in at least 2 samples
keep_rows_N <- rowSums(counts(dds_N) >= 10) >= smallestGroupSize
dds_N <- dds_N[keep_rows_N,]
dim(dds_N)

            ####################### 
        #Transform the data for visualization
            #######################

#ONLY FOR VISUALIZATION WE TRANSFORM THE DATA. FOR THE STATISTICAL ANALYSIS WE USE THE RAW DATA

#Wich transformation to choose?
#vst is good for medium-to-large datasets (numOfSamples > 30)
#rlog is good for small datasets (numOfSamples < 30)

#blind=True if the design of the experiment should be ignored (good for first analysis)

#Variance Stabilizing Transformation (vst)
vsd_P <- vst(dds_P, blind=TRUE)
vsd_N <- vst(dds_N, blind=TRUE)
#Regularized-Logarithm (rlog)
rld_P <- rlog(dds_P, blind=TRUE)
rld_N <- rlog(dds_N, blind=TRUE)
# this gives log2(n + 1) without variance stabilization for the comparison
ntd_P <- normTransform(dds_P)
ntd_N <- normTransform(dds_N)

#Comparing the transformations (the goal is to have a low sd over all genes)
#y-Axis: standard deviation (sd)
#x-Axis: genes sorted by mean of expression (from low to high)
meanSdPlot(assay(ntd_P))
meanSdPlot(assay(vsd_P))
meanSdPlot(assay(rld_P)) # this one gave the best flat curve --> we use this in the following plots

meanSdPlot(assay(ntd_N))
meanSdPlot(assay(vsd_N))
meanSdPlot(assay(rld_N)) # this one gave the best flat curve --> we use this in the following plots

#Write the plots for rld in png-files
png("Figure_general/meanSdPlot_rld_P.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld_P))
dev.off()
png("Figure_general/meanSdPlot_rld_N.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld_N))
dev.off()

            ####################### 
    #First visualizations to have an overview
            ####################### 

#Heatmaps of sample-to-sample distances

#Primed
#calculate the distances between the samples
sampleDists_P <- dist(t(assay(rld_P)))
sampleDistMatrix_P <- as.matrix(sampleDists_P)
annotation_P <- data.frame(Group = rld_P$group)
rownames(annotation_P) <- colnames(rld_P)

#Write the plot in a png file
png("Figure_general/heatmap_P.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix_P,
         clustering_distance_rows = sampleDists_P,
         clustering_distance_cols = sampleDists_P,
         annotation_col = annotation_P,
         annotation_row = annotation_P,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

#Not primed
#calculate the distances between the samples
sampleDists_N <- dist(t(assay(rld_N)))
sampleDistMatrix_N <- as.matrix(sampleDists_N)
annotation_N <- data.frame(Group = rld_N$group)
rownames(annotation_N) <- colnames(rld_N)

#Write the plot in a png file
png("Figure_general/heatmap_N.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix_N,
         clustering_distance_rows = sampleDists_N,
         clustering_distance_cols = sampleDists_N,
         annotation_col = annotation_N,
         annotation_row = annotation_N,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

#Principal Components Analysis (PCA) plot

#Primed
pcaData_P <- plotPCA(rld_P, intgroup=c("group"), returnData=TRUE)
percentVar_P <- round(100 * attr(pcaData_P, "percentVar"))

png("Figure_general/pca_P.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_P, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_P[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_P[2],"% variance")) + 
  coord_fixed()
dev.off()

#Not primed
pcaData_N <- plotPCA(rld_N, intgroup=c("group"), returnData=TRUE)
percentVar_N <- round(100 * attr(pcaData_N, "percentVar"))

png("Figure_general/pca_N.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_N, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_N[2],"% variance")) + 
  coord_fixed()
dev.off()

            ####################### 
#Differential expression analysis (WITH THE NO-NORMALIZED DATA!!!)
            ####################### 

#Analysis done by DESeq2
dds_P <- DESeq(dds_P)
dds_N <- DESeq(dds_N)

#Dispersion Estimate 
#smooth curve that approaches zero with increasing mean of normalized counts is good
png("Figure_general/dispersion_estimate_P.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds_P)
dev.off()

png("Figure_general/dispersion_estimate_N.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds_N)
dev.off()

#Exploring the analysis
#Intercept is the log2 for the reference group
resultsNames(dds_P)
resultsNames(dds_N)

#Set alpha and l2fc
alpha <- 0.05
l2fc <- 1

#Creating results by comparing the groups

#LFC is positive if the gene is upregulated in PT/T
#LFC is negative if the gene is downregulated in PT/T
#set alpha to create correct MA-plots
res_P <- results(dds_P, name = "group_PT_vs_P",alpha=alpha)

res_N <- results(dds_N, name = "group_T_vs_N",alpha=alpha)

#The MA plot shows the mean of the normalized counts versus the l2fc 
#for all genes tested. The genes that are significant (padj < alpha) are 
#colored to be easily identified.

#Primed
png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/MA_plots_P.png",sep=""), width = 18, height = 16, units = "cm", res=300)
plotMA(res_P, ylim=c(-2.5,2.5))
abline(h=c(-1,1),col="red",lwd=2)
dev.off()

#Not primed
png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/MA_plots_N.png",sep=""), width = 18, height = 16, units = "cm", res=300)
plotMA(res_N, ylim=c(-2.5,2.5))
abline(h=c(-1,1),col="red",lwd=2)
dev.off()

#Identify significant differential expressed genes

#Primed
signif_P_all <- res_P[(res_P$padj < alpha & !is.na(res_P$padj)) & (res_P$log2FoldChange > l2fc | res_P$log2FoldChange < -l2fc),]
signif_P_all <- sort(as.character(rownames(signif_P_all)), decreasing = TRUE)
length(signif_P_all)

#Not primed
signif_N_all <- res_N[(res_N$padj < alpha & !is.na(res_N$padj)) & (res_N$log2FoldChange > l2fc | res_N$log2FoldChange < -l2fc),]
signif_N_all <- sort(as.character(rownames(signif_N_all)), decreasing = TRUE)
length(signif_N_all)

#All significant differential expressed genes
signif_all <- union(signif_P_all,signif_N_all)
length(signif_all)

#Write all significant differential expressed genes in a file
write.csv(signif_all,file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/signif_all_liu.csv",sep=""),row.names = FALSE)
#save the data from the analysis
write.csv2(as.data.frame(res_P),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_P.csv",sep=""))
write.csv2(as.data.frame(res_N),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_N.csv",sep=""))
