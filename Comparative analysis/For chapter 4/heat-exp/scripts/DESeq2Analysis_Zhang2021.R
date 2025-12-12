#Differential Gene Expression Analysis of Zhang et al. 2021 with DESeq2
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
setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\Zhang2021\\DESeq2")

#Read Information about the samples (group)
#Group 0: heat-treatment for 0 hours, Group 1: heat-treatment for 3 hours 
#Group 2: heat-treatment for 6 hours, Group 3: heat-treatment for 9 hours
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
      #Create DESeqDataSet and pre-filter it
            #######################

#create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ group)
#set 0 as the standard of the factor group
dds$group <- relevel(dds$group,ref="0")
dim(dds)
#Pre-filter
#for each group we have 2 biological replicates
smallestGroupSize <- 2
#at least 10 reads in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)

            ####################### 
      #Transform the data for visualization
            #######################

#ONLY FOR VISUALIZATION WE TRANSFORM THE DATA. FOR THE STATISTICAL ANALYSIS WE USE THE RAW DATA

#Wich transformation to choose?
#vst is good for medium-to-large datasets (numOfSamples > 30)
#rlog is good for small datasets (numOfSamples < 30)

#blind=True if the design of the experiment should be ignored (good for first analysis)

#Variance Stabilizing Transformation (vst)
vsd <- vst(dds, blind=TRUE)
#Regularized-Logarithm (rlog)
rld <- rlog(dds, blind=TRUE)
# this gives log2(n + 1) without variance stabilization for the comparison
ntd <- normTransform(dds)

#Comparing the transformations (the goal is to have a low sd over all genes)
#y-Axis: standard deviation (sd)
#x-Axis: genes sorted by mean of expression (from low to high)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) 
meanSdPlot(assay(rld)) #this one gave the best flat curve --> we use this in the following plots

#Write the plot for rld in png-file
png("Figure_general/meanSdPlot_rld.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld))
dev.off()

            ####################### 
    #First visualizations to have an overview
            ####################### 

#Heatmap of sample-to-sample distances

#calculate the distances between the samples
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
annotation <- data.frame(Group = rld$group)
rownames(annotation) <- colnames(rld)

#Write the plot in a png file
png("Figure_general/heatmap.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = annotation,
         annotation_row = annotation,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

#Principal Components Analysis (PCA) plot

pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("Figure_general/pca.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

           ####################### 
#Differential expression analysis (WITH THE NO-NORMALIZED DATA!!!)
            ####################### 

#Analysis done by DESeq2
dds <- DESeq(dds)

#Dispersion Estimate 
#smooth curve that approaches zero with increasing mean of normalized counts is good
png("Figure_general/dispersion_estimate.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds)
dev.off()

#Exploring the analysis
#Intercept is the log2 for the reference group (group 0)
resultsNames(dds)

#Set alpha and l2fc
alpha <- 0.05
l2fc <- 1

#Creating results by comparing the groups

#LFC is positive if the gene is upregulated in group1/group2/group3
#LFC is negative if the gene is downregulated in group1/group2/group3
res_group1 <- results(dds, name="group_1_vs_0",alpha=alpha)
res_group2 <- results(dds, name="group_2_vs_0",alpha=alpha)
res_group3 <- results(dds, name="group_3_vs_0",alpha=alpha)

#The MA plot shows the mean of the normalized counts versus the l2fc 
#for all genes tested. The genes that are significantly (padj < alpha) are 
#colored to be easily identified.
png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/MA_plots.png",sep=""), width = 18, height = 16, units = "cm", res=300)
par(mfrow = c(1, 3), mar=c(2,2,2,2))
plotMA(res_group1, ylim=c(-2.5,2.5), main="Group 1")
abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_group2, ylim=c(-2.5,2.5), main="Group 2")
abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_group3, ylim=c(-2.5,2.5), main="Group 3")
abline(h=c(-1,1),col="red",lwd=2)
dev.off()

#Identify significant differential expressed genes

#Group 1
signif_group1_all <- res_group1[res_group1$padj < alpha & !is.na(res_group1$padj) & (res_group1$log2FoldChange > l2fc | res_group1$log2FoldChange < -l2fc), ]
signif_group1_all <- sort(as.character(rownames(signif_group1_all)), decreasing = TRUE)
length(signif_group1_all)

#Group 2
signif_group2_all <- res_group2[res_group2$padj < alpha & !is.na(res_group2$padj) & (res_group2$log2FoldChange > l2fc | res_group2$log2FoldChange < -l2fc), ]
signif_group2_all <- sort(as.character(rownames(signif_group2_all)), decreasing = TRUE)
length(signif_group2_all)

#Group 3
signif_group3_all <- res_group3[res_group3$padj < alpha & !is.na(res_group3$padj) & (res_group3$log2FoldChange > l2fc | res_group3$log2FoldChange < -l2fc), ]
signif_group3_all <- sort(as.character(rownames(signif_group3_all)), decreasing = TRUE)
length(signif_group3_all)

#all
signif_all <- union(signif_group1_all,union(signif_group2_all,signif_group3_all))
length(signif_all)

#Write all significant differential expressed genes in a file
write.csv(signif_all,file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/signif_all_zhang2021.csv",sep=""),row.names = FALSE)

#save the data from the analysis
write.csv2(as.data.frame(res_group1),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group1.csv",sep=""))
write.csv2(as.data.frame(res_group2),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group2.csv",sep=""))
write.csv2(as.data.frame(res_group3),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group3.csv",sep=""))
