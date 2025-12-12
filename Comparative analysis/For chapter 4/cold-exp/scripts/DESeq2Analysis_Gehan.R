#Differential Gene Expression Analysis of Gehan et al. 2015 with DESeq2
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
setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Gehan\\DESeq2")

#Read Information about the samples (sample, repetition, group, accession)
#Group 0: warm, Group 1: 1-week-cold Group 2: 2-week-cold
#Repetition indicates the repetition of the experiment (the whole experiment was repeated 3 times)
#Accession indicates whether the plant is from Italy or Sweden
samples <- read.csv("sample_for_DESeq2.csv",header=TRUE)

#factorize repetition
samples$repetition <- factor(samples$repetition)
#factorize group
samples$group <- factor(samples$group)
#factorize accession
samples$accession <- factor(samples$accession)

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
#We do the analysis separately for the IT and SW accession to find all differential expressed genes

#IT
keep <- samples$accession == "IT"
dds_IT <- DESeqDataSetFromTximport(txi = list(abundance = txi$abundance[,keep],counts = txi$counts[,keep], length = txi$length[,keep], countsFromAbundance = txi$countsFromAbundance),
                                   colData = samples[keep, ],
                                   design = ~ group)
dim(dds_IT)
#Pre-filter
#for each group we have 3 repetitions
smallestGroupSize <- 3
#at least 10 reads in at least 3 samples
keep_rows <- rowSums(counts(dds_IT) >= 10) >= smallestGroupSize
dds_IT <- dds_IT[keep_rows,]
dim(dds_IT)


#SW
keep_SW <- samples$accession == "SW"
dds_SW <- DESeqDataSetFromTximport(txi = list(abundance = txi$abundance[,keep_SW],counts = txi$counts[,keep_SW], length = txi$length[,keep_SW], countsFromAbundance = txi$countsFromAbundance),
                                   colData = samples[keep_SW, ],
                                   design = ~ group)
dim(dds_SW)
#Pre-filter
#for each group we have 3 repetitions
smallestGroupSize <- 3
#at least 10 reads in at least 3 samples
keep_rows_SW <- rowSums(counts(dds_SW) >= 10) >= smallestGroupSize
dds_SW <- dds_SW[keep_rows_SW,]
dim(dds_SW)


            ####################### 
      #Transform the data for visualization
            #######################

#ONLY FOR VISUALIZATION WE TRANSFORM THE DATA. FOR THE STATISTICAL ANALYSIS WE USE THE RAW DATA

#Wich transformation to choose?
#vst is good for medium-to-large datasets (numOfSamples > 30)
#rlog is good for small datasets (numOfSamples < 30)

#blind=True if the design of the experiment should be ignored (good for first analysis)

#Variance Stabilizing Transformation (vst)
vsd_IT <- vst(dds_IT, blind=TRUE)
vsd_SW <- vst(dds_SW, blind=TRUE)
#Regularized-Logarithm (rlog)
rld_IT <- rlog(dds_IT, blind=TRUE)
rld_SW <- rlog(dds_SW, blind=TRUE)
# this gives log2(n + 1) without variance stabilization for the comparison
ntd_IT <- normTransform(dds_IT)
ntd_SW <- normTransform(dds_SW)

#Comparing the transformations (the goal is to have a low sd over all genes)
#y-Axis: standard deviation (sd)
#x-Axis: genes sorted by mean of expression (from low to high)
meanSdPlot(assay(ntd_IT))
meanSdPlot(assay(vsd_IT))
meanSdPlot(assay(rld_IT)) # this one gave the best flat curve --> we use this in the following plots

meanSdPlot(assay(ntd_SW))
meanSdPlot(assay(vsd_SW))
meanSdPlot(assay(rld_SW)) # this one gave the best flat curve --> we use this in the following plots

#Write the plots for rld in png-files
png("Figure_general/meanSdPlot_rld_IT.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld_IT))
dev.off()
png("Figure_general/meanSdPlot_rld_SW.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld_SW))
dev.off()

            ####################### 
      #First visualizations to have an overview
            ####################### 

#Heatmaps of sample-to-sample distances

#IT
#calculate the distances between the samples
sampleDists_IT <- dist(t(assay(rld_IT)))
sampleDistMatrix_IT <- as.matrix(sampleDists_IT)
annotation_IT <- data.frame(Group = rld_IT$group)
rownames(annotation_IT) <- colnames(rld_IT)

#Write the plot in a png file
png("Figure_general/heatmap_IT.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix_IT,
         clustering_distance_rows = sampleDists_IT,
         clustering_distance_cols = sampleDists_IT,
         annotation_col = annotation_IT,
         annotation_row = annotation_IT,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

#SW
#calculate the distances between the samples
sampleDists_SW <- dist(t(assay(rld_SW)))
sampleDistMatrix_SW <- as.matrix(sampleDists_SW)
annotation_SW <- data.frame(Group = rld_SW$group)
rownames(annotation_SW) <- colnames(rld_SW)

#Write the plot in a png file
png("Figure_general/heatmap_SW.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix_SW,
         clustering_distance_rows = sampleDists_SW,
         clustering_distance_cols = sampleDists_SW,
         annotation_col = annotation_SW,
         annotation_row = annotation_SW,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()


#Principal Components Analysis (PCA) plot

#IT
pcaData_IT <- plotPCA(rld_IT, intgroup=c("group"), returnData=TRUE)
percentVar_IT <- round(100 * attr(pcaData_IT, "percentVar"))

png("Figure_general/pca_IT.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_IT, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_IT[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_IT[2],"% variance")) + 
  coord_fixed()
dev.off()

#SW
pcaData_SW <- plotPCA(rld_SW, intgroup=c("group"), returnData=TRUE)
percentVar_SW <- round(100 * attr(pcaData_SW, "percentVar"))

png("Figure_general/pca_SW.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_SW, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_SW[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_SW[2],"% variance")) + 
  coord_fixed()
dev.off()

            ####################### 
#Differential expression analysis (WITH THE NO-NORMALIZED DATA!!!)
            ####################### 

#Analysis done by DESeq2
dds_IT <- DESeq(dds_IT)
dds_SW <- DESeq(dds_SW)

#Dispersion Estimate 
#smooth curve that approaches zero with increasing mean of normalized counts is good
png("Figure_general/dispersion_estimate_IT.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds_IT)
dev.off()

png("Figure_general/dispersion_estimate_SW.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds_SW)
dev.off()

#Exploring the analysis
#Intercept is the log2 for the reference group (group 0)
resultsNames(dds_IT)
resultsNames(dds_SW)

#Set alpha and l2fc
alpha <- 0.05
l2fc <- 1

#Creating results by comparing the groups

#LFC is positive if the gene is upregulated in group1/group2
#LFC is negative if the gene is downregulated in group1/group2
#set alpha to create correct MA-plots
res_group1_IT <- results(dds_IT, name = "group_1_vs_0",alpha=alpha)
res_group2_IT <- results(dds_IT, name = "group_2_vs_0",alpha=alpha)

res_group1_SW <- results(dds_SW, name = "group_1_vs_0",alpha=alpha)
res_group2_SW <- results(dds_SW, name = "group_2_vs_0",alpha=alpha)


#The MA plot shows the mean of the normalized counts versus the l2fc 
#for all genes tested. The genes that are significant (padj < alpha) are 
#colored to be easily identified.

#IT
png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/MA_plots_IT.png",sep=""), width = 18, height = 16, units = "cm", res=300)
par(mfrow = c(1, 2), mar=c(2,2,2,2))
plotMA(res_group1_IT, ylim=c(-2.5,2.5), main="Group 1")
abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_group2_IT, ylim=c(-2.5,2.5), main="Group 2")
abline(h=c(-1,1),col="red",lwd=2)
dev.off()

#SW
png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/MA_plots_SW.png",sep=""), width = 18, height = 16, units = "cm", res=300)
par(mfrow = c(1, 2), mar=c(2,2,2,2))
plotMA(res_group1_SW, ylim=c(-2.5,2.5), main="Group 1")
abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_group2_SW, ylim=c(-2.5,2.5), main="Group 2")
abline(h=c(-1,1),col="red",lwd=2)
dev.off()

#Identify significant differential expressed genes

#Group 1
#IT
signif_group1_IT_all <- res_group1_IT[(res_group1_IT$padj < alpha & !is.na(res_group1_IT$padj)) & (res_group1_IT$log2FoldChange > l2fc | res_group1_IT$log2FoldChange < -l2fc),]
signif_group1_IT_all <- sort(as.character(rownames(signif_group1_IT_all)), decreasing = TRUE)
length(signif_group1_IT_all)

#SW
signif_group1_SW_all <- res_group1_SW[(res_group1_SW$padj < alpha & !is.na(res_group1_SW$padj)) & (res_group1_SW$log2FoldChange > l2fc | res_group1_SW$log2FoldChange < -l2fc),]
signif_group1_SW_all <- sort(as.character(rownames(signif_group1_SW_all)), decreasing = TRUE)
length(signif_group1_SW_all)

#Group 2
#IT
signif_group2_IT_all <- res_group2_IT[(res_group2_IT$padj < alpha & !is.na(res_group2_IT$padj)) & (res_group2_IT$log2FoldChange > l2fc | res_group2_IT$log2FoldChange < -l2fc),]
signif_group2_IT_all <- sort(as.character(rownames(signif_group2_IT_all)), decreasing = TRUE)
length(signif_group2_IT_all)

#SW
signif_group2_SW_all <- res_group2_SW[(res_group2_SW$padj < alpha & !is.na(res_group2_SW$padj)) & (res_group2_SW$log2FoldChange > l2fc | res_group2_SW$log2FoldChange < -l2fc),]
signif_group2_SW_all <- sort(as.character(rownames(signif_group2_SW_all)), decreasing = TRUE)
length(signif_group2_SW_all)


#All significant differential expressed genes
signif_all <- union(signif_group1_IT_all,union(union(signif_group1_SW_all,signif_group2_IT_all),signif_group2_SW_all))
#Significant differential expressed genes in IT plants
signif_all_IT <- union(signif_group1_IT_all,signif_group2_IT_all)
#Significant differential expressed genes in SW plants
signif_all_SW <- union(signif_group1_SW_all,signif_group2_SW_all)

#Write all significant differential expressed genes in a file
write.csv(signif_all,file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/signif_all_gehan.csv",sep=""),row.names = FALSE)

#save the data from the analysis
write.csv2(as.data.frame(res_group1_IT),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group1_IT.csv",sep=""))
write.csv2(as.data.frame(res_group2_IT),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group2_IT.csv",sep=""))
write.csv2(as.data.frame(res_group1_SW),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group1_SW.csv",sep=""))
write.csv2(as.data.frame(res_group2_SW),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group2_SW.csv",sep=""))
