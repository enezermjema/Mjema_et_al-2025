#Differential Gene Expression Analysis of Calixto et al. 2018 with DESeq2
#TASK: Find the genes that are differential expressed

#Packages
library(DESeq2)
library(readr)
library(tximport)
library(vsn)
library(pheatmap)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(org.At.tair.db)

            #######################
#Read and prepare file with sample information, quant.sf files and tx2gene file
            #######################
setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Calixto\\DESeq2")

#Read Information about the samples (sample,timepoint,repetition,group)
#Group 1: 20°C, Group 2: 4°C Day 1, Group 3: 4°C Day 4
#Repetition indicates the repetition of the experiment (the whole experiment was repeated 3 times)
samples <- read.csv("sample_for_DESeq2.csv",header=TRUE)

#factorize timepoint 
samples$timepoint <- factor(samples$timepoint)
#factorize repetition
samples$repetition <- factor(samples$repetition)
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

#create DESeqDataSet (the factor influencing the counts is the group, that´s why it´s in the design-formula)
dds <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ group)
#set 1 (20°C) as the standard of the factor group
dds$group <- relevel(dds$group,ref=1)

dim(dds)

#Pre-filter
#for each timepoint we have 9 samples because the experiment was repeated 3 times with 3 biological replicates each
smallestGroupSize <- 9
#at least 10 reads in at least 9 samples
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
# this gives log2(n + 1) without variance stabilization for the comparison
ntd <- normTransform(dds)

#Comparing the transformations (the goal is to have a low sd over all genes)
#y-Axis: standard deviation (sd)
#x-Axis: genes sorted by mean of expression (from low to high)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) # this one gave the best flat curve --> we use this in the following plots

#Write the plot for vsd in png-file
png("Figure_general/meanSdPlot_vsd.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd))
dev.off()

            ####################### 
        #First visualizations to have an overview
            ####################### 

#Heatmap of sample-to-sample distances

#calculate the distances between the samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
annotation <- data.frame(Group = vsd$group)
rownames(annotation) <- colnames(vsd)

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

#For group1
pcaData_group1 <- plotPCA(vsd[,vsd$group == 1], intgroup=c("timepoint"), returnData=TRUE)
percentVar_group1 <- round(100 * attr(pcaData_group1, "percentVar"))
#For group2
pcaData_group2 <- plotPCA(vsd[,vsd$group == 2], intgroup=c("timepoint"), returnData=TRUE)
percentVar_group2 <- round(100 * attr(pcaData_group2, "percentVar"))
#For group 3
pcaData_group3 <- plotPCA(vsd[,vsd$group == 3], intgroup=c("timepoint"), returnData=TRUE)
percentVar_group3 <- round(100 * attr(pcaData_group3, "percentVar"))

#group 1
png("Figure_general/pca_group1.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_group1, aes(PC1, PC2, color=timepoint)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_group1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_group1[2],"% variance")) + 
  coord_fixed()
dev.off()

#group 2
png("Figure_general/pca_group2.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_group2, aes(PC1, PC2, color=timepoint)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_group2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_group2[2],"% variance")) + 
  coord_fixed()
dev.off()

#group 3
png("Figure_general/pca_group3.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData_group3, aes(PC1, PC2, color=timepoint)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar_group3[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_group3[2],"% variance")) + 
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
#Intercept is the log2 for the reference group (group 1)
resultsNames(dds)

#Set alpha and l2fc
alpha <- 0.05
l2fc <- 1

#Creating results by comparing the groups

#LFC is positive if the gene is upregulated in group2/group3
#LFC is negative if the gene is downregulated in group2/group3
res_group2 <- results(dds, name="group_2_vs_1",alpha=alpha)
res_group3 <- results(dds, name="group_3_vs_1",alpha=alpha)

#The MA plot shows the mean of the normalized counts versus the l2fc 
#for all genes tested. The genes that are significantly (padj < alpha) are 
#colored to be easily identified.
png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/MA_plots.png",sep=""), width = 18, height = 16, units = "cm", res=300)
par(mfrow = c(1, 2), mar=c(2,2,2,2))
plotMA(res_group2, ylim=c(-2.5,2.5), main="Group 2")
abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_group3, ylim=c(-2.5,2.5), main="Group 3")
abline(h=c(-1,1),col="red",lwd=2)
dev.off()

#Identify significant differential expressed genes

#Group 2
signif_group2_all <- res_group2[res_group2$padj < alpha & !is.na(res_group2$padj) & (res_group2$log2FoldChange > l2fc | res_group2$log2FoldChange < -l2fc), ]
signif_group2_all <- sort(as.character(rownames(signif_group2_all)), decreasing = TRUE)
length(signif_group2_all)

#Group 3
signif_group3_all <- res_group3[res_group3$padj < alpha & !is.na(res_group3$padj) & (res_group3$log2FoldChange > l2fc | res_group3$log2FoldChange < -l2fc), ]
signif_group3_all <- sort(as.character(rownames(signif_group3_all)), decreasing = TRUE)
length(signif_group3_all)

#all
signif_all <- union(signif_group2_all,signif_group3_all)
length(signif_all)

#Write all significant differential expressed genes in a file
write.csv(signif_all,file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/signif_all_calixto.csv",sep=""),row.names = FALSE)

#save the data from the analysis
write.csv2(as.data.frame(res_group2),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group2.csv",sep=""))
write.csv2(as.data.frame(res_group3),file=paste("Genes_alpha=",alpha,"_l2fc=",l2fc,"/res_group3.csv",sep=""))

            ####################### 
            #GO Enrichment Analysis
            ####################### 

####Group 2####

#Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes_group2 <- sort(as.character(rownames(res_group2)), decreasing = TRUE)

#Extract significant results: padj < alpha and log2FoldChange > l2fc or < -l2fc
signif_res_group2_up <- res_group2[(res_group2$padj < alpha & !is.na(res_group2$padj) & res_group2$log2FoldChange > l2fc), ]
signif_res_group2_up <- sort(as.character(rownames(signif_res_group2_up)), decreasing = TRUE)

signif_res_group2_down <- res_group2[(res_group2$padj < alpha & !is.na(res_group2$padj) & res_group2$log2FoldChange < -l2fc), ]
signif_res_group2_down <- sort(as.character(rownames(signif_res_group2_down)), decreasing = TRUE)

#Run GO enrichment analysis UP
ego_group2_up <- enrichGO(gene = signif_res_group2_up,
                       universe = all_genes_group2,
                       OrgDb = "org.At.tair.db",
                       keyType = "TAIR",
                       ont = "BP", #biological process
                       pAdjustMethod = "BH",
                       readable = TRUE)

dotplot(ego_group2_up, showCategory = 20)
upsetplot(ego_group2_up, n = 10)

png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/enrich_group2_up.png",sep=""), width = 18, height = 16, units = "cm", res=300)
barplot(ego_group2_up, showCategory = 15, title = "Group 2 up")
dev.off()

#Run GO enrichment analysis DOWN  
ego_group2_down <- enrichGO(gene = signif_res_group2_down,
                         universe = all_genes_group2,
                         OrgDb = "org.At.tair.db",
                         keyType = "TAIR",
                         ont = "BP", #biological process
                         pAdjustMethod = "BH",
                         readable = TRUE)

dotplot(ego_group2_down, showCategory = 20)
upsetplot(ego_group2_down, n = 10)

png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/enrich_group2_down.png",sep=""), width = 18, height = 16, units = "cm", res=300)
barplot(ego_group2_down, showCategory = 15, title = "Group 2 down")
dev.off()

####Group 3####

#Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes_group3 <- sort(as.character(rownames(res_group3)), decreasing = TRUE)

#Extract significant results: padj < 0.05 and log2FoldChange > 0.5 or < -0.5
signif_res_group3_up <- res_group3[(res_group3$padj < alpha & !is.na(res_group3$padj) & res_group3$log2FoldChange > l2fc), ]
signif_res_group3_up <- sort(as.character(rownames(signif_res_group3_up)), decreasing = TRUE)

signif_res_group3_down <- res_group3[(res_group3$padj < alpha & !is.na(res_group3$padj) & res_group3$log2FoldChange < -l2fc), ]
signif_res_group3_down <- sort(as.character(rownames(signif_res_group3_down)), decreasing = TRUE)

#Run GO enrichment analysis UP
ego_group3_up <- enrichGO(gene = signif_res_group3_up,
                          universe = all_genes_group3,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_group3_up, showCategory = 20)
upsetplot(ego_group3_up, n = 10)

png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/enrich_group3_up.png",sep=""), width = 18, height = 16, units = "cm", res=300)
barplot(ego_group3_up, showCategory = 15, title = "Group 3 up")
dev.off()

#Run GO enrichment analysis DOWN  
ego_group3_down <- enrichGO(gene = signif_res_group3_down,
                            universe = all_genes_group3,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_group3_down, showCategory = 20)
upsetplot(ego_group3_down, n = 10)

png(paste("Figure_alpha=",alpha,"_l2fc=",l2fc,"/enrich_group3_down.png",sep=""), width = 18, height = 16, units = "cm", res=300)
barplot(ego_group3_down, showCategory = 15, title = "Group 3 down")
dev.off()