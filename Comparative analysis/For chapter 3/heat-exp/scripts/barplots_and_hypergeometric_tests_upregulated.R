#TASK: create bar plots that show the percentage of the genes found among all upregulated genes for each publication
#TASK: calculate the hypergeometric probabilities
setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress")

#read files
liu <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\Liu\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_liu.csv",head=TRUE)
zhang2021 <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\Zhang2021\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_zhang2021.csv",head=TRUE)
blair <- read.csv("diff_reg_genes_blair.csv",head=FALSE)
liWang <- read.csv("diff_reg_genes_LiWang.csv",head=FALSE)
zhang2017 <- read.csv("all_diff_reg_genes_Zhang2017.csv",head=FALSE)

#print number of differential expressed genes for every publication
length(liu[[1]])
length(zhang2021[[1]])
length(blair[[1]])
length(liWang[[1]])
length(zhang2017[[1]])

#read file with upregulated genes
up <- read.csv2("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/Job/RNASeq/DEGsEnezer/Genes/UP_regulated.csv")
length(up[[1]])

setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/Job/RNASeq/DEGsEnezer/plots/up-heat")

genes_found_in_all_publications = intersect(liu[[1]],intersect(zhang2021[[1]],intersect(blair[[1]],intersect(liWang[[1]],zhang2017[[1]]))))
length(genes_found_in_all_publications)
genes_found_in_all_publications = data.frame(Genes = genes_found_in_all_publications)

genes_found_in_at_least_one = union(liu[[1]],union(zhang2021[[1]],union(blair[[1]],union(liWang[[1]],zhang2017[[1]]))))
length(genes_found_in_at_least_one)
genes_found_in_at_least_one = data.frame(Genes = genes_found_in_at_least_one)

#define function to calculate the mean for the hypergeometric probability
mean_hyper <- function(k,m,N){
  x <- k * m/N
}
#define function to calculate the standard deviation for the hypergeometric probability
sd_hyper <- function(k,m,N){
  sqrt(k * (m/N) * ((N-m)/N) * ((N-k)/(N-1)))
}

my_list <- list(Liu = liu, Zhang2021 = zhang2021, Blair = blair, LiWang = liWang, Zhang2017 = zhang2017, all_publications = genes_found_in_all_publications, at_least_one_publication = genes_found_in_at_least_one)

for(name in names(my_list)){
  #find overlapping genes
  sim <- intersect(my_list[[name]][[1]],up$Gene)
  #Calculate the percentage of the genes we found from all upregulated genes
  percent <- (length(sim) / length(up$Gene)) * 100
  #define parameters for hypergeometric test
  x <- length(sim)
  m <- length(my_list[[name]][[1]])
  n <- 32833-length(my_list[[name]][[1]]) #32833 is the number of genes in Arabidopsis thaliana based on the tx2gene
  k <- 808 #number of upregulated genes
  #Create bar plot
  png(paste("up_all_",name,".png",sep=""), width = 18, height = 16, units = "cm", res=300)
  h <- barplot(c(percent,NA), xlab = NULL, ylab = "Found Genes in %",main=paste("Percentage of genes found among all upregulated genes\n","Number of genes found in ",name," et al.: ",length(my_list[[name]][[1]]),"\nNumber of upregulated genes: ",length(up$Gene), sep=""),cex.main=1,col="lightblue")
  text(h[1], labels = paste(round(percent, digits=1)," %","\n(",length(sim)," Genes)",sep=""), pos = 3)
  #hypergeometric test
  mtext(
    paste("Hypergeometric probability P(X=",length(sim),") = ",dhyper(x = x, m = m, k = k, n = n),sep="", 
          "\nProbability < 0.01 = ", dhyper(x = x, m = m, k = k, n = n) < 0.01,
          "\nMean = ",round(mean_hyper(k=k,m=m,N=32833),digits=2),
          "\nStandard deviation = ", round(sd_hyper(k=k,m=m,N=32833),digits=2)),
    side = 1, 
    line = 3,     
    cex = 0.9     
  )
  dev.off()
}
