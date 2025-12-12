#TASK: create bar plots that show the percentage of the genes found among all predicted genes for each publication
#TASK: calculate the hypergeometric probabilities

setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\Comparative analysis\\For chapter 4\\cold-exp\\results")

#read files
calixto <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Calixto\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_calixto.csv")
gehan <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Gehan\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_gehan.csv")
park <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Park\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_park.csv")
jia <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\signif_all_jia.csv", head=FALSE)
tiwari <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\signif_all_tiwari.csv", head=FALSE)

#print number of differential expressed genes for every publication
length(calixto[[1]])
length(gehan[[1]])
length(park[[1]])
length(jia[[1]])
length(tiwari[[1]])

all_predicted_ML <- read.csv2("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\Comparative analysis\\For chapter 4\\cold-exp\\data\\genes_winter.csv")
length(all_predicted_ML[[1]])

genes_found_in_all_publications = intersect(calixto[[1]],intersect(gehan[[1]],intersect(park[[1]],intersect(jia[[1]],tiwari[[1]]))))
length(genes_found_in_all_publications)
genes_found_in_all_publications = data.frame(Genes = genes_found_in_all_publications)

genes_found_in_at_least_one = union(calixto[[1]],union(gehan[[1]],union(park[[1]],union(jia[[1]],tiwari[[1]]))))
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

my_list <- list(Calixto = calixto, Gehan = gehan, Park = park, Jia = jia, Tiwari = tiwari, all_publications = genes_found_in_all_publications, at_least_one_publication = genes_found_in_at_least_one)

for(name in names(my_list)){
  #find overlapping genes
  sim <- intersect(my_list[[name]][[1]],all_predicted_ML$Genes)
  #Calculate the percentage of the genes we found from all predicted genes from Enezer
  percent <- (length(sim) / length(all_predicted_ML$Genes)) * 100
  #define parameters for hypergeometric test
  x <- length(sim)
  m <- length(my_list[[name]][[1]])
  n <- 32833-length(my_list[[name]][[1]]) #32833 is the number of genes in Arabidopsis thaliana based on the tx2gene
  k <- 195 #number of genes predicted for winter
  #Create bar plot
  png(paste("predicted_all_",name,".png",sep=""), width = 18, height = 16, units = "cm", res=300)
  h <- barplot(c(percent,NA), xlab = NULL, ylab = "Found Genes in %",main=paste("Percentage of genes found among all predicted genes (alpha=0.05 and l2fc=1)\n","Number of genes found in ",name," et al.: ",length(my_list[[name]][[1]]),"\nNumber of genes predicted: ",length(all_predicted_ML$Genes), sep=""),cex.main=1,col="lightblue")
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
