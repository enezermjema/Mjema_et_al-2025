            #############################
#create a table, that shows which genes from all predicted genes (for winter) were found in each publication
            #############################

setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\Comparative analysis\\For chapter 4\\cold-exp\\resulting data")

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

#read file from Eneza with all predicted genes for winter
all_predicted <- read.csv2("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\Comparative analysis\\For chapter 4\\cold-exp\\data\\genes_winter.csv")
length(all_predicted[[1]])

#build table
all_predicted$Calixto <- is.element(all_predicted$Genes,calixto$x)
sum(all_predicted$Calixto)

all_predicted$Gehan <- is.element(all_predicted$Genes,gehan$x)
sum(all_predicted$Gehan)

all_predicted$Park <- is.element(all_predicted$Genes,park$x)
sum(all_predicted$Park)

all_predicted$Jia <- is.element(all_predicted$Genes,jia$V1)
sum(all_predicted$Jia)

all_predicted$Tiwari <- is.element(all_predicted$Genes,tiwari$V1)
sum(all_predicted$Tiwari)

#write table
write.csv(all_predicted,file="winter_table.csv",row.names = FALSE)

            #############################
#extract lists of genes from the table, that contain the predicted genes, that were found in 5, 4, 3, 2, 1 or 0 of the publications
            #############################
my_list <- list(winter_in_0 = 0,winter_in_1 = 1,winter_in_2 = 2, winter_in_3 = 3, winter_in_4 = 4, winter_in_5 = 5)
for(name in names(my_list)){
  x <- name
  #filter genes
  name <- all_predicted[rowSums(all_predicted[2:6]) == my_list[[name]],][,1]
  #convert to data frame to name the column
  name <- data.frame(Gene=name)
  #write file
  write.csv(name,file=paste(x,"csv",sep="."),row.names=FALSE,quote=FALSE)
}

winter_in_at_least_1 <- all_predicted[rowSums(all_predicted[2:6]) > 0,][,1]
winter_in_at_least_1 <- data.frame(Gene = winter_in_at_least_1)
write.csv(winter_in_at_least_1,file="winter_in_at_least_1.csv",row.names=FALSE,quote = FALSE)


