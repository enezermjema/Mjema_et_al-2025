            #############################
#create a table, that shows which genes from all predicted genes (for spring) were found in each publication
            #############################

setwd("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\Comparative analysis\\For chapter 4\\heat-exp\\resulting data")
#read files
liu <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\Liu\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_liu.csv",head=TRUE)
zhang2021 <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\Zhang2021\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_zhang2021.csv",head=TRUE)
blair <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\diff_reg_genes_blair.csv",head=FALSE)
liWang <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\diff_reg_genes_LiWang.csv",head=FALSE)
zhang2017 <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Heat_Stress\\all_diff_reg_genes_Zhang2017.csv",head=FALSE)

#print number of differential expressed genes for every publication
length(liu[[1]])
length(zhang2021[[1]])
length(blair[[1]])
length(liWang[[1]])
length(zhang2017[[1]])

#read file from Eneza with all predicted genes for spring
all_predicted <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\Comparative analysis\\For chapter 4\\heat-exp\\data\\genes_spring.csv")
length(all_predicted[[1]])

all_predicted$Liu <- is.element(all_predicted$Genes,liu$x)
sum(all_predicted$Liu)

all_predicted$Zhang2021 <- is.element(all_predicted$Genes,zhang2021$x)
sum(all_predicted$Zhang2021)

all_predicted$Blair <- is.element(all_predicted$Genes,blair$V1)
sum(all_predicted$Blair)

all_predicted$LiWang <- is.element(all_predicted$Genes,liWang$V1)
sum(all_predicted$LiWang)

all_predicted$Zhang2017 <- is.element(all_predicted$Genes,zhang2017$V1)
sum(all_predicted$Zhang2017)

write.csv(all_predicted[c(2,4:8)],file="spring_table.csv",row.names=FALSE)

            #############################
#extract lists of genes from the table, that contain the predicted genes, that were found in 5, 4, 3, 2, 1 or 0 of the publications
            #############################
my_list <- list(spring_in_0 = 0,spring_in_1 = 1,spring_in_2 = 2, spring_in_3 = 3, spring_in_4 = 4, spring_in_5 = 5)
for(name in names(my_list)){
  x <- name
  #filter genes
  name <- all_predicted[rowSums(all_predicted[4:8]) == my_list[[name]],][,2]
  #convert to data frame to name the column
  name <- data.frame(Gene=name)
  #write file
  write.csv(name,file=paste(x,".csv",sep=""),row.names=FALSE,quote=FALSE)
}

spring_in_at_least_1 <- all_predicted[rowSums(all_predicted[4:8]) > 0,][,2]
spring_in_at_least_1 <- data.frame(Gene = spring_in_at_least_1)
write.csv(spring_in_at_least_1,file="spring_in_at_least_1.csv",row.names=FALSE,quote = FALSE)


