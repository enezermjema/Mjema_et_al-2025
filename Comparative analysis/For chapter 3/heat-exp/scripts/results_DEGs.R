setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/Job/RNASeq/DEGsEnezer")

            #############################
#create a table, that shows which DEGs (downregulated) were found in each publication that examined cold treatment
            #############################
down <- read.csv2("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/Job/RNASeq/DEGsEnezer/Genes/DOWN_regulated.csv")

#read files
calixto <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Calixto\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_calixto.csv")
gehan <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Gehan\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_gehan.csv")
park <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\Park\\DESeq2\\Genes_alpha=0.05_l2fc=1\\signif_all_park.csv")
jia <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\all_diff_reg_genes_jia.csv", head=FALSE)
tiwari <- read.csv("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\EnezerML\\Cold_Stress\\all_diff_reg_genes_tiwari.csv", head=FALSE)

#print number of differential expressed genes for every publication
length(calixto[[1]])
length(gehan[[1]])
length(park[[1]])
length(jia[[1]])
length(tiwari[[1]])      

#create columns that contain the information if the DEG was found in the specific publication (cold-treatment)
down$Calixto <- is.element(down$Gene,calixto$x)
sum(down$Calixto)

down$Gehan <- is.element(down$Gene,gehan$x)
sum(down$Gehan)

down$Park <- is.element(down$Gene,park$x)
sum(down$Park)

down$Jia <- is.element(down$Gene,jia$V1)
sum(down$Jia)

down$Tiwari <- is.element(down$Gene,tiwari$V1)
sum(down$Tiwari)

write.csv(down[c(1,5:9)],file="down_table.csv",row.names = FALSE)

            #############################
#extract lists of genes from the table (downregulated), that contain the predicted genes, that were found in 5, 4, 3, 2, 1 or 0 of the publications
            #############################
my_list <- list(down_in_0 = 0,down_in_1 = 1,down_in_2 = 2, down_in_3 = 3, down_in_4 = 4, down_in_5 = 5)
for(name in names(my_list)){
  x <- name
  #filter genes
  name <- down[rowSums(down[5:9]) == my_list[[name]],][,1]
  #convert to data frame to name the column
  name <- data.frame(Gene=name)
  #write file
  write.csv(name,file=paste("Results/",x,".csv",sep=""),row.names=FALSE,quote=FALSE)
}

down_in_at_least_1 <- down[rowSums(down[5:9]) > 0,][,1]
down_in_at_least_1 <- data.frame(Gene = down_in_at_least_1)
write.csv(down_in_at_least_1,file="Results/down_in_at_least_1.csv",row.names=FALSE,quote = FALSE)




            #############################
#create a table, that shows which DEGs (upregulated) were found in each publication that examined heat treatment
            #############################
up <- read.csv2("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/Job/RNASeq/DEGsEnezer/Genes/UP_regulated.csv")

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

#create columns that contain the information if the DEG was found in the specific publication (heat-treatment)
up$Liu <- is.element(up$Gene,liu$x)
sum(up$Liu)

up$Zhang2021 <- is.element(up$Gene,zhang2021$x)
sum(up$Zhang2021)

up$Blair <- is.element(up$Gene,blair$V1)
sum(up$Blair)

up$LiWang <- is.element(up$Gene,liWang$V1)
sum(up$LiWang)

up$Zhang2017 <- is.element(up$Gene,zhang2017$V1)
sum(up$Zhang2017)

write.csv(up[c(1,5:9)],file="up_table.csv",row.names=FALSE)

            #############################
#extract lists of genes from the table (upregulated), that contain the predicted genes, that were found in 5, 4, 3, 2, 1 or 0 of the publications
            #############################
my_list <- list(up_in_0 = 0,up_in_1 = 1,up_in_2 = 2, up_in_3 = 3, up_in_4 = 4, up_in_5 = 5)
for(name in names(my_list)){
  x <- name
  #filter genes
  name <- up[rowSums(up[5:9]) == my_list[[name]],][,1]
  #convert to data frame to name the column
  name <- data.frame(Gene=name)
  #write file
  write.csv(name,file=paste("Results/",x,".csv",sep=""),row.names=FALSE,quote=FALSE)
}

up_in_at_least_1 <- up[rowSums(up[5:9]) > 0,][,1]
up_in_at_least_1 <- data.frame(Gene = up_in_at_least_1)
write.csv(up_in_at_least_1,file="Results/up_in_at_least_1.csv",row.names=FALSE,quote = FALSE)
