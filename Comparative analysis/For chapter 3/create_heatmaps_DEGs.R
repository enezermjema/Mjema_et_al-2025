library(ggplot2)
library(patchwork)

make_heatmap <- function(file, title, colors, bar_width = 1) {
  #read file
  res <- read.csv(file)
  #calculate the group for every gene predicted
  res$Group <- rowSums(res[2:6])
  #create matrix that contains the number of genes for every group
  mat_new <- as.matrix(table(res$Group))
  
  #convert matrix to data frame
  gruppen <- as.numeric(rownames(mat_new))
  haeufigkeit <- mat_new[,1]
  
  data <- data.frame(
    Gruppe = factor(gruppen),
    Haeufigkeit = haeufigkeit
  )
  #like this the label '1.23 %' fits in the box
  data$Haeufigkeit[data$Gruppe == "5"] <- data$Haeufigkeit[data$Gruppe == "5"] * 2
  #change order of levels
  data$Gruppe <- factor(data$Gruppe, levels = rev(levels(data$Gruppe)))
  #calculate y position for every field
  data$y_position <- cumsum(data$Haeufigkeit) - (data$Haeufigkeit / 2)
  
  #calculate the percentage of genes for every group
  percentage <- floor(mat_new[,1] / length(res[[1]]) * 100 * 100) / 100
  
  #create the plot (based on the theme from Eneza)
  #it was not possible to use one of the custom themes because no axes are needed
  ggplot(data, aes(x = title, y = y_position)) +
    geom_tile(aes(fill = Gruppe, height = Haeufigkeit, width = bar_width), color = "black") +
    geom_text(aes(label = paste(percentage, "%")), color = "black", size = 2.5, family="Arial") +
    scale_fill_manual(values = colors) +
    labs(title = title, fill = "Group") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 10 - 1, family = "Arial"),
      legend.title = element_text(size = 10, face = "bold", family = "Arial"),
      # Titles
      plot.title = element_text(size = 10 + 2, face = "bold", hjust = 0.5, family = "Arial", colour = "black"),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

down_plot <- make_heatmap(
  "C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\DEGsEnezer\\Results\\down_table.csv", "Winter/Down",
  c("skyblue4", "skyblue3", "skyblue2", "skyblue1", "skyblue", "lightgrey"),
)

up_plot <- make_heatmap(
  "C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\DEGsEnezer\\Results\\up_table.csv", "Spring/Up",
  c("coral4", "coral3", "coral2", "coral1","coral", "lightgrey"),
)

#plot the two plots side by side with separate legends
combined <- down_plot + up_plot + plot_layout(nrow = 1)

ggsave("C:\\Users\\cajus\\OneDrive\\Dokumente\\Bildung\\Studium\\Job\\RNASeq\\DEGsEnezer\\Results\\heatmap_up_and_down.svg", combined, width = 12, height = 12, units = "cm", dpi = 300)

