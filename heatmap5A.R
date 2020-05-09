library(gplots)

ilc_expression_data <- read.csv("data/ilc_expression_df.csv", header = TRUE)
ilc_expression_df <- data.frame(ilc_expression_data)
ilc_expression_matrix = as.matrix(ilc_expression_df)
cell_colors <- colorpanel(1000,"blue","white","red")

ilc_classification_data <- read.csv("data/ilc_classification_df.csv", header = TRUE)
ilc_classification_df <- data.frame(ilc_classification_data)
ilc_classification_matrix = as.matrix(ilc_classification_df[1] + 1)
col_colors <- colorpanel(3,"red","green","black")
col_color_mapping_function = function(x) col_colors[x]
col_colors_matrix <- apply(ilc_classification_matrix,1,col_color_mapping_function)

heatmap.2(ilc_expression_matrix, scale="row", ColSideColors=col_colors_matrix,
          col=cell_colors, trace="none", density.info="none", labRow=FALSE, labCol=FALSE,
          margin=c(6,3), lhei=c(1,5), lwid=c(1,4), dendrogram="row", Colv=FALSE, main="SAM FDR=0 n=1276, subtypes marked above", xlab="sample", ylab="expression level")
