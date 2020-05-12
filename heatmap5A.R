library(gplots)
library("devtools")
library(RColorBrewer)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


ilc_expression_data <- read.csv("data/ilc_expression_df.csv", header = TRUE)
ilc_expression_df <- data.frame(ilc_expression_data)
ilc_expression_matrix = as.matrix(ilc_expression_df)
cell_colors <- colorpanel(1000,"blue","white","red")

ilc_classification_data <- read.csv("data/ilc_classification_df.csv", header = TRUE)
ilc_classification_df <- data.frame(ilc_classification_data)
ilc_pam50_matrix = as.matrix(ilc_classification_df[1])
ilc_subtype_matrix = as.matrix(ilc_classification_df[4])
ilc_purity_matrix = as.matrix(ilc_classification_df[2])
ilc_proliferation_matrix = as.matrix(ilc_classification_df[3])
bluecols = brewer.pal(9, 'Oranges')
newcol <- colorRampPalette(bluecols)
ncols <- 100
bluecols2 <- newcol(ncols)
proliferation_colors <- bluecols2
proliferation_color_mapping_function = function(x) proliferation_colors[x]
proliferation_colors_matrix <- apply(ilc_proliferation_matrix,1,proliferation_color_mapping_function)
subtype_colors <- colorpanel(3,"red","green","black")
subtype_color_mapping_function = function(x) subtype_colors[x]
subtype_colors_matrix <- apply(ilc_subtype_matrix,1,subtype_color_mapping_function)
pam50_colors <- c("brown","blue","cyan","pink","yellow")
pam50_color_mapping_function = function(x) pam50_colors[x]
pam50_colors_matrix <- apply(ilc_pam50_matrix,1,pam50_color_mapping_function)
clab=cbind(proliferation_colors_matrix, pam50_colors_matrix, subtype_colors_matrix)
colnames(clab)=c("Proliferation", "PAM50", "ILC Class")

pdf("img/heatmap5A.pdf",width=5,height=10)
heatmap.3(ilc_expression_matrix, scale="row", ColSideColors=clab,
          col=cell_colors, trace="none", density.info="none", labRow=FALSE, labCol=FALSE,
          margin=c(6,3), lhei=c(1,2), lwid=c(1,4), dendrogram="row", Colv=FALSE,
          xlab="sample", ylab="SAM FDR=0 n=1276",
          distfun=function(x) as.dist(1-cor(t(x))),hclustfun=function(x) hclust(x, method="ward.D2") )

legend("topright",legend=c("Immune-related","Proliferative","Reactive-like","", "Basal","Her2","LumA","LumB","Normal"),
       fill=c("red","green","black", "white", "brown","blue","cyan","pink","yellow"), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

dev.off()

