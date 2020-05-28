#https://f1000research.com/articles/5-1408/v3

library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)   # contains the heatmap.2 package
library(tidyr)
library(stringi)

library(limma) 
library(edgeR) 

# RNAseq data:
rnaseq <- fread("data/data_RNA_Seq_v2_expression_median.txt")
rnaseq$Hugo_Symbol <- NULL # because I will use only rnaseq$Entrez_Gene_Id as gene ID in this analysis

rnaseq_genes <- rnaseq$Entrez_Gene_Id  

# Gene lengths for all 20440 IDs
gene_lenght <- fread("data/entrez_gene_length_all_20440ids.txt")
gene_lenght <- gene_lenght[match(rnaseq_genes, gene_lenght$entrez),] # order by rnaseq_genes

#get metadata about samples 
proliferation_df <- fread("data/proliferation_df.csv")

proliferation_samples <-  proliferation_df$Case.ID

rnaseq <-rnaseq[,(names(rnaseq) %in% proliferation_samples ) | (names(rnaseq)=="Entrez_Gene_Id"),with=FALSE] 

rnaseq_mat <- as.matrix(rnaseq[,grepl( "TCGA" , names( rnaseq ) ),with=FALSE]) 

row.names(rnaseq_mat) <- rnaseq$Entrez_Gene_Id

colnames(rnaseq_mat) 

sample_types <- data.table(colnames(rnaseq_mat))

sample_types[sample_types$V1 %in% proliferation_samples ,"type"] <- "proliferation"

#If the counts from all samples were stored in a single file, [LIKE IN Ciriello 2015 cancer paper]
#the data can be read into R and then converted into a DGEList-object using the DGEList function!!!

x <- DGEList(counts=rnaseq_mat,genes=gene_lenght,samples=sample_types, group=sample_types$type )
class(x)

dim(x)

samplenames <- colnames(x) #substring(colnames(x), 1, nchar(colnames(x))-3)  # if you want shorter names


group <- x$samples$group

library(Homo.sapiens) 
geneid <- rownames(x) 
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
dim(genes)

head(genes)

genes <- genes[!duplicated(genes$ENTREZID),]

#In this example, the gene order is the same in both the annotation and the data object. 
#If this is not the case due to missing and/or rearranged gene IDs, 
#the match function can be used to order genes correctly. 
#The data frame of gene annotations is then added to the data object and neatly packaged in a DGEList-object 
#containing raw count data with associated sample information and gene annotations.

##########################
#Data pre-processing
##########################

#Popular transformations include counts per million (rpkm), log2-counts per million (log-rpkm), 
#reads per kilobase of transcript per million (RPKM), and fragments per kilobase of transcript per million (FPKM).

#In our analyses, rpkm and log-rpkm transformations are used regularly although 
#they do not account for gene length differences as RPKM and FPKM values do. 

#gene lengths remain constant for comparisons of interest and any observed differences 
#are a result of changes in condition rather than changes in gene length.

#rpkm function in edgeR
rpkm <- rpkm(x, log=TRUE)

###### How to find "library size"?
#tmp = read.delim("GSE63310_RAW/GSM1545545_JMS9-P8c.txt.gz") #
#sum(tmp$Count) # this is called "library size"
#####

#A rpkm value of 1 for a gene equates to having 20 counts in the sample library size ≈20 million

#log-rpkm values will be used for exploratory plots. When log=TRUE, the rpkm function adds an offset to the rpkm values 
#before converting to the log2-scale. By default, the offset is 2/L where 2 is the “prior count” and 
#L is the average library size in millions, so the log-rpkm values are related to the rpkm values by log2(rpkm + 2/L). 

#The prior count avoids taking the logarithm of zero, and also reduces spurious variability for genes with 
#very low counts by shrinking all the inter-sample log-fold-changes towards zero, something that is helpful for exploratory plotting. 

#For TUTORIAL f1000 dataset, the average library size is about 45.5 million, so L ≈ 45.5 and the minimum log-rpkm value for each sample becomes log2(2/45.5) = −4.51. 
#In other words, a count of zero for this data maps to a log-rpkm value of −4.51 after adding the prior count or offset:

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(rpkm)

#Log-rpkm values are also used in downstream linear modeling via limma’s voom function, 
#although voom recomputes its own log-rpkm values internally with a smaller prior count.


################
#Removing genes that are lowly expressed
################

#filterByExpr function in the edgeR
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#for Fig1 from Tutorial
rpkm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(rpkm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-rpkm")
abline(v=rpkm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(rpkm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
rpkm <- rpkm(x, log=TRUE)
plot(density(rpkm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-rpkm")
abline(v=rpkm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(rpkm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
dim(rpkm)
#legend("topright", samplenames, text.col=col, bty="n")

#

###############
#Normalising gene expression distributions
###############

#calcNormFactors function in edgeR

#When working with DGEList-objects, these normalisation factors are automatically stored in x$samples$norm.factors

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#########

# fig 2
# too many samples for this plot, so that I select a subset of samples 


# diff are amplified for better visual demonstration
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(oma=c(6,2,2,2),mfrow=c(1,2),pch=16) # this setup to fit labels in figures 

rpkm <- rpkm(x2, log=TRUE)
boxplot(rpkm[,1:10], las=2, main="") #boxplot(rpkm, las=2, col=col, main="")
title(main="A. Sample unnormalized data", ylab="Log-rpkm")

x2 <- calcNormFactors(x2)
x2$samples$norm.factors
dim(rpkm)

rpkm <- rpkm(x2, log=TRUE)
boxplot(rpkm[,1:10], las=2, main="") #boxplot(rpkm, las=2, col=col, main="")
title(main="B. Sample normalized data", ylab="Log-rpkm")
###############
save(rpkm, file="data/norm.Rda")
columns <- colnames(rpkm)
save(columns, file="data/colnames.Rda")
rows <- rownames(rpkm)
save(rows, file="data/rownames.Rda")
