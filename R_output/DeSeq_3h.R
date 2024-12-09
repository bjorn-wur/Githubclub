if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("devtools")
BiocManager::install("pachterlab/sleuth")
BiocManager::install('EnhancedVolcano')
install.packages(pheatmap)
install.packages('dplyr')
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("glue")
#load packages
library(DESeq2)
library(pheatmap)
library(dplyr)
#generate PCA
library("FactoMineR")
library("factoextra")
#vulcanoplot
library(EnhancedVolcano)
library(glue)
#load dataset
expression_data=read.table("gene_count_matrix_3h.csv", row.names=1,
                           header=TRUE, sep =",", stringsAsFactors=FALSE)

colnames(expression_data)
#subset for Xag8ra 3 hour treatment
Xag8ra_treatment<- c("SRR24630893","SRR24630894","SRR24630896","SRR24630912","SRR24630913","SRR24630914")
#subset for Xanthomonas camperstris 3 hour treatment
XCV3_treatment<- c("SRR24630931","SRR24630932","SRR24630933","SRR24630977","SRR24630978","SRR24630979")
#subset for TMV 4 hours, using Xag8ra mock
TMV_treatment<- c("SRR24630893","SRR24630894","SRR24630896","SRR25031738","SRR25031739","SRR25031740")
#subset for p.capsici
p_capsici_treatment<- c("SRR24630947","SRR24630948","SRR24630949","SRR24630967","SRR24630968","SRR24630969")


#subset data
Xag8ra_subset <- expression_data[Xag8ra_treatment]
p_capsici_subset <- expression_data[p_capsici_treatment]
TMV_subset <- expression_data[TMV_treatment]
XCV3_subset <- expression_data[XCV3_treatment]


colnames(Xag8ra_subset)<- c("mock_1","mock_2","mock_3","Xag8ra_1","Xag8ra_2","Xag8ra_3")
colnames(XCV3_subset)<- c("mock_1","mock_2","mock_3","XCV3_1","XCV3_2","XCV3_3")
colnames(p_capsici_subset)<- c("mock_1","mock_2","mock_3","P_capsici_1","P_capsici_2","P_capsici_3")
colnames(TMV_subset) <-c("mock_1","mock_2","mock_3","TMV_1","TMV_2","TMV_3")


#check data
dim(Xag8ra_subset)
dim(XCV3_subset)
dim(p_capsici_subset)
dim(TMV_subset)

#remove all na values
Xag8ra_subset <- na.omit(Xag8ra_subset)
XCV3_subset <- na.omit(XCV3_subset)
p_capsici_subset<- na.omit(p_capsici_subset)
TMV_subset<- na.omit(TMV_subset)

#filter out zero values
Xag8ra_subset <- filter_if(Xag8ra_subset, is.numeric, all_vars((.) != 0))
XCV3_subset <- filter_if(XCV3_subset, is.numeric, all_vars((.) != 0))
p_capsici_subset<- filter_if(p_capsici_subset, is.numeric, all_vars((.) != 0))
TMV_subset<- filter_if(TMV_subset, is.numeric, all_vars((.) != 0))

dim(Xag8ra_subset)
dim(XCV3_subset)
dim(p_capsici_subset)
dim(TMV_subset)
#Run these commands to create a DESeqDataSet data object
condition = factor(c("mock","mock","mock","sample","sample","sample"),c("mock","sample"))
col_data = data.frame(condition)
#add conditions


treatments <- list(Xag8ra_subset,XCV3_subset,TMV_subset,p_capsici_subset)
names(treatments) <- c("Xag8ra_3h","XCV3_3h","TMV_4h","p_capsici_4h")
for (i in seq_along(treatments)) {
  subset <- treatments[[i]]
  subset_name <- names(treatments)[i]
  print(subset_name)
  class(subset_name)
  }

for (i in seq_along(treatments)) {
  subset <- treatments[[i]]
  subset_name <- names(treatments)[i]
  print(subset_name)
  class(subset_name)
  name <- glue("{subset_name}")
  dds = DESeqDataSetFromMatrix(subset, col_data, ~condition)
  #Normalization, calculate the (linear) correction factors for each sample:
  dds = estimateSizeFactors(dds)
  #The correction (size) factors can then be retrieved using:
  sizeFactors(dds)
  #t transform the counts to get a more normal distribution;
  rld = rlog(dds)
  #check data after log transformation:
  plot(density(assay(rld)[,1]), main="log counts")
  #To estimate the dispersions, run:
  dds = estimateDispersions(dds)
  plotDispEsts(dds)
  dds = nbinomWaldTest(dds)
  #To get a table with differential expression values for the genes, type:
  res = results(dds)
  #For a small number of genes no padj could be calculated, these have no
  #value which will cause problems later, so set them to 1 with this command:
  res$padj = ifelse(is.na(res$padj), 1, res$padj)
  #MA plot
  plotMA(res, main="MA plot",ylim=c(-8,8),alpha=0.01)
  #Enhanced vulcano bacteria
  print(EnhancedVolcano(res,
                  title = name,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'padj',
                  FCcutoff = 1,
                  pCutoff = 0.05))
  #turn res into a dataframe
  output_table<- as.data.frame(res)
  #filter for significant genes <0.05
  output_filtered <- output_table[output_table$padj<0.05,]
  #filter for fold change
  output_filtered_up <- output_filtered[output_filtered$log2FoldChange>1,]
  output_filtered_down <-output_filtered[output_filtered$log2FoldChange < -1,]
  #order data of upregulated genes by adjusted p-value
  output_up_sorted <- output_filtered_up[order(output_filtered_up$padj),]
  output_down_sorted <- output_filtered_down[order(output_filtered_down$padj),]
  #write bacteria result to tsv file
  filename_up= glue("{subset_name}_upregulated.tsv")
  filename_down = glue("{subset_name}_downregulated.tsv")
  write.table(output_filtered_up, col.names=NA, row.names=T, file=as.character(filename_up), sep ="\t")
  write.table(output_filtered_down, col.names=NA, row.names=T, file =as.character(filename_down), sep="\t")
}



