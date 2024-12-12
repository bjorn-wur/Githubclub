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
expression_data=read.table("gene_count_matrix.csv", row.names=1,
                           header=TRUE, sep =",", stringsAsFactors=FALSE)

colnames(expression_data)
#subset for Xag8ra 3 hour treatment
Xag8ra_treatment<- c("SRR24630893","SRR24630894","SRR24630896","SRR24630912","SRR24630913","SRR24630914")
#subset for Xanthomonas camperstris strain 3 (XCV3) 3 hour treatment
XCV3_treatment<- c("SRR24630931","SRR24630932","SRR24630933","SRR24630977","SRR24630978","SRR24630979")
#subset for Xanthomonas camperstris strain 1 (XCV1) 3 hour treatment
XCV1_treatment<- c("SRR24630931","SRR24630932","SRR24630933","SRR24630881","SRR24630882","SRR24630883")
#subset for TMV 4 hours, using Xag8ra mock
TMV_treatment<- c("SRR24630893","SRR24630894","SRR24630896","SRR25031738","SRR25031739","SRR25031740")
#subset for p.capsici 4 hours
p_capsici_treatment<- c("SRR24630947","SRR24630948","SRR24630949","SRR24630967","SRR24630968","SRR24630969")
#subset for mocks to check virus 
mock_treatment <-c("SRR24630893","SRR24630894","SRR24630896","SRR24630931","SRR24630932","SRR24630933","SRR24630947","SRR24630948","SRR24630949")

#subset data
Xag8ra_subset <- expression_data[Xag8ra_treatment]
p_capsici_subset <- expression_data[p_capsici_treatment]
TMV_subset <- expression_data[TMV_treatment]
XCV3_subset <- expression_data[XCV3_treatment]
XCV1_subset <- expression_data[XCV1_treatment]
mock_subset <- expression_data[mock_treatment]


colnames(Xag8ra_subset)<- c("mock_1","mock_2","mock_3","Xag8ra_1","Xag8ra_2","Xag8ra_3")
colnames(XCV3_subset)<- c("mock_1","mock_2","mock_3","XCV3_1","XCV3_2","XCV3_3")
colnames(XCV1_subset)<- c("mock_1","mock_2","mock_3","XCV1_1","XCV1_2","XCV1_3")
colnames(p_capsici_subset)<- c("mock_1","mock_2","mock_3","P_capsici_1","P_capsici_2","P_capsici_3")
colnames(TMV_subset)<-c("mock_1","mock_2","mock_3","TMV_1","TMV_2","TMV_3")
colnames(mock_subset)<-c("Xag8ra_m1","Xag8ra_m2","Xag8ra_m3","XCV3_m1","XCV3_m2","XCV3_m3","p_capsici_m1","p_capsici_m2","p_capsici_m3")

#check data
dim(Xag8ra_subset)
dim(XCV3_subset)
dim(p_capsici_subset)
dim(TMV_subset)
dim(mock_subset)

#remove all na values
Xag8ra_subset <- na.omit(Xag8ra_subset)
XCV3_subset <- na.omit(XCV3_subset)
XCV1_subset <- na.omit(XCV1_subset)
p_capsici_subset<- na.omit(p_capsici_subset)
TMV_subset<- na.omit(TMV_subset)
expressiondata_nona <- na.omit(expression_data)

#filter out zero values
Xag8ra_subset <- filter_if(Xag8ra_subset, is.numeric, all_vars((.) != 0))
XCV3_subset <- filter_if(XCV3_subset, is.numeric, all_vars((.) != 0))
XCV1_subset <- filter_if(XCV1_subset, is.numeric, all_vars((.) != 0))
p_capsici_subset<- filter_if(p_capsici_subset, is.numeric, all_vars((.) != 0))
TMV_subset<- filter_if(TMV_subset, is.numeric, all_vars((.) != 0))
mock_subset<- filter_if(mock_subset, is.numeric, all_vars((.) != 0))

#add conditions
condition = factor(c("mock","mock","mock","sample","sample","sample"),c("mock","sample"))
col_data = data.frame(condition)

treatments <- list(Xag8ra_subset,XCV3_subset,XCV1_subset,TMV_subset,p_capsici_subset)
names(treatments) <- c("Xag8ra_3h","XCV3_3h","XCV1_3h","TMV_4h","p_capsici_4h")
for (i in seq_along(treatments)) {
  subset <- treatments[[i]]
  subset_name <- names(treatments)[i]
  print(subset_name)
  class(subset_name)
  }
#Loop to create a DESeqDataSet data object
for (i in seq_along(treatments)) {
  subset <- treatments[[i]]
  subset_name <- names(treatments)[i]
  print(subset_name)
  dds = DESeqDataSetFromMatrix(subset, col_data, ~condition)
  #Normalization, calculate the (linear) correction factors for each sample:
  dds = estimateSizeFactors(dds)
  #The correction (size) factors can then be retrieved using:
  sizeFactors(dds)
  #t transform the counts to get a more normal distribution;
  rld = rlog(dds)
  #check data after log transformation:
  png(glue("{subset_name}_density.png"), width = 800, height = 600)
  plot(density(assay(rld)[,1]), main=glue("log counts {subset_name}"))
  dev.off()
  #To estimate the dispersions, run:
  dds = estimateDispersions(dds)
  png(glue("{subset_name}_dispersion.png"), width = 800, height = 600)
  plotDispEsts(dds)
  title(main = glue("Dispersion Estimates {subset_name}"))
  dev.off()
  dds = nbinomWaldTest(dds)
  #To get a table with differential expression values for the genes, type:
  res = results(dds)
  #For a small number of genes no padj could be calculated, these have no
  #value which will cause problems later, so set them to 1 with this command:
  res$padj = ifelse(is.na(res$padj), 1, res$padj)
  #MA plot
  png(glue("{subset_name}_MA_plot.png"), width = 800, height = 600)
  plotMA(res, main=glue("{subset_name} MA plot"),ylim=c(-8,8),alpha=0.01)
  dev.off()
  #Enhanced vulcano bacteria
  name <- glue("{subset_name}")
  png(glue("{subset_name}_vulcano.png"), width = 800, height = 600)
  print(EnhancedVolcano(res,
                  title = name,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'padj',
                  FCcutoff = 1,
                  pCutoff = 0.05))
  dev.off()
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
  upregulated_genes <- glue("{subset_name}_up")
  downregulated_genes <- glue("{subset_name}_down")
  #save all output dataframes to global environment
  assign(upregulated_genes,output_up_sorted)
  assign(downregulated_genes,output_down_sorted)
  #set filenames
  filename_up = glue("{subset_name}_upregulated.tsv")
  filename_down = glue("{subset_name}_downregulated.tsv")
  write.table(output_up_sorted, col.names=NA, row.names=T, file=as.character(filename_up), sep ="\t")
  write.table(output_down_sorted, col.names=NA, row.names=T, file =as.character(filename_down), sep="\t")
}
#check amount upregulated for each sample
treatments_up <- list(Xag8ra_3h_up,XCV3_3h_up,XCV1_3h_up,TMV_4h_up,p_capsici_4h_up)
names(treatments_up) <- c("Xag8ra_3h_up","XCV3_3h_up","XCV1_3h_up","TMV_4h_up","p_capsici_4h_up")
for (i in seq_along(treatments_up)) {
  up <- treatments_up[[i]]
  print(dim(up))
}
#check amount downregulated for each sample
treatments_down <- list(Xag8ra_3h_down,XCV3_3h_down,XCV1_3h_down,TMV_4h_down,p_capsici_4h_down)
names(treatments_up) <- c("Xag8ra_3h_down","XCV3_3h_down","XCV1_3h_down","TMV_4h_down","p_capsici_4h_down")
for (i in seq_along(treatments_down)) {
  down <- treatments_down[[i]]
  print(dim(down))
}

#prepare data for heatmap
expressiondata_nona <- na.omit(expression_data)
#remove all zero values from dataframe
expr_data_filtered<- filter_if(expressiondata_nona, is.numeric, all_vars((.) != 0))

#create dataframe with average for each treatment
average_df<-as.data.frame(rowMeans(expr_data_filtered[,c("SRR24630912","SRR24630913","SRR24630914")]))
names(average_df)[1] <-c("Xag8ra_3h")
average_df$XCVA3_3h<-rowMeans(expr_data_filtered[,c("SRR24630977","SRR24630978","SRR24630979")])
average_df$XCVA1_3h<-rowMeans(expr_data_filtered[,c("SRR24630881","SRR24630882","SRR24630883")])
average_df$TMV_4h<-rowMeans(expr_data_filtered[,c("SRR25031738","SRR25031739","SRR25031740")])
average_df$p_capscici_4h<-rowMeans(expr_data_filtered[,c("SRR24630967","SRR24630968","SRR24630969")])

#logtransform data
log_average = log(average_df)
#Compare avg of treatments
pheatmap(log_average,show_rownames=FALSE)
#clustering treatments
correlation_distance=as.dist(1-cor(log_average))
hc=hclust(correlation_distance,method="complete")
plot(hc)

#create dataframe with average for each mocktreatment
average_mock_df<-as.data.frame(rowMeans(expr_data_filtered[,c("SRR24630893","SRR24630894","SRR24630896")]))
names(average_mock_df)[1] <-c("Xag8ra_3h_mock")
average_mock_df$XCVA_3h_mock<-rowMeans(expr_data_filtered[,c("SRR24630931","SRR24630932","SRR24630933")])
average_mock_df$p_capscici_4h_mock<-rowMeans(expr_data_filtered[,c("SRR24630947","SRR24630948","SRR24630949")])
#logtransform data
log_avg_mock = log(average_mock_df)
#Compare avg of mock treatments
pheatmap(log_avg_mock,show_rownames=FALSE)
#clustering treatments
correlation_distance=as.dist(1-cor(log_avg_mock))
hc=hclust(correlation_distance,method="complete")
plot(hc)

avg_expressionset <- merge(log_avg_mock,log_average,by='row.names')
rownames(avg_expressionset) <- avg_expressionset$Row.names
avg_expressionset <- avg_expressionset[,-1]
#Compare avg and mock treatments
pheatmap(avg_expressionset,show_rownames=FALSE)


