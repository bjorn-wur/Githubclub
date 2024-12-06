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

#load packages
library(DESeq2)
library(pheatmap)
library(dplyr)
#generate PCA
library("FactoMineR")
library("factoextra")
#vulcanoplot
library(EnhancedVolcano)
#load dataset
expression_data=read.table("gene_count_matrix.csv", row.names=1,
  header=TRUE, sep =",", stringsAsFactors=FALSE)


#subset for the bacteria
Xag8ra_treatment<- c("SRR24630900","SRR24630901","SRR24630902","SRR24630918","SRR24630919","SRR24630920")
#subset for the virus
virus_treatment<- c("SRR25031747","SRR25031741","SRR25031748")
#subset for the fungus
fungus_treatment<- c("SRR24630954","SRR24630955","SRR24630956","SRR24630895","SRR24630906","SRR24630917")

#subset mocksamples
mock_treatment <- c("SRR24630954","SRR24630955","SRR24630956","SRR24630900","SRR24630901","SRR24630902")

#subset data
bacteria_subset <- expression_data[Xag8ra_treatment]
fungus_subset <- expression_data[fungus_treatment]
mock_subset <- expression_data[mock_treatment]

colnames(bacteria_subset)<- c("mock_1","mock_2","mock_3","Xag8ra_1","Xag8ra_2","Xag8ra_3")
colnames(fungus_subset)<- c("mock_1","mock_2","mock_3","P_capsici_1","P_capsici_2","P_capsici_3")
colnames(mock_subset)<- c("mock_bacteria_1","mock_bacteria_2","mock_bacteria_3","mock_fungus_1","mock_fungus_2","mock_fungus_3")
#check data
dim(bacteria_subset)
colnames(bacteria_subset)
rownames(bacteria_subset)

#remove all na values
bacteria_subset <- na.omit(bacteria_subset)
fungus_subset <- na.omit(fungus_subset)
mock_subset<- na.omit(mock_subset)
#remove genes that have in none of the samples more than 10 counts
mx = apply(bacteria_subset, 1, max)
mx = apply(fungus_subset, 1, max)

# Next, we will make a new table that only contains rows for which the maximum count is greater than 10
bacteria_subset = bacteria_subset[mx>0,]
fungus_subset = fungus_subset[mx>0,]
#filter-out zero values
bacteria_subset <- filter_if(bacteria_subset, is.numeric, all_vars((.) != 0))
fungus_subset <- filter_if(fungus_subset, is.numeric, all_vars((.) != 0))
mock_subset<- filter_if(mock_subset, is.numeric, all_vars((.) != 0))

dim(bacteria_subset)
dim(fungus_subset)
dim(mock_subset)
#Run these commands to create a DESeqDataSet data object
condition = factor(c("mock","mock","mock","sample","sample","sample"),c("mock","sample"))
col_data = data.frame(condition)
#add conditions
dds_b = DESeqDataSetFromMatrix(bacteria_subset, col_data, ~condition)
dds_f = DESeqDataSetFromMatrix(fungus_subset, col_data, ~condition)


summary(bacteria_subset)
summary(fungus_subset)

#Generate heatmap of dataset
pheatmap(bacteria_subset,show_rownames=FALSE)
pheatmap(fungus_subset,show_rownames=FALSE)

#logtransform
bacteria_log = log(bacteria_subset)
fungus_log = log(fungus_subset)
mock_log = log(mock_subset)
#Generate heatmap of log transformed dataset
pheatmap(bacteria_log,show_rownames=FALSE)
pheatmap(fungus_log,show_rownames=FALSE)

#meanscale & variance normalize data
bacteria_scaled = scale(bacteria_log,center=TRUE,scale=TRUE)
fungus_scaled = scale(fungus_log,center=TRUE,scale=TRUE)

#Generate heatmap of scaled dataset
pheatmap(bacteria_scaled,show_rownames=FALSE,main = "Xag8ra")
pheatmap(fungus_scaled,show_rownames=FALSE,main ="P.capsici")

#hierarchical clustering, pearson correlation distance
#correlation based distance matrix
#pearson correlation based distance, complete linkage (calculates max distance)
#Cluster samples bacteria:
correlation_distance=as.dist(1-cor(bacteria_scaled))
hc=hclust(correlation_distance,method="complete")
plot(hc)

#Cluster samples fungus:
correlation_distance=as.dist(1-cor(fungus_scaled))
hc=hclust(correlation_distance,method="complete")
plot(hc)

correlation_distance=as.dist(1-cor(mock_log))
hc=hclust(correlation_distance,method="complete")
plot(hc)

kmeans(2, centers, iter.max = 10, nstart = 1)

#Normalization, calculate the (linear) correction factors for each sample:
dds_b = estimateSizeFactors(dds_b)
dds_f = estimateSizeFactors(dds_f)
#check data after normalization:
#norm_versus_non_norm( dds, 2, 4, left = 2, right = 8 )
#The correction (size) factors can then be retrieved using:
sizeFactors(dds_b)
sizeFactors(dds_f)
#t transform the counts to get a more normal distribution;
rld_b = rlog(dds_b)
rld_f = rlog(dds_f)
#check data after log transformation:
plot(density(assay(rld_b)[,1]), main="log counts")
plot(density(assay(rld_f)[,1]), main="log counts")
#To estimate the dispersions, run:
dds_b = estimateDispersions(dds_b)
dds_f = estimateDispersions(dds_f)
plotDispEsts(dds_b)
plotDispEsts(dds_f)
dds_b = nbinomWaldTest(dds_b)
dds_f = nbinomWaldTest(dds_f)
#To get a table with differential expression values for the genes, type:
res_b = results(dds_b)
res_f = results(dds_f)
#To see the top rows from the differential expression table, type the following command in
#the console:
head(res_b)
head(res_f)
#For a small number of genes no padj could be calculated, these have no
#value which will cause problems later, so set them to 1 with this command:
res_b$padj = ifelse(is.na(res_b$padj), 1, res_b$padj)
res_f$padj = ifelse(is.na(res_f$padj), 1, res_f$padj)
#To end this exercise, we will make an MA plot. You can make the plot by typing
plotMA(res_b, main="MA plot",ylim=c(-8,8),alpha=0.01)
plotMA(res_f, main="MA plot",ylim=c(-8,8),alpha=0.01)

#Enhanced vulcano bacteria
EnhancedVolcano(res_b,
                lab = rownames(res_b),
                title = 'Xag8ra',
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05)
head(res_b)

#Enhanced vulcano fungus
EnhancedVolcano(res_f,
                lab = rownames(res_f),
                title = 'P.capsici',
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05)
head(res_f)

#turn res into a dataframe
output_table_b <- as.data.frame(res_b)
output_table_f <- as.data.frame(res_f)

#filter for significant genes <0.05
output_filtered_b <- output_table_b[output_table_b$padj<0.05,]
output_filtered_f <- output_table_f[output_table_f$padj<0.05,]
#filter for fold change
output_filtered_up_b <- output_filtered_b[output_filtered_b$log2FoldChange>1,]
output_filtered_down_b <-output_filtered_b[output_filtered_b$log2FoldChange < -1,]
output_filtered_up_f <- output_filtered_f[output_filtered_f$log2FoldChange>1,]
output_filtered_down_f <-output_filtered_f[output_filtered_f$log2FoldChange < -1,]

#order data of upregulated genes by adjusted p-value
output_up_f_sorted <- output_filtered_up_f[order(output_filtered_up_f$padj),]
output_up_b_sorted <- output_filtered_up_b[order(output_filtered_up_b$padj),]

#order data of upregulated genes by adjusted p-value
output_down_f_sorted <- output_filtered_down_f[order(output_filtered_down_f$padj),]
output_down_b_sorted <- output_filtered_down_b[order(output_filtered_down_b$padj),]

#Combine up/down in the same dataframe
diff_expr_f <- rbind(output_up_f_sorted,output_down_f_sorted)
diff_expr_b <- rbind(output_up_b_sorted,output_down_b_sorted)

#


dim(output_filtered_up_b)
dim(output_filtered_down_b)

dim(output_filtered_up_f)
dim(output_filtered_down_f)

#write bacteria result to tsv file
write.table(output_filtered_up_b, col.names=NA, row.names=T, file ="bacteria_upregulated.tsv", sep
            ="\t")
write.table(output_filtered_down_b, col.names=NA, row.names=T, file ="bacteria_downregulated.tsv", sep
            ="\t")
#write fungus result to tsv file
write.table(output_filtered_up_f, col.names=NA, row.names=T, file ="fungus_upregulated.tsv", sep
            ="\t")
write.table(output_filtered_down_f, col.names=NA, row.names=T, file ="fungus_downregulated.tsv", sep
            ="\t")


#PCA
pca_exp_data <- PCA(scaled_expressiondata)
plot.PCA(pca_exp_data,choix=c("ind"))
fviz_pca_ind(pca_exp_data, label="none")
fviz_screeplot(pca_exp_data, ncp=10, addlabels='TRUE',xlab='principal components',main='scree plot')
#habbilage




