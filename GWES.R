library(GEOquery)
library(MetaIntegrator)
library(dplyr)
library(metafor)
library(robumeta)
library(RankProd)
 #Download the data from GEO
human_gse <- getGEO("GSE209609", GSEMatrix = TRUE)
 #Extract the expression matrix
human_exprs <- exprs(human_gse[[1]])
 #Extract the phenotype data
human_pheno <- pData(phenoData(human_gse[[1]]))
 #Extract the feature data
human_feature <- fData(human_gse[[1]])
 #Get class labels from phenotype data
human_cl <- rep(1,ncol(human_exprs))
#Extract only the follwing samples from the study (GSM6380487,GSM6380488,GSM6380489,GSM6380490) from the expression matrix
human_exprs <- human_exprs[,c(1,2,3,4)]
 #Run MetaIntegrator
human_cl <- c(0,1,0,1)
#Run RankProd
human_results <- RankProducts(human_exprs, human_cl)
#Export the results to a excel file where each item in the list is a sheet
library(xlsx)
write.xlsx(human_results$RPs, "results.xlsx", sheetName = "RankProd",row.names = TRUE)

#Get new GEO data
axo_gse <- getGEO("GSE67118", GSEMatrix = TRUE)
#Extract the expression matrix
axo_exprs <- exprs(axo_gse[[1]])
#Extract the phenotype data
axo_pheno <- pData(phenoData(axo_gse[[1]]))
#Extract the feature data
axo_feature <- fData(axo_gse[[1]])


#Select from the expression matrix only the following column names (GSM1639337,GSM1639338,GSM1639339,GSM1639340 and GSM1639526,GSM1639527,GSM1639528,GSM1639529)
samples <- c("GSM1639337","GSM1639338","GSM1639339","GSM1639340","GSM1639526","GSM1639527","GSM1639528","GSM1639529")
#Use the column names to select the columns from the expression matrix
axo_exprs <- axo_exprs[,samples]
#Get the class labels
axo_cl <- c(0,0,0,0,1,1,1,1)
#Run RankProd
axo_results <- RankProducts(axo_exprs, axo_cl)
plotRP(axo_results,cutoff = 0.05)

genestats <- topGene(axo_results, cutoff = 0.05,method= "pval")

#Get Top 5 genes from Table 1
top5_under <- genestats$Table1[1:5]
top5_under[4] <- genestats$Table1[6]
top5_under[5] <- genestats$Table1[9]
probes_under <- axo_feature$ID[top5_under]
genes_under <- axo_feature$GeneSymbol[top5_under]
genes_under[2] <- "sodefrin precursor-like factor [Desmognathus monticola]"
top5_over <- genestats$Table2[1:5]
top5_over[5] <- genestats$Table2[6]
probes_over <- axo_feature$ID[top5_over]
genes_over <- axo_feature$GeneSymbol[top5_over]

data_under <- genestats$Table1[1:5,]
data_under[4] <- genestats$Table1[6,]
data_under[5] <- genestats$Table1[9,]
data_over <- genestats$Table2[1:5,]
data_over[5] <- genestats$Table2[6,]
#Combine the data into a dataframe vertically
data <- rbind(data_under,data_over)
genes <- c(genes_under,genes_over)
#Add the genes to the dataframe as a column at the beginning
data <- cbind(genes,data)
#Drop gene.id column
data <- data[,-2]
#Add a column to indicate if the gene is up or down regulated
regulation <- c("Down","Down","Down","Down","Down","Up","Up","Up","Up","Up")
data <- cbind(regulation,data)
data <- data <- data[,c(2,6,1,3,4,5)]
#Save the data to a csv file
write.csv(data,"data.csv",row.names = FALSE)

path <- file.choose()
data <- read.csv(path,header = TRUE)
#Convert data into a dataframe
data <- as.data.frame(data)
#Change name of FC column to FoldChange
names(data)[5] <- "FoldChange"
names(data)[1] <- "Gene.Symbol"

#Import corrected data
data_corrected <- read.csv("data_proc.csv",header = TRUE)
data_corrected <- as.data.frame(data_corrected)

#Generate a plot for the genes and the P-values and Fold Change
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggExtra)
library(ggcorrplot)

ggplot(data_corrected, aes(y = P.value, x = FoldChange, color = regulation)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Gene.Symbol), size = 5)+
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 1, color = "orange") +
  annotate("text", x = 5.5, y = 0.0000000001, label = "Upregulated", color = "blue") +
  annotate("text", x = 0.05, y = 0.0000000001, label = "Downregulated", color = "red") +
  ggtitle("Most Significant Genes")

#Save the plot
ggsave("plot.png")
