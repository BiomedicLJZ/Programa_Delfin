#Libraries
library(biomaRt)
library(annotables)
library(org.Hs.eg.db)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

DF <- read.csv("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90078639/GCST90078639_ENSEMBL.tsv",sep = '\t')
Genes <- DF$ensembl_gene_stable_id
listEnsembl()
listMarts()
ensembl <- useEnsembl(biomart = "genes")
list <- listDatasets(ensembl)
ensembl.con <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl.con)
filter <- listFilters(ensembl.con)
DF2 <- getBM(attributes = c("hgnc_id","external_gene_name", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values= Genes,
      mart = ensembl.con)

#DF <- read.csv("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90083797_sig_P_rsID.tsv",sep = '\t')
names(DF)[10] <-"ensembl_gene_id"
DFF <- merge(DF,DF2,by="ensembl_gene_id",)
write.table(DFF,file="D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90078639_GENES.tsv",sep ='\t',col.names=TRUE,row.names=FALSE)
