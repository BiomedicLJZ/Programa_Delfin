#START OF SCRIPT 2
library("biomaRt")
library("tidyr")
library("dbplyr")
library("dplyr")
DF <- read.table("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90078639_sig_P_rsID.tsv",sep="\t",header=TRUE)
mart.snp <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp")
getENSG <- function(rs = "", mart = mart.snp) {
  results <- getBM(attributes = c("refsnp_id","ensembl_gene_stable_id"),
                   filters    = "snp_filter",
                   values = rs,
                   mart = mart)
  return(results)
}
ensg <- data.frame()
for (i in 1:100) {
  rs <- DF[i,1]
  temp <- getENSG(rs)
  ensg <- rbind(ensg,temp)
}
#DF$ensembl_gene_id <- ensg$ensembl_gene_stable_id
#write.table(DF,file="D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90083793_sig_P_rsID_ENSG.tsv",sep="\t",row.names = FALSE,col.names = TRUE)
ensg[ensg==""]<-NA
ensg %>% drop_na(ensembl_gene_stable_id) ->ensg
ensg <- unique(ensg)
names(ensg)[1]<-"RefSNP_id"
final <- merge(DF,ensg,by="RefSNP_id")
dir.create("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90078639")
write.table(final,"D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90078639/GCST90078639_ENSEMBL.tsv",sep='\t')