
library("BSgenome")
setwd("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed")
library("SNPlocs.Hsapiens.dbSNP144.GRCh38")
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
#chrsnps <- snpsBySeqname(snps, "1")
#create a cicle for each chromosome
dir.create("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90083797")
dir.create("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90083797/rsID")
for (i in 1:22) {
  chr <- i
  seq <- toString(chr)
  sprintf("Processing chromosome: %s", seq)
  chrsnps <- snpsBySeqname(snps,seq)
  path <- paste("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/GCST90083797/GCST90083797_chr",seq,".tsv",sep="")
  gwas <- read.table(path,sep="\t",header=TRUE)
  gr <- makeGRangesFromDataFrame(
  gwas,
  seqnames.field="chromosome",
  start.field="base_pair_location",
  end.field="base_pair_location",
  keep.extra.columns = TRUE,
  na.rm = TRUE
)
qHits <- queryHits(findOverlaps(query = gr, subject = chrsnps,type="within"))
subHits <- subjectHits(findOverlaps(query = gr, subject = chrsnps,type="within"))
overlaps <- data.frame(gwas[qHits,],chrsnps[subHits,])
head(overlaps)
destination <- paste("D:/Documentos/CUCEI/Delfin/GWAS/Anxiety/Processed/GCST90083797/rsID/GCST90083797_chr",seq,"_rsID.tsv",sep="")
write.table(overlaps,file=destination,sep="\t",row.names = FALSE,col.names = TRUE)
sprintf("Chromosome %s finished", seq)
rm(chrsnps,gr,qHits,subHits,overlaps,gwas)
}
#write table to file
print("Finished")
#END OF SCRIPT 1
