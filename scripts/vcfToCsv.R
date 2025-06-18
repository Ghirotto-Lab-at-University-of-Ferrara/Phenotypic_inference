library(matrixStats)
library(tibble)
library(tidyr)
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)

sampleid<-as.character(args[1])
allele_rs<-read.csv(as.character(args[2]),header=T) #Allele_rs_hs37d5.csv
geno_indel<-as.character(args[3])
geno_vcf<-as.character(args[4])

options(digits = 6)

###### CREATE CSV FOR HIRISPLEX
Hps_SNPs<-allele_rs
final_csv_df <- data.frame(matrix(ncol = 42, nrow = 1))
colnames(final_csv_df)[1] = "sampleid"
colnames(final_csv_df)[2:42]<-Hps_SNPs$Rs_allele
final_csv_df$sampleid <- sampleid

###### INSERT INDEL GENOTYPE (0/1/2/NA) ######
final_csv_df$rs312262906_A <- geno_indel

###### INSERT DIRECT CALL GENOTYPES (0/1/2/NA) ######
if (geno_vcf != 0) {
  geno_vcf<-read.table(geno_vcf)
  colnames(geno_vcf)<-c("Chr","Position","Ref","Alt","GT","Rs_allele")
  rs_vcf<-c()
  for (i in 1:length(as.vector(geno_vcf[,"Rs_allele"]))){
  rs_vcf<-as.vector(geno_vcf[i,"Rs_allele"])
  final_csv_df[,rs_vcf]<-geno_vcf[i,"GT"]
  }
}


###### WRITE CSV FILE FOR HIRISPLEX WEBSITE
write.csv(final_csv_df, file = paste(sampleid,"_directCall.csv",sep=""), row.names = FALSE, quote=FALSE)
