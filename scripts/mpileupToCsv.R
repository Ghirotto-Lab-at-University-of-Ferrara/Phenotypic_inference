library(matrixStats)
library(tibble)
library(tidyr)
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)

sampleid<-as.character(args[1])
mpileup_file<-as.character(args[2])
allele_rs<-read.csv(as.character(args[3]),header=T) #Allele_rs_hs37d5.csv
geno_indel<-as.character(args[4])
geno_imputed<-as.character(args[5])

options(digits = 6)

############ FORMAT MPILEUP FILE ############
mpileup <- read.table(mpileup_file, sep="\t", header=FALSE, stringsAsFactors=FALSE, col.names=c("chr", "pos", "ref", "depth", "bases", "qualities_ASCII"))

# REMOVE SYMBOLS 
for (i in 1:nrow(mpileup)) {
  mpileup$bases[i] <- gsub("\\$", "", mpileup$bases[i])
  mpileup$bases[i] <- gsub("\\^[[:graph:]]", "", mpileup$bases[i])
  mpileup$bases[i] <- gsub("\\<|\\>", "", mpileup$bases[i])
  mpileup$bases[i] <- gsub("\\.", mpileup$ref[i], mpileup$bases[i])
  mpileup$bases[i] <- gsub("\\,", mpileup$ref[i], mpileup$bases[i])
  mpileup$bases[i] <- sapply(mpileup$bases[i], function(x) toupper(x))
}

# REMOVE MISSING BASES
mpileup <- mpileup[mpileup$bases != "*", ]

if(nrow(mpileup)!=0){

	# REMOVE BASES WITH QUALITY BELLOW 30
	mpileup$bases <- sapply(1:nrow(mpileup), function(i) {
	  qualities <- strsplit(mpileup[i, "qualities_ASCII"], "")[[1]]
	  low_quality_indices <- which((as.integer(charToRaw(paste(qualities, collapse = ""))) - 33) < 30)
	  chars <- strsplit(mpileup[i, "bases"], "")[[1]]
	  filtered_chars <- chars[!(seq_along(chars) %in% low_quality_indices)]
	  paste(filtered_chars, collapse = "")
	})

	mpileup$qualities_ASCII <- sapply(1:nrow(mpileup), function(i) {
	  qualities <- strsplit(mpileup[i, "qualities_ASCII"], "")[[1]]
	  low_quality_indices <- which((as.integer(charToRaw(paste(qualities, collapse = ""))) - 33) < 30)
	  chars <- strsplit(mpileup[i, "qualities_ASCII"], "")[[1]]
	  filtered_chars <- chars[!(seq_along(chars) %in% low_quality_indices)]
	  paste(filtered_chars, collapse = "")
	})

	# REMOVE BASES WITH QUALITY BELLOW 30 (CHECK)
	mpileup <- subset(mpileup, bases != "")

	#COMPUTING DEPTH
	mpileup$depth <- sapply(mpileup$bases, nchar)

	#CONVERSION ASCII SYMBOLS INTO NUMBERS
	mpileup$qualities <- lapply(mpileup$qualities_ASCII, function(x) as.numeric(charToRaw(iconv(x, to="UTF-8")))-33)

	#COMPUTE LOG ERROR RATE
	mpileup$log10_err <- lapply(mpileup$qualities, function(x) -x/10)

	#FORMAT BASES
	mpileup$bases <- sapply(strsplit(as.character(mpileup$bases), ""), function(x) paste(x, collapse = ","))


	############ COMPUTE GENOTYPE LIKELIHOODS ############

	P_BmatchA <- function(log10_err) {
	  risultato <- log10(1 - 10^log10_err)
	  return(risultato)
	}

	P_BmismatchA <- function(log10_err) {
	  risultato <- log10_err - log10(3)
	  return(risultato)
	}

	logLike_match_match <- function(log10_err) {
	  logP_BmatchA <- P_BmatchA(log10_err)
	  risultato <- matrixStats::logSumExp(c(logP_BmatchA))
	  return(risultato)
	}

	logLike_match_mismatch <- function(log10_err) {
	  logP_BmatchA <- P_BmatchA(log10_err) - log10(2)
	  logP_BmismatchA <- P_BmismatchA(log10_err) - log10(2)
	  risultato <- matrixStats::logSumExp(c(logP_BmatchA, logP_BmismatchA)) 
	  return(risultato)
	}

	logLike_mismatch_mismatch <- function(log10_err) {
	  logP_BmismatchA <- P_BmismatchA(log10_err)
	  risultato <- matrixStats::logSumExp(c(logP_BmismatchA))
	  return(risultato)
	}


	# CREATE GENOTYPES LIKELIHOOD VECTORS
	genotypes <- c("AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT")
	risultati_list <- list()

	# COMPUTING LOG ENOTYPE LIKELIHOOD
	for (i in 1:nrow(mpileup)) {
	  bases <- unlist(strsplit(mpileup$bases[i], ","))
	  log10_err <- unlist(mpileup$log10_err[i])
	  qualities <- unlist(mpileup$qualities[i])
	  
	  risultati <- matrix(nrow = length(bases), ncol = length(genotypes), dimnames = list(bases, genotypes))
	  
	  for (j in 1:length(genotypes)) {
		allele <- unlist(strsplit(genotypes[j], ""))
		A1 <- allele[1]
		A2 <- allele[2]
		
		for (k in 1:length(bases)) {
		  if ((bases[k] == A1) && (bases[k] == A2)) {
			risultati[k, j] <- logLike_match_match(log10_err[k])
		  } else if ((bases[k] == A1) && (bases[k] != A2)) {
			risultati[k, j] <- logLike_match_mismatch(log10_err[k])
		  } else if ((bases[k] != A1) && (bases[k] == A2)) {
			risultati[k, j] <- logLike_match_mismatch(log10_err[k])
		  } else {
			risultati[k, j] <- logLike_mismatch_mismatch(log10_err[k])
		  }
		}
		risultati_list[[i]] <- risultati
	  }
	}

	# CREATE NEW TABLE WITH THE 10 POSSIBLE GENOTYPES
	genLogLikelihood <- as.data.frame(matrix(NA, nrow = nrow(mpileup), ncol = length(genotypes)))
	colnames(genLogLikelihood) <- genotypes
	rownames(genLogLikelihood) <- paste(mpileup$chr, mpileup$pos, sep = ":")

	# NORMALIZE LOG GENOTYPE LIKELIHOOD VALUES
	for (i in 1:nrow(genLogLikelihood)) {
	  for (j in 1:ncol(genLogLikelihood)) {
		genLogLikelihood[i, j] <- sum(risultati_list[[i]][, j])
	  }
	}
	genLogLikelihood_normalized <- t(apply(genLogLikelihood, 1, function(row) {
	  sumExp <- matrixStats::logSumExp(row, na.rm = TRUE)
	  row - sumExp
	}))

	# WRITE FILE genLogLikelihood
	write.table(genLogLikelihood_normalized, paste(sampleid,"_genLogLikelihood.txt",sep=""), sep = "\t")


	############ COMPUTING GENOTYPES POSTERIOR WITH PRIOR HET 1/1000 ############

	# CREATE EMPTY DATAFRAME
	genPosteriors <- genLogLikelihood_normalized
	col_names <- colnames(genLogLikelihood_normalized)

	# PRIOR FOR HOM AND HET GENOTYPES
	prior_hom <- log10(0.999)
	prior_het <- log10(0.001)

	# COMPUTE GENOTYPE POSTERIORS
	for (col_name in col_names) {
	  genotype <- substring(col_name, first = 1, last = 2)
	  if (genotype %in% c("AA", "CC", "GG", "TT")) {
		genPosteriors[, col_name] <- genPosteriors[, col_name] + prior_hom
	  } else {
		genPosteriors[, col_name] <- genPosteriors[1,1] + prior_het
	  }
	}

	# NORMALIZE GENOTYPE POSTERIORS
	for (i in 1:nrow(genPosteriors)) {
	  row <- genPosteriors[i, ]
	  sumExp <- matrixStats::logSumExp(row)
	  row <- row - sumExp
	  row <- 10^row
	  row <- row / sum(row)
	  genPosteriors[i, ] <- row
	}

	# WRITE FILE genPosteriors
	write.table(genPosteriors, paste(sampleid,"_genPosteriors.txt",sep=""), sep = "\t") 
	  

	############ SAMPLING 1000 GENOTYPES FROM POSTERIORS AND COUNT HIRISPLEX-S ALLELES ############

	genPosteriors <- as.data.frame(genPosteriors)
	genPosteriors <- rownames_to_column(genPosteriors, var = "Position")
	Hps_SNPs<-allele_rs
	genPosteriors <- merge(genPosteriors, Hps_SNPs, by = "Position", all = TRUE)
	genPosteriors <- genPosteriors[order(genPosteriors$Num), ]

	sampling_list <- list()

	for (j in 1:1000) {
	  sampling <- data.frame(matrix(ncol = 1, nrow = nrow(genPosteriors)))
	  row.names(sampling) <- rownames(genPosteriors)
	  colnames(sampling) <- "Sampling Genotype"
	  colnames(sampling) <- "Sampled"
	  sampling$Num <- genPosteriors$Num
	  sampling$Rs_allele <- genPosteriors$Rs_allele
	  sampling$Strand <- genPosteriors$Strand
	  sampling$Position <- genPosteriors$Position
	  for (i in 1:nrow(genPosteriors)) {
		allele_interesse <- genPosteriors[i, "Strand"]
		probabilita <- genPosteriors[i, 2:11]
		probabilita[is.na(probabilita)] <- 0 
		if (all(probabilita == 0)) {
		  genotipo_i <- NA
		} else {
		  genotipo_i <- sample(names(genPosteriors)[2:11], 1, prob = probabilita)
		}
		sampling[i, "Sampling Genotype"] <- genotipo_i
		if (is.na(genotipo_i)) {
		  sampling[i, "Sampled"] <- NA
		} else if (genotipo_i == paste0(allele_interesse, allele_interesse)) {
		  sampling[i, "Sampled"] <- 2
		} else if (grepl(allele_interesse, genotipo_i)) {
		  sampling[i, "Sampled"] <- 1
		} else {
		  sampling[i, "Sampled"] <- 0
		}
	  }
	  sampling_list[[j]] <- sampling
	  sampling_list[[j]][is.na(sampling_list[[j]])] <- "NA"
	}


	for (j in 1:length(sampling_list)) {
	  names(sampling_list[[j]]) <- paste(names(sampling_list[[j]]), j, sep="_")
	}

	# CREATE DATAFRAME
	sampling_df_check <- do.call(cbind.data.frame, sampling_list)
	sampled_column <- grep("^Sampled", names(sampling_df_check), value = TRUE)
	sampling_df <- sampling_df_check[, sampled_column]

	# FORMAT AND CREATE CSV FILE TO USE AS INPUT IN THE HIriplex-S WEBSITE
	sampling_df <- data.frame(t(sampling_df))
	colnames(sampling_df) = Hps_SNPs$Rs_allele

	###### INSERT INDEL GENOTYPE (0/1/2/NA) ######
	sampling_df$rs312262906_A <- geno_indel

	###### INSERT IMPUTED GENOTYPES (0/1/2/NA) ######
	if (geno_imputed != 0) {
	  geno_imputed<-read.table(geno_imputed)
	  colnames(geno_imputed)<-c("Chr","Position","Ref","Alt","GT","Rs_allele")
	  rs_imputed<-c()
	  for (i in 1:length(as.vector(geno_imputed[,"Rs_allele"]))){
	  rs_imputed<-as.vector(geno_imputed[i,"Rs_allele"])
	  sampling_df[,rs_imputed]<-geno_imputed[i,"GT"]
	  }
	}

	row.names(sampling_df) <- paste(sampleid, row.names(sampling_df), sep = "_")
	output_df <- rownames_to_column(sampling_df)
	colnames(output_df)[1] = "sampleid"

	# WRITE CSV FILE FOR HIRISPLEX WEBSITE
	write.csv(output_df, file = paste(sampleid,"_1000sampling.csv",sep=""), row.names = FALSE, quote=FALSE)

} else {

	###### CREATE CSV FOR HIRISPLEX
 	#Hps_SNPs<-allele_rs
	output_df <- data.frame(matrix(ncol = 42, nrow = 1))
	colnames(output_df)[1] = "sampleid"
	colnames(output_df)[2:42]<-allele_rs$Rs_allele
	output_df$sampleid <- sampleid

	###### INSERT INDEL GENOTYPE (0/1/2/NA) ######
	output_df$rs312262906_A <- geno_indel

	###### INSERT IMPUTED GENOTYPES (0/1/2/NA) ######
	if (geno_imputed != 0) {
	  geno_imputed<-read.table(geno_imputed)
	  colnames(geno_imputed)<-c("Chr","Position","Ref","Alt","GT","Rs_allele")
	  rs_imputed<-c()
	  for (i in 1:length(as.vector(geno_imputed[,"Rs_allele"]))){
	  rs_imputed<-as.vector(geno_imputed[i,"Rs_allele"])
	  output_df[,rs_imputed]<-geno_imputed[i,"GT"]
	  }
	}
	write.csv(output_df, file = paste(sampleid,"_onlyimputedcalls.csv",sep=""), row.names = FALSE, quote=FALSE)
}

