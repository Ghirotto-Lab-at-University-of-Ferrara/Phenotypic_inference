# Phenotypic_inference

Pipeline for the inference of phenotypic traits (eye, hair and skin colours) from ancient low-coverage data.

**Requirements:**
-------------------------------------------------------------------------------------
* samtools
* bcftools
* VarScan
* Java
* R
* R libraries:
  * matrixStats
  * tibble
  * tidyr
  * dplyr
  * ggplot2
  * cowplot
  * ggrepel


**This pipeline is divided in two main parts that need to be run separately:**
-------------------------------------------------------------------------------------

**1. Convert BAM file to the input required by the hirisplex-System (function 1):**

  a) run samtools mpileup for the 41 hirisplex positions;
  
  b) quality control of the 41 positions present in the mpileup file;
  
  c) check for missing positions in the mpileup file. If there are missing positions (except the indel position, which you cannot impute) it will need to use the imputed VCF/BCF file for the sample under analysis and it will use the genotype imputed;
	
  d) check if indel position is covered. If yes, run Varscan;
	
  e) calculate genotype likelihood for the covered positions;
  
  f) convert GL/imputed genotypes to the hirisplex input format.

**2. Phenotypic prediction and plot results for the three phenotypic traits - eye, hair and skin colours - (function 2 or 3):**
	
  a1) uses as input the file generated in the hirisplex website (GL - function 2);

  a2) uses as input the file generated in the hirisplex website (direct call - function 3).

**Files required (the files available in this repository are for the hg19 genome reference. If you use a different genome reference you need to edit these files):**
-------------------------------------------------------------------------------------
- hiris_strand.list
- indel_hs37d5.list
- Allele_rs_hs37d5.csv

**Usage:**
-------------------------------------------------------------------------------------
```
./phenoPrediction.sh -ref /PATH/TO/GENOME_REF -pos /PATH/TO/positions_hirisplex-s_hs37d5.list -bam /PATH/TO/INPUT/BAM -id SAMPLE_ID -prog_func INT -impute /PATH/TO/INPUT/IMPUTED/VCF-BCF -hiris_strand_imputed /PATH/TO/hiris_strand.list -indel /PATH/TO/indel_hs37d5.list -allele_rs /PATH/TO/Allele_rs_hs37d5.csv -hirisplex_out /PATH/TO/HIRISPLEX/OUTPUT/FILE -output_folder /PATH/TO/OUTPUT/FOLDER
 
 -h,--help                           display this message
 -ref                                REFPATH: path to the reference genome
 -pos                                POSPATH: path to the list Hirisplex positions
 -bam                                BAMPATH: path to the BAM file
 -id                                 ID: Sample ID
 -prog_func                          PROGFUNC: 1-GLs Computation; 2-Prediction with GLs; 3-Prediction with direct calls
 -impute                             IMPUTEPATH: path to the VCF/BCF for imputed missing Hirisplex positions
 -hiris_strand_imputed               HIRISSTRANDPATH: path to the Hirisplex informative alleles list (check strand for Imputed VCF/BCF file)
 -indel                              INDELPATH: path to the INDEL position file (10bp upstream and downstream) for VarScan call
 -allele_rs                          ALLELERSPATH: path to the Hirisplex Informative alleles list (check strand for Rscript)
 -hirisplex_out                      OUTHIRISPATH: path to the Hirisplex output file
 -output_folder                      OUTPATH: path to the output folder
```

After run function 1, you have a csv file that can be used as input in the hirisplex website. After you upload the csv in the hirisplex website and process the data (data processing will take more or less 20 minutes), the csv file with the results will be automatically downloaded into your computer. You will use this file to run function 2 or 3 (depending on your data type) that will do the phenotypic prediction for eye, hair, and skin colours, and will give you the final results as a csv file and in pdf format with three plots.

**Example:**
-------------------------------------------------------------------------------------

**Function 1:**
```
./phenoPrediction.sh -ref ref.fasta -pos positions_hirisplex-s_hs37d5.list -bam sample.bam -id sample_ID -prog_func 1 -impute imputed.bcf  -hiris_strand_imputed hiris_strand.list -indel indel_hs37d5.list -allele_rs Allele_rs_hs37d5.csv -output_folder output
```
**Function 2:**
```
./phenoPrediction.sh -ref ref.fasta -pos positions_hirisplex-s_hs37d5.list -bam sample.bam -VCF sample.vcf.gz -id sample_ID -prog_func 2 -hiris_strand_imputed hiris_strand.list -indel indel_hs37d5.list -allele_rs Allele_rs_hs37d5.csv -output_folder output
```
**Function 3:**
```
./phenoPrediction.sh -id sample_ID -prog_func 3 -hirisplex_out Result.csv -output_folder output
```

**Function 4:**
```
./phenoPrediction.sh -id sample_ID -prog_func 4 -hirisplex_out Result.csv -output_folder output
```

**Output Example:**
-------------------------------------------------------------------------------------

```
Output GL
```


```
Output Direct Call
```
