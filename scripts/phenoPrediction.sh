#!/bin/bash

######################## ARGS ########################

ref=${1}                        #Reference genome
pos=${2}                        #list Hirisplex positions
bam=${3}                        #BAM file
VCF=${4}                        #VCF file
id=${5}                         #Sample ID
prog_func=${6}                  #Program Function: 1-GLs Computation; 2-Prediction with GLs; 3-Prediction with direct calls
impute=${7}                     #VCF/BCF for imputed missing Hirisplex positions
hiris_strand_imputed=${8}       #Informative alleles Hirisplex (check strand for Imputed VCF/BCF file)          
indel=${9}                      #INDEL position with 10bp upstream and downstream for VarScan call
allele_rs=${10}                  #Informative alleles Hirisplex (check strand for Rscript)
hirisplex_out=${11}             #Hirisplex output file
plot=${12}                      #Create plot for phenotype prediction
output_folder=${13}             #Path to output folder

######################################################



######################## MODULES #####################

source /opt/miniconda3/bin/activate r-env
module load samtools-1.11
module load bcftools-1.11
varscan=/jarvis/scratch/usr/perretti/software/VarScan.v2.3.9.jar
scriptsWD=/jarvis/scratch/usr/santos/pheno_lowCov

######################################################

################ VARIABLES' ASSIGNMENT ###############

while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo " "
                        echo "Command Line Arguments:"
                        echo " "
                        echo "Example: ./global_test.sh -ref /PATH/TO/GENOME_REF -pos /PATH/TO/positions_hirisplex-s_hs37d5.list -bam /PATH/TO/INPUT/BAM -VCF /PATH/TO/INPUT/VCF -id SAMPLE_ID -prog_func INT -impute /PATH/TO/INPUT/IMPUTED/VCF-BCF -hiris_strand_imputed /PATH/TO/hiris_strand.list -indel /PATH/TO/indel_hs37d5.list -allele_rs /PATH/TO/Allele_rs_hs37d5.csv -hirisplex_out /PATH/TO/HIRISPLEX/OUTPUT/FILE -output_folder /PATH/TO/OUTPUT/FOLDER"
                        echo " "                 
                        echo " -h,--help                           Display this message"
                        echo " -ref                                REFPATH: path to the reference genome"
                        echo " -pos                                POSPATH: path to the list Hirisplex positions"
                        echo " -bam                                BAMPATH: path to the BAM file"
                        echo " -VCF                                VCFPATH: path to the VCF file"
                        echo " -id                                 ID: Sample ID"
                        echo " -prog_func                          PROGFUNC: 1-GLs Computation; 2-Prediction with GLs; 3-Prediction with direct calls"
                        echo " -impute                             IMPUTEPATH: path to the VCF/BCF for imputed missing Hirisplex positions"
                        echo " -hiris_strand_imputed               HIRISSTRANDPATH: path to the Hirisplex informative alleles list (check strand for Imputed VCF/BCF file)"
                        echo " -indel                              INDELPATH: path to the INDEL position file (10bp upstream and downstream) for VarScan call"
                        echo " -allele_rs                          ALLELERSPATH: path to the Hirisplex Informative alleles list (check strand for Rscript)"
                        echo " -hirisplex_out                      OUTHIRISPATH: path to the Hirisplex output file"
                        echo " -output_folder                      OUTPATH: path to the output folder"
                        echo " "
                        exit 0
                        ;;
                -ref)
                        shift
                        if test $# -gt 0; then
                                 ref=${1}
                        else
                                echo "[PREDICTION]---[NO REFPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;
                -pos)
                        shift
                        if test $# -gt 0; then
                                 pos=${1}
                        else
                                echo "[PREDICTION]---[NO POSPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;
                -bam)
                        shift
                        if test $# -gt 0; then
                                 bam=${1}
                        else
                                echo "[PREDICTION]---[NO BAMPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;
                -VCF)
                        shift
                        if test $# -gt 0; then
                                 VCF=${1}
                        else
                                echo "[PREDICTION]---[NO VCFPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;
                -id)
                        shift
                        if test $# -gt 0; then
                                 id=${1}
                        else
                                echo "[PREDICTION]---[NO ID SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;    
                -prog_func)
                        shift
                        if test $# -gt 0; then
                                 prog_func=${1}
                        else
                                echo "[PREDICTION]---[NO PROGFUNC SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;  
                -impute)
                        shift
                        if test $# -gt 0; then
                                 impute=${1}
                        else
                                echo "[PREDICTION]---[NO IMPUTEPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;   
                -hiris_strand_imputed)
                        shift
                        if test $# -gt 0; then
                                 hiris_strand_imputed=${1}
                        else
                                echo "[PREDICTION]---[NO HIRISSTRANDPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;   
                -indel)
                        shift
                        if test $# -gt 0; then
                                 indel=${1}
                        else
                                echo "[PREDICTION]---[NO INDELPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;            
                -allele_rs)
                        shift
                        if test $# -gt 0; then
                                 allele_rs=${1}
                        else
                                echo "[PREDICTION]---[NO ALLELERSPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;  
                -hirisplex_out)
                        shift
                        if test $# -gt 0; then
                                 hirisplex_out=${1}
                        else
                                echo "[PREDICTION]---[NO OUTHIRISPATH SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;  
                -output_folder)
                        shift
                        if test $# -gt 0; then
                                 output_folder=${1}
                        else
                                echo "[PREDICTION]---[NO OUTPUT FOLDER SPECIFIED]"
                                exit 1
                        fi
                        shift
                        ;;  
        esac
done

######################################################

echo -e "\n--------------------------> EXECUTING CMD LINE: ./global_test.sh -ref $ref -pos $pos -bam $bam -VCF $VCF -id $id -prog_func $prog_func -impute $impute -hiris_strand_imputed $hiris_strand_imputed -indel $indel -allele_rs $allele_rs -hirisplex_out $hirisplex_out -output_folder $output_folder \n"

#################### PROG FUNC 1 #####################

if [[ $prog_func == "1" ]]; then
  echo -e "--------------------------> Running program function $prog_func: \n"
  echo -e "--------------------------> Compute Genotype Likelihood and Imputation for missing positions (if present) for sample: $id \n"

  # Create ouput ID folder
  mkdir -p $output_folder/$id
  cd $output_folder/$id
  
  # Create mpileup file
  samtools mpileup -f $ref -l $pos -a -B -Q 30 $bam -o $id.mpileup
  
  # Create report QC positions
  file=$id.mpileup
  
  echo -e "--------------------------" > QC_report.log
  echo -e "report QC positions" >> QC_report.log
  
  # Check for positions with DP==0
  if awk '$4 == 0' "$file" | grep -q '.'; then
  
    echo -e "--> Missing positions:"  >> QC_report.log
    awk '$4 == 0 {print "chr="$1, "pos="$2, "DP="$4; count++} END {print "Positions count:", count}' "$file" >> QC_report.log
  
  else
    echo -e "--> Missing positions: 0" >> QC_report.log
    awk '$4 == 0 {print "chr="$1, "pos="$2, "DP="$4; count++} END {print "Positions count: 0"}' "$file" >> QC_report.log
  fi
  # Compute minimum DP for covered positions
    echo -n "--> Minimum coverage: " >> QC_report.log
    min=$(awk '$4 != 0 && (min=="" || $4 < min) {min=$4} END {print min}' "$file") >> QC_report.log
    echo $min >> QC_report.log
  
  # Report for minimum DP covered positions
  if [ -n "$min" ]; then
    awk -v min="$min" '$4 == min {
      print "chr="$1, "pos="$2, "DP="$4 
      }' "$file" | LC_ALL=C sort -t= -k2,2n >> QC_report.log
    echo -e "Total Positions with minimum DP: $(awk -v min="$min" '$4 == min {count++} END {print count}' "$file")" >> QC_report.log
  fi
  
    echo -e "--------------------------" >> QC_report.log

  # Check for missing position to be imputed and INDEL status
  QC_file=$output_folder/$id/QC_report.log
  missing_pos=$(cat $QC_file | grep 'Positions count' | awk -F': ' '{print $2}')
  indel_pos=$(cat $file | grep -wE "16[[:space:]]89985753" | awk '{print $4}')
  
  if [[ $missing_pos -gt 1 ]]; then
    if [[ $impute == "0" ]]; then 
      echo -e "\n--------------------------> Missing positions detected, check QC report! \n"
      exit 1
    else 
      echo -e "\n--------------------------> Imputation file detected! \n"
      echo -e "--------------------------> Check missing positions on the QC report and extract corresponding GT from the imputed VCF/BCF \n"
      cat $QC_file | grep "DP=0" | awk -F'=' '{print $2,$3}' | awk -F' ' '{print $1"\t"$3}' > missing_pos_toImpute.list
      bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n' -R missing_pos_toImpute.list $impute  > imputed_genoMissing_noStrand.geno
      awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print $3,$4, a[$1$2]}' imputed_genoMissing_noStrand.geno $hiris_strand_imputed | awk '{print $3,$4,$5,$6,$7,$1,$2}' | awk '{if ($4 == $7) {print $5, $6} else if ($3 == $7) {if ($5 == "0|0") {print "1|1", $6} else if ($5 == "1|1") {print "0|0", $6} else {print $5, $6}} else {print $5, $6}}' | paste imputed_genoMissing_noStrand.geno - | awk '{print $1, $2, $3, $4, $6, $7}' | sed -e "s/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g" > imputed_genoMissing.geno
      echo -e "--------------------------> Check if INDEL position is covered in the sample BAM file \n"
      if [[ $indel_pos -eq 0 ]]; then
        echo -e "--------------------------> INDEL position is missing, run Pipeline with INDEL as missing position: NA! \n"
        geno_indel="NA"
      else 
        echo -e "--------------------------> INDEL position is covered, run VarScan and Pipeline! \n"
        samtools mpileup -f $ref -l $indel -B $bam | java -Xmx20g -jar $varscan mpileup2indel --min-coverage 2 --min-reads2 2 --output-vcf > "$id"_VarScan.vcf
        ref_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $5}' )
        alt_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $6}' )
        if [[ $ref_DP_indel -eq '' ]]; then 
          geno_indel=0
        else
          if [[ $ref_DP_indel -eq $alt_DP_indel ]]; then
            geno_indel=1
          else
            geno_indel=2
          fi
        fi    
      fi
      Rscript --vanilla $scriptsWD/mpileupToCsv.R $id $file $allele_rs $geno_indel imputed_genoMissing.geno
    fi
  else 
    if [[ $missing_pos -eq 1 ]] && [[ $indel_pos -eq 0 ]]; then
      echo -e "\n--------------------------> Only INDEL position is missing, run Pipeline with INDEL as missing position: NA! \n"
      geno_indel="NA"
      Rscript --vanilla $scriptsWD/mpileupToCsv.R $id $file $allele_rs $geno_indel 0
    else 
        if [[ $missing_pos -eq 1 ]] && [[ $indel_pos -ne 0 ]]; then
          if [[ $impute == "0" ]]; then 
            echo -e "\n--------------------------> Missing position detected, check QC report! \n"
            exit 1
          else 
            echo -e "\n--------------------------> Imputation file detected, run VarScan and Pipeline! \n"
            echo -e "--------------------------> Check missing positions on the QC report and extract corresponding GT from the imputed VCF/BCF \n"
            cat $QC_file | grep "DP=0" | awk -F'=' '{print $2,$3}' | awk -F' ' '{print $1"\t"$3}' > missing_pos_toImpute.list
            bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n' -R missing_pos_toImpute.list $impute  > imputed_genoMissing_noStrand.geno
            awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print $3,$4, a[$1$2]}' imputed_genoMissing_noStrand.geno $hiris_strand_imputed | awk '{print $3,$4,$5,$6,$7,$1,$2}' | awk '{if ($4 == $7) {print $5, $6} else if ($3 == $7) {if ($5 == "0|0") {print "1|1", $6} else if ($5 == "1|1") {print "0|0", $6} else {print $5, $6}} else {print $5, $6}}' | paste imputed_genoMissing_noStrand.geno - | awk '{print $1, $2, $3, $4, $6, $7}' | sed -e "s/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g" > imputed_genoMissing.geno
            echo -e "--------------------------> INDEL position is covered, run VarScan and Pipeline! \n"
            samtools mpileup -f $ref -l $indel -B $bam | java -Xmx20g -jar $varscan mpileup2indel --min-coverage 2 --min-reads2 2 --output-vcf > "$id"_VarScan.vcf
            ref_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $5}' )
            alt_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $6}' )
            if [[ $ref_DP_indel -eq '' ]]; then 
              geno_indel=0
            else
              if [[ $ref_DP_indel -eq $alt_DP_indel ]]; then
                geno_indel=1
              else
                geno_indel=2
              fi
            fi            
            Rscript --vanilla $scriptsWD/mpileupToCsv.R $id $file $allele_rs $geno_indel imputed_genoMissing.geno
          fi
        fi
    fi
  fi
  if [[ $missing_pos -eq 0 ]]; then 
    echo -e "\n--------------------------> No missing positions detected, run VarScan and Pipeline! \n"
    samtools mpileup -f $ref -l $indel -B $bam | java -Xmx20g -jar $varscan mpileup2indel --min-coverage 2 --min-reads2 2 --output-vcf > "$id"_VarScan.vcf
    ref_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $5}' )
    alt_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $6}' )
    if [[ $ref_DP_indel -eq '' ]]; then 
      geno_indel=0
    else
      if [[ $ref_DP_indel -eq $alt_DP_indel ]]; then
        geno_indel=1
      else
        geno_indel=2
      fi
    fi
    Rscript --vanilla $scriptsWD/mpileupToCsv.R $id $file $allele_rs $geno_indel 0
  fi
fi

######################################################

#################### PROG FUNC 2 #####################

if [[ $prog_func == "2" ]]; then
  echo -e "--------------------------> Running program function $prog_func: \n"
  echo -e "--------------------------> Extract genotypes from VCF direct call for sample: $id \n"

  # Create ouput ID folder
  mkdir -p $output_folder/$id
  cd $output_folder/$id

  echo -e "--------------------------> Extract 41 Hirisplex positions using BCFtools and check strand: \n"

  bcftools mpileup -Ou --gvcf 0 -R $hiris_strand_imputed -f $ref $bam | bcftools call --gvcf 0 -m -Ov - > DirectCall.vcf
  /opt/software/ngs/angsd-0.932/htslib/bgzip DirectCall.vcf
  /opt/software/ngs/angsd-0.932/htslib/tabix -p vcf DirectCall.vcf.gz
  bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n' -R $hiris_strand_imputed DirectCall.vcf.gz > DirectCall.geno
  awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print $3,$4, a[$1$2]}' DirectCall.geno $hiris_strand_imputed | awk '{print $3,$4,$5,$6,$7,$1,$2}' | awk '{if ($4 == $7) {print $5, $6} else if ($3 == $7) {if ($5 == "0/0") {print "1/1", $6} else if ($5 == "1/1") {print "0/0", $6} else {print $5, $6}} else {print $5, $6}}' | paste DirectCall.geno - | awk '{print $1, $2, $3, $4, $6, $7}' | sed -e "s/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g" > DirectCall_final.geno

  echo -e "--------------------------> run VarScan for INDEL position! \n"
  samtools mpileup -f $ref -l $indel -B $bam | java -Xmx20g -jar $varscan mpileup2indel --min-coverage 2 --min-reads2 2 --output-vcf > "$id"_VarScan.vcf
  ref_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $5}' )
  alt_DP_indel=$(cat "$id"_VarScan.vcf |  grep -wE "16[[:space:]]89985753" | awk '{print $10}' | awk -F':' '{print $6}' )
  if [[ $ref_DP_indel -eq '' ]]; then 
    geno_indel=0
  else
    if [[ $ref_DP_indel -eq $alt_DP_indel ]]; then
      geno_indel=1
    else
      geno_indel=2
    fi
  fi
  Rscript --vanilla $scriptsWD/vcfToCsv.R $id $allele_rs $geno_indel DirectCall_final.geno 
fi

######################################################

#################### PROG FUNC 3 #####################

if [[ $prog_func == "3" ]]; then
  echo -e "--------------------------> Running program function $prog_func: \n"
  echo -e "--------------------------> Prediction of the hirisplex phenotypic traits and plot the obtained results for GL formatted file - for sample: $id \n"

  cd $output_folder/$id

  Rscript --vanilla $scriptsWD/phenotypicPrediction.R $id $hirisplex_out
  #Rscript --vanilla $scriptsWD/plotPhenoGL.R $id

fi

######################################################

#################### PROG FUNC 4 #####################

if [[ $prog_func == "4" ]]; then
  echo -e "--------------------------> Running program function $prog_func: \n"
  echo -e "--------------------------> Prediction of the hirisplex phenotypic traits and plot the obtained results for Direct Call formatted file for sample: $id \n"

  cd $output_folder/$id

  Rscript --vanilla $scriptsWD/phenotypicPrediction.R $id $hirisplex_out
  #Rscript --vanilla $scriptsWD/plotPhenoDC.R $id

fi

######################################################

