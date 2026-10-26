#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression
#SBATCH --job-name=ldsc-Diet
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/02.00.Ldsc_regression_Diet_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/02.00.Ldsc_regression_Diet_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=03:00:00
#SBATCH --array=1-94


# Libraries

# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate ldsc


# Variables
wd=$(pwd)
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/Diet_project/ldsc_diet/ldsc_input
path_packages=/home/csic/gpo/olg/ccs/packages
out_dir=${main}/Diet_project/ldsc_diet/ldsc_regression_Diet
log_file=ldsc_analysis_logs_Diet_GWAS.txt

mkdir -p ${out_dir}
cd ${out_dir} # change to outdir


# Log file that will store errors
if ! [[ -f ${out_dir}/${log_file} ]]
then
  echo "Creating an error/log file summarizing issues during the pairwise analyses performed with the current array of GWAS"
  touch ${out_dir}/${log_file}
fi


# Slurm array variables
array_file=${wd}/Example_arrays/array_Diet_ldsc_regression_ok.tsv
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
subdirectory=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
# The UKB files have different extension to remove. The idea would be to use only the code of the phenotype
if [[ ${subdirectory} == "UKB_ldsc_sumstats_input" ]]
then
  filename=${file%.gwas.imputed_v3.both_sexes*_rsID_merged.sumstats.gz}
else
  filename=$(basename ${file} _rsID_merged.sumstats.gz)
fi


# Find which files go after the current one, so the correlations has not yet previously computed between the following pairs
# Print lines after file1 and get filenames and subdirectories
files_to_iterate=( $(awk -F'\t' -v f="$file" '$1==f{seen=1;next} seen' "$array_file" | cut -f1) )
#( $(grep -A 600 ${file} ${array_file} | grep -v -w ${file} | cut -f1) )
subdirectories_to_iterate=( $(awk -F'\t' -v f="$file" '$1==f{seen=1;next} seen' "$array_file" | cut -f2) )
# Declare arrays to store grep output in a iterable way
declare -p files_to_itarate
declare -p subdirectories_to_iterate

# Get the length of the arrays
num_elements=${#files_to_iterate[@]}
echo ${num_elements}

# Iterate over all files except the one we're using
#for f s in ${files_to_itarate}
for ((i = 0; i < num_elements; i++))
do  
  f="${files_to_iterate[i]}"
  s="${subdirectories_to_iterate[i]}"
  
  echo ${f}
  echo ${s}

  # Get file name 2
  if [[ ${s} == "UKB_ldsc_sumstats_input" ]]
  then
    filename2=${f%.gwas.imputed_v3.both_sexes*_rsID_merged.sumstats.gz}
  else
    filename2=$(basename ${f} _rsID_merged.sumstats.gz)
  fi
  
  # Run munge function from ldsc package
  ${path_packages}/ldsc/ldsc.py \
  --out ${filename}.${filename2} \
  --rg ${input_dir}/${subdirectory}/${file},${input_dir}/${s}/${f} \
  --w-ld-chr ${path_packages}/eur_w_ld_chr/ \
  --ref-ld-chr ${path_packages}/eur_w_ld_chr/


  # After computing ldsc regression, check if the heritability estimate of the trait is negative. Because, otherwise does not make sense to continue. 
  # There won't be any genetic correlation to compute. 
  # Be sure that the pairwise comparison had a minimum number of SNPs, otherwise this is the only analysis where this happens  
  # 1. But first check if any error occured 
  ERROR_during_estimates=$(cat ${out_dir}/${filename}.${filename2}.log | grep "ERROR")

  if [[ -z ${ERROR_during_estimates} ]]
  then
    # 2. If there are no errors, check out if the heritability estimate is negative for the trait 1 (the file of the array that is iterating)
    h2_trait1=$(cat ${out_dir}/${filename}.${filename2}.log | grep "Total Observed scale h2" | cut -d " " -f5 | head -n1)
    h2_trait1_sign=$(cat ${out_dir}/${filename}.${filename2}.log | grep "Total Observed scale h2" | cut -d " " -f5 | head -n1 | grep "-")
    # Get Number of SNPs used
    N_SNPs=$(cat ${out_dir}/${filename}.${filename2}.log | grep "valid" | cut -d " " -f1)
    # If N of SNPs is >600K but heritability is negative, then BREAK the pairwise analysis for the current GWAS
    if [[ $h2_trait1_sign == $h2_trait1 ]] && [[ ${N_SNPs} -gt 599999 ]]
    then
      echo "Heritability estimate with >600K SNPs is negative for trait 1 from ${file}. Ldsc analysis has been interrupted because there are not genetic correlations to compute" >> ${log_file}
      break
    # If N of SNPs does not reach the minimum for enough statistical power and the heritibality estimate is negative, then print a warning because it can happen maybe just in this analysis
    elif [[ $h2_trait1_sign == $h2_trait1 ]] && [[ ${N_SNPs} -lt 599999 ]]
    then
      echo "Warning! Heritability estimate is negative for trait 1 when using "${N_SNPs}". It can be caused by the lower number of SNPs found between the pair "${filename}" and "${filename2} >> ${log_file}
    fi
  # Report if an error is reported by ldsc in the .log file for the current pair of GWAS
  else
    echo "An error occured during ldsc regression analysis between "${filename}" and "${filename2} >> ${log_file}
  fi
done

echo "Check out file ${log_file} with the summary of issues"
echo "Done"

