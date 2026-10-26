#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression
#SBATCH --job-name=munge-GWAS
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/01.04.Munge_sumstats_meta-analyses_by_GenomicCoord_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/01.04.Munge_sumstats_meta-analyses_by_GenomicCoord_%A_%a.out
#SBATCH --mem=20000
#SBATCH --time=04:00:00
#SBATCH --array=22


# Libraries

# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate ldsc


# Slurm array variables
wd=$(pwd)
array_file=${wd}/Example_arrays/array_meta-analyses_Diet_munge_by_GenomicCoords_ok.tsv
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
filename=$(basename ${file} .tsv.gz)
sample_size_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
a1_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f3)
a2_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f4)
sumstat_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f5)
null_value=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f6)
pval_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f7)
col_ID_ignore=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f8)


# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/Diet_project/ldsc_diet/ldsc_input/gwas_meta-analyses_parsed_and_with_ldsc_rsID
out_dir=${main}/Diet_project/ldsc_diet/ldsc_input/GWAS_meta-analyses_munge_sumstats_ldsc_input
path_packages=/home/csic/gpo/olg/ccs/packages

mkdir -p ${out_dir}
cd ${out_dir} # change to outdir


# Regex to detect numbers
re='^[0-9]+$'


# Initializate log
log_file=log_munge_sumstats_array.tsv
if [ ! -f ${log_file} ]
then
    echo "file SNPs_nonmissing pval_out_bounds pval_very_small" | tr " " "\t" > ${log_file}
fi


# IMPORTANT TO KNOW...
# Munge uses the --signed-sumstats argument to know which is a signed effect column (directionality included), and it also ask for the null value (0 or 1)
# It refers to the type of statistic being used
# From github code:
# null_values = {
#
#    'LOG_ODDS': 0,
#    'BETA': 0,
#    'OR': 1,
#    'Z': 0
#}


# BUG ALERT!!!!
# WARNING!! 
# There is a bug in munge_sumstats.py. Extremely low p-values lead to the wrong type (non-numeric) in columns['P']. Check out manual the fix provided in IMPORTANT_munge_sumstats_bug.txt file 

# Check if there are columns that must be ignored
if [ -z "${col_ID_ignore}" ]
then
    # Check if sample size is actually a number instead of a col
    if [[ ${sample_size_col} =~ $re ]]
    then
        # Run munge function from ldsc package, but using a default N size for all SNPs
        ${path_packages}/ldsc/munge_sumstats.py \
        --sumstats ${input_dir}/${file} \
        --N ${sample_size_col} \
        --snp ldsc_rsID \
        --a1 ${a1_col} \
        --a2 ${a2_col} \
        --signed-sumstats ${sumstat_col},${null_value} \
        --p ${pval_col} \
        --chunksize 500000 \
        --out ${filename} \
        --merge-alleles ${path_packages}/eur_w_ld_chr/w_hm3.snplist
    else
        # Run munge function from ldsc package, but with col sample size avialable. This is not usually the same for all the SNPs
        ${path_packages}/ldsc/munge_sumstats.py \
        --sumstats ${input_dir}/${file} \
        --N-col ${sample_size_col} \
        --snp ldsc_rsID \
        --a1 ${a1_col} \
        --a2 ${a2_col} \
        --signed-sumstats ${sumstat_col},${null_value} \
        --p ${pval_col} \
        --chunksize 500000 \
        --out ${filename} \
        --merge-alleles ${path_packages}/eur_w_ld_chr/w_hm3.snplist
    fi
else
    # Check if sample size is actually a number instead of a col
    if [[ ${sample_size_col} =~ $re ]]
    then
        # Run munge function from ldsc package, but using a defalt N size for all SNPs
        ${path_packages}/ldsc/munge_sumstats.py \
        --sumstats ${input_dir}/${file} \
        --N ${sample_size_col} \
        --ignore ${col_ID_ignore} \
        --snp ldsc_rsID \
        --a1 ${a1_col} \
        --a2 ${a2_col} \
        --signed-sumstats ${sumstat_col},${null_value} \
        --p ${pval_col} \
        --chunksize 500000 \
        --out ${filename} \
        --merge-alleles ${path_packages}/eur_w_ld_chr/w_hm3.snplist
    else
        # Run munge function from ldsc package, but with col sample size avialable. This is not usually the same for all the SNPs
        ${path_packages}/ldsc/munge_sumstats.py \
        --sumstats ${input_dir}/${file} \
        --N-col ${sample_size_col} \
        --ignore ${col_ID_ignore} \
        --snp ldsc_rsID \
        --a1 ${a1_col} \
        --a2 ${a2_col} \
        --signed-sumstats ${sumstat_col},${null_value} \
        --p ${pval_col} \
        --chunksize 500000 \
        --out ${filename} \
        --merge-alleles ${path_packages}/eur_w_ld_chr/w_hm3.snplist
    fi
fi

# Add important info to a logfile summarizing all the arrays
n_out=$(cat ${filename}.log | grep "bounds" | cut -d " " -f2)
nonmissing=$(cat ${filename}.log | grep "nonmissing" | cut -d "(" -f2 | cut -d " " -f1)
n_small=$(cat ${filename}.log | grep "floating" | cut -d ":" -f2 | cut -d " " -f2)
if [[ -z ${n_small} ]]
then
        n_small=0
fi

# Write
echo "${filename} ${nonmissing} ${n_out} ${n_small}" | tr " " "\t" >> ${log_file}