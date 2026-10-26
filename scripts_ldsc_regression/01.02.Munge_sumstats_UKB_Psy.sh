#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/
#SBATCH --job-name=munge-UKB
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/01.02.Munge_sumstats_UKB_Psy_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/01.02.Munge_sumstats_UKB_Psy_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=02:00:00
#SBATCH --array=1-6


# Libraries

# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate ldsc


# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/Diet_project/ldsc_diet/ldsc_input/UKB_parsed_and_with_ldsc_rsID/UKB_Psychiatric
out_dir=${main}/Diet_project/ldsc_diet/ldsc_input/UKB_ldsc_sumstats_input
path_packages=/home/csic/gpo/olg/ccs/packages

mkdir -p ${out_dir}
cd ${out_dir} # change to outdir


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


# Slurm array variables
array_file=${input_dir}/array_UKB_Psychiatric_with_rsID.txt
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file})
filename=$(basename ${file} .tsv.gz)


# Run munge function from ldsc package
# --ignore variant \
${path_packages}/ldsc/munge_sumstats.py \
--sumstats ${input_dir}/${file} \
--N-col n_complete_samples \
--snp rsID \
--a1 alt \
--a2 ref \
--chunksize 500000 \
--out ${filename} \
--merge-alleles ${path_packages}/eur_w_ld_chr/w_hm3.snplist


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

