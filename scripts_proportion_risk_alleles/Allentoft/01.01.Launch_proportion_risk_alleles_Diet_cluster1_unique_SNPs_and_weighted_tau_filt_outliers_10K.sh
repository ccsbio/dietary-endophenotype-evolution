#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=prop-risk-allele
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_proportion_risk_alleles/Allentoft
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_proportion_risk_alleles/Allentoft/err/01.01.Launch_proportion_risk_alleles_Diet_cluster1_unique_SNPs_and_weighted_tau_filt_outliers_10K_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_proportion_risk_alleles/Allentoft/out/01.01.Launch_proportion_risk_alleles_Diet_cluster1_unique_SNPs_and_weighted_tau_filt_outliers_10K_%A_%a.out
#SBATCH --mem=40000
#SBATCH --time=04:00:00
#SBATCH --array=3-9

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Main paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
path_scripts=$(pwd)
pval_intermediate_files=1e-03
cluster=cluster1

# Array variables
array_file=${path_scripts}/array_Diet_${cluster}_unique_SNPs_trend_analysis.tsv
SNP_dict_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
info_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)

# Variables
python_dir=${path_scripts}/python_scripts
SNP_dict_path=${HOME}/ccs/Diet_project/scripts_SNP_ascertainment/SNP_dictionaries/${cluster}
vcf_file=Diet_SNPs_${cluster}_${pval_intermediate_files}_outer_merge_aDNA_Allentoft_afterQC_all_chr_EUR_samples.sorted.vcf.gz
tabix=${vcf_file}.tbi
path_vcf=${in_path}/Diet_project/aDNA_trend_analysis/Allentoft_VCF_afterQC/outer_merge_pval_${pval_intermediate_files}/${cluster}/sorted/
out_dir=${in_path}/Diet_project/aDNA_trend_analysis/Allentoft_results_proportion_risk_alleles/${cluster}/
log_file=${out_dir}/log_observed_tau_res.txt
metadata_file=${path_scripts}/metadata_samples_trend_analysis/Metadata_Allentoft_EUR_samples_400-15KyearsBP_946indv_diploid_imputed_noCanary.tsv
y_var=Age.average
b_size=500
b_start=400
b_end=15000

mkdir -p ${out_dir}


# Initializate log
if [ ! -f ${log_file} ]
then
    echo "trait tau_obs_res pval_tau n_indep_markers bool_weights percentage_ties filter_strategy" | tr " " "\t" > ${log_file}
fi


# Launch vegeta
python ${python_dir}/calculate_proportion_risk_alleles_Allentoft_and_weighted_or_regular_tau_using_ancient_DNA_and_SNP_dict.py \
--in_dir ${path_vcf} \
--file ${vcf_file} \
--tabix ${tabix} \
--SNP_dict ${SNP_dict_path}/${SNP_dict_file} \
--out_dir ${out_dir} \
--info_outfile ${info_file} \
--metadata_samples ${metadata_file} \
--y_variable ${y_var} \
--log_filename ${log_file} \
--weight_y_var \
--bin_start ${b_start} \
--bin_end ${b_end} \
--bin_size ${b_size} \
--filt_outlier_bins_by_y 10000