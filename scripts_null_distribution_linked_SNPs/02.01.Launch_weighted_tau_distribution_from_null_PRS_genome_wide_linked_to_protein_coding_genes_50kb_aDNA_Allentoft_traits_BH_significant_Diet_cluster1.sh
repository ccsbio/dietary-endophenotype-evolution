#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=tau-null-PRS
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs
#SBATCH --cpus-per-task=20
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs/err/02.01.Launch_weighted_tau_distribution_from_null_PRS_genome_wide_linked_to_protein_coding_genes_50kb_aDNA_Allentoft_traits_BH_significant_Diet_cluster1_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs/out/02.01.Launch_weighted_tau_distribution_from_null_PRS_genome_wide_linked_to_protein_coding_genes_50kb_aDNA_Allentoft_traits_BH_significant_Diet_cluster1_%j.out
#SBATCH --mem=20000
#SBATCH --time=05:00:00

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Main paths
in_path=${LUSTRE}/ccs
path_scripts=$(pwd)
python_dir=${path_scripts}/python_scripts


# Variables
cluster=cluster1
null_PRS_file=${in_path}/Diet_project/Null_PRS_Allentoft_protein_coding/${cluster}/processed_df/Null_PRS_1000replicates_allentoft_SNPs_50Kb_togenes_${cluster}.tsv
observed_tau_log=${in_path}/Diet_project/aDNA_trend_analysis/Allentoft_results_proportion_risk_alleles/${cluster}/permutation_test_Tau_filter_outlier_bins_10K/log_permutation_tau_res_${cluster}.txt
cpus=17
out_dir=${in_path}/Diet_project/Null_PRS_Allentoft_protein_coding/${cluster}/weighted_tau_from_null_PRS_BH_significant/
info_outfile=Diet_${cluster}
log_file=${out_dir}/log_weighted_tau_from_null_PRS_genome_wide_${cluster}.txt
metadata_file=${path_scripts}/metadata_samples_trend_analysis/Metadata_Allentoft_EUR_samples_400-15KyearsBP_946indv_diploid_imputed_noCanary.tsv
y_var=Age.average
b_size=500
b_start=400
b_end=15000

mkdir -p ${out_dir}


# Initializate log
if [ ! -f ${log_file} ]
then
    echo "trait tau_obs_res pval_greater pval_lower n_rep test_type filter_strategy" | tr " " "\t" > ${log_file}
fi


# Launch vegeta
python ${python_dir}/From_null_PRS_to_weighted_tau_distribution.py \
--null_PRS_cluster_results ${null_PRS_file} \
--n_cpus ${cpus} \
--out_dir ${out_dir} \
--info_outfile ${info_outfile} \
--metadata_samples ${metadata_file} \
--y_variable ${y_var} \
--log_filename ${log_file} \
--tau_log_observed ${observed_tau_log} \
--bin_start ${b_start} \
--bin_end ${b_end} \
--bin_size ${b_size} \
--filt_outlier_bins_by_y 10000
