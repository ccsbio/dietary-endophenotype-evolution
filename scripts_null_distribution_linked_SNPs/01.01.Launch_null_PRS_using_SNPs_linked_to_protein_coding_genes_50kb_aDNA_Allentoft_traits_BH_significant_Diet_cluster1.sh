#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=null-PRS
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs
#SBATCH --cpus-per-task=20
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs/err/01.01.Launch_null_PRS_using_SNPs_linked_to_protein_coding_genes_50kb_aDNA_Allentoft_traits_BH_significant_Diet_cluster1_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs/out/01.01.Launch_null_PRS_using_SNPs_linked_to_protein_coding_genes_50kb_aDNA_Allentoft_traits_BH_significant_Diet_cluster1_%A_%a.out
#SBATCH --mem=60000
#SBATCH --time=05:00:00
#SBATCH --array=1-3

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Main paths
in_path=${LUSTRE}/ccs
path_scripts=$(pwd)
python_dir=${path_scripts}/python_scripts
SNP_sets_with_AF_path=${path_scripts}/files_SNPs_sets_with_AF
vcf_path=${in_path}/SNPs_linked_to_genes_Allentoft/Allentoft_data_final_set_protein_coding/sorted

# Array variables
cluster=cluster1
array_file=${path_scripts}/arrays_clusters/array_Diet_${cluster}_unique_SNPs_trend_analysis_BH_significant_traits.tsv
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
trait_clean="$(echo ${trait} | sed -r 's/_UKB.+//' | sed -r 's/_([0-9]+){2,}_.+//' | sed -r 's/_Nat.+//' | sed -r 's/_GIANT.+//' | sed -r 's/_GLGC.+//' | sed 's/_BioPsych//' | sed -r 's/indpMarkers.+//')"
n_prot="$(echo "${trait}" | tr "_" "\n" | grep "indpMarkers" | sed -r 's/indpMarkers//')"
trait_label_merged_GWAS="$(echo "${trait}" | sed -r 's/_[0-9]+indpMarkers//')"
echo ${trait_label_merged_GWAS}

# Variables
vcf_file=${vcf_path}/SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_all_chr_EUR_samples_sorted.vcf.gz
tabix=${vcf_path}/SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_all_chr_EUR_samples_sorted.vcf.gz.tbi
AF_prototype_SNPs=${SNP_sets_with_AF_path}/${cluster}/SNPs_${trait_clean}_${n_prot}_indpMarkers_with_EUR_AF_risk_alleles.tsv
null_set_SNPs=${SNP_sets_with_AF_path}/SNPs_linked_to_protein_coding_genes_Allentoft_AF_EUR.tsv
span=0.05
rep=1000
cpus=17
out_dir=${in_path}/Diet_project/Null_PRS_Allentoft_protein_coding/${cluster}/
info_outfile=${trait_clean}
seed=1234

# Filter SNPs associated using GWAS info
gwas=${in_path}/Diet_project/GWAS_matching_aDNA_data_1e-03pval/Merged_GWAS_meta_analyses_Diet_Cluster_1_logTG_prot_outer_pval_1e-03_common_and_fixed_Allentoft_aDNA.tsv.gz
pval_thres=1e-05

mkdir -p ${out_dir}


# Launch vegeta
python ${python_dir}/Null_distribution_PRS_SNPs_linked_to_genes_matching_AF.py \
--null_set_file ${null_set_SNPs} \
--prototype_SNPs_file ${AF_prototype_SNPs} \
--vcf_file ${vcf_file} \
--tabix_file ${tabix} \
--n_replicates ${rep} \
--n_cpus ${cpus} \
--seed ${seed} \
--AF_bin_span ${span} \
--out_dir ${out_dir} \
--info_outfile ${info_outfile} \
--remove_snps_from_gwas \
--gwas ${gwas} \
--pval_threshold ${pval_thres} \
--trait_label_merged_GWAS ${trait_label_merged_GWAS}
