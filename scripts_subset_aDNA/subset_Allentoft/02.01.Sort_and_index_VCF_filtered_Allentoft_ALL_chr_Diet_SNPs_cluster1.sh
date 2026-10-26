#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft
#SBATCH --job-name=vcf-sort-index
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/err/02.01.Sort_and_index_VCF_filtered_Allentoft_ALL_chr_Diet_SNPs_cluster1_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/out/02.01.Sort_and_index_VCF_filtered_Allentoft_ALL_chr_Diet_SNPs_cluster1_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=05:00:00
#SBATCH --array=1-22

# Load modules


# CESGA
module load bcftools/1.19


# Array
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
wd=$(pwd)
cluster=cluster1
pval_intermediate=1e-03
array_file=${wd}/chromosomes.txt
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# Paths
vcf_path=${in_path}/Diet_project/aDNA_trend_analysis/Allentoft_VCF_afterQC/outer_merge_pval_${pval_intermediate}/${cluster}
vcf_original=Diet_SNPs_aDNA_Allentoft_afterQC_chr${chr}_${pval_intermediate}_outer_merge.vcf.gz
vcf_name=$(basename ${vcf_original} .vcf.gz)
out_path=${vcf_path}/sorted

mkdir -p ${out_path}

# 1. Sort
bcftools sort ${vcf_path}/${vcf_original} -Oz -o ${out_path}/${vcf_name}_sorted.vcf.gz

# 2. Index
tabix -p vcf ${out_path}/${vcf_name}_sorted.vcf.gz