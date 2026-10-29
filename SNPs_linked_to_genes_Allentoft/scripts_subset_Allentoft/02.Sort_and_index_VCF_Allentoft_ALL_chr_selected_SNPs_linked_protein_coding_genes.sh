#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft
#SBATCH --job-name=vcf-sort-index
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/err/02.Sort_and_index_VCF_Allentoft_ALL_chr_selected_SNPs_linked_protein_coding_genes_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/out/02.Sort_and_index_VCF_Allentoft_ALL_chr_selected_SNPs_linked_protein_coding_genes_%A_%a.out
#SBATCH --mem=20000
#SBATCH --time=06:00:00
#SBATCH --array=1-22

# Load modules


# CESGA
module load bcftools/1.19


# Array
in_path=${LUSTRE}/ccs
wd=$(pwd)
array_file=${wd}/chromosomes.txt
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# Paths
vcf_path=${in_path}/SNPs_linked_to_genes_Allentoft/Allentoft_data_final_set_protein_coding
vcf_original=SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_chr${chr}.vcf.gz
vcf_name=$(basename ${vcf_original} .vcf.gz)
out_path=${vcf_path}/sorted

mkdir -p ${out_path}

# 1. Sort
bcftools sort ${vcf_path}/${vcf_original} -Oz -o ${out_path}/${vcf_name}_sorted.vcf.gz

# 2. Index
tabix -p vcf ${out_path}/${vcf_name}_sorted.vcf.gz