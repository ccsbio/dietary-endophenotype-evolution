#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G
#SBATCH --job-name=vcf-sort-index
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/err/02.03.Sort_and_index_VCF_1000G_phase3_hg19_ALL_chr_Diet_SNPs_cluster3_with_GT_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/out/02.03.Sort_and_index_VCF_1000G_phase3_hg19_ALL_chr_Diet_SNPs_cluster3_with_GT_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=05:00:00
#SBATCH --array=1-22

# Load modules


# CESGA
module load bcftools/1.19


# Array
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
wd=$(pwd)
cluster=cluster3
array_file=${wd}/chromosomes.txt
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# Paths
vcf_path=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}
vcf_original=SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_chr${chr}.vcf.gz
vcf_name=$(basename ${vcf_original} .vcf.gz)
out_path=${vcf_path}/sorted

mkdir -p ${out_path}

# 1. Sort
bcftools sort ${vcf_path}/${vcf_original} -Oz -o ${out_path}/${vcf_name}_sorted.vcf.gz

# 2. Index
tabix -p vcf ${out_path}/${vcf_name}_sorted.vcf.gz