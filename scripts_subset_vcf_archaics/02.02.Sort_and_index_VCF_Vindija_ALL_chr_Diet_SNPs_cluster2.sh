#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics
#SBATCH --job-name=vcf-sort-index
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics/err/02.02.Sort_and_index_VCF_Vindija_ALL_chr_Diet_SNPs_cluster2_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics/out/02.02.Sort_and_index_VCF_Vindija_ALL_chr_Diet_SNPs_cluster2_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00
#SBATCH --array=1-22

# Load modules


# CESGA 
module load bcftools/1.19


# Array
in_path=${LUSTRE}/ccs
wd=$(pwd)
cluster=cluster2
array_file=${wd}/chromosomes.txt
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# Paths
vcf_path=${in_path}/Diet_project/Vindija_VCFs_Diet_clusters/${cluster}
vcf_original=Vindija_SNPs_Diet_${cluster}_chr${chr}.vcf.gz
vcf_name=$(basename ${vcf_original} .vcf.gz)
out_path=${vcf_path}/sorted

mkdir -p ${out_path}


# 1. Sort
bcftools sort ${vcf_path}/${vcf_original} -Oz -o ${out_path}/${vcf_name}_sorted.vcf.gz

# 2. Index
tabix -p vcf ${out_path}/${vcf_name}_sorted.vcf.gz