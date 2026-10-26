#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment
#SBATCH --job-name=data-merge
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/err/01.01.Launch_merging_alleles_GWAS_and_1000G_Diet_SNPs_cluster1_1e-03_fix_VCF_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/out/01.01.Launch_merging_alleles_GWAS_and_1000G_Diet_SNPs_cluster1_1e-03_fix_VCF_%j.out
#SBATCH --mem=30000
#SBATCH --time=04:00:00

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last

# Main paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
path_scripts=$(pwd)
#path_scripts=${wd}/scripts_SNP_ascertainment

# Variables
pval=1e-03
cluster=cluster1
cluster_name_gwas=Cluster_1_logTG_prot
gwas_path=${in_path}/Diet_project/GWAS_matching_aDNA_data_${pval}pval/
gwas_file=${gwas_path}Merged_GWAS_meta_analyses_Diet_${cluster_name_gwas}_outer_pval_${pval}_common_and_fixed_Allentoft_aDNA.tsv.gz
vcf_path=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted/European_populations/
vcf=SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.EUR.recalc.vcf.gz
python_dir=${path_scripts}/python_scripts
filename_out=Diet_${cluster}_SNPs_${pval}_1000G_ph3_merge_GWAS_and_aDNA.EUR.recalc


# Launch vegeta
python ${python_dir}/Merging_alleles_GWAS_1000G_data_fix_VCF.py --gwas_data_file ${gwas_file} --in_dir ${vcf_path} --vcf ${vcf} --out_dir ${vcf_path} --info_outfile ${filename_out} 


# Add VCF header
#zcat ${vcf_path}/${vcf} | grep "#" > ${vcf_path}/header.txt
cat <(zcat ${vcf_path}/${vcf} | grep "#") <(zcat ${vcf_path}/${filename_out}.vcf.gz | tail -n+2) | bgzip > ${vcf_path}/${filename_out}.header.vcf.gz
rm ${vcf_path}/${filename_out}.vcf.gz


# Load modules for latest steps
module --force purge
module load cesga/2020 bcftools/1.19
# 1. Sort
bcftools sort ${vcf_path}/${filename_out}.header.vcf.gz -Oz -o ${vcf_path}/${filename_out}.sorted.vcf.gz
rm ${vcf_path}/${filename_out}.header.vcf.gz

# 2. Index
tabix -p vcf ${vcf_path}/${filename_out}.sorted.vcf.gz