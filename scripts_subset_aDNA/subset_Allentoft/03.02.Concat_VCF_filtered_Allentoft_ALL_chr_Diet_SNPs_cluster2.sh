#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft
#SBATCH --job-name=vcf-concat
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/err/03.02.Concat_VCF_filtered_Allentoft_ALL_chr_Diet_SNPs_cluster2_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/out/03.02.Concat_VCF_filtered_Allentoft_ALL_chr_Diet_SNPs_cluster2_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules


# CESGA
module load bcftools/1.19

# Paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
wd=$(pwd)
cluster=cluster2
pval_intermediate=1e-03
vcf_path=${in_path}/Diet_project/aDNA_trend_analysis/Allentoft_VCF_afterQC/outer_merge_pval_${pval_intermediate}/${cluster}/sorted
file_with_vcf_filenames_list_to_concat=${vcf_path}/vcf_filtering_name_list_sorted.txt # it is expected that each line depicts a vcf filename
out_filename=Diet_SNPs_${cluster}_${pval_intermediate}_outer_merge_aDNA_Allentoft_afterQC_all_chr_EUR_samples

# vcf filename list stored in one line 
cd ${vcf_path}
ls *.vcf.gz | sort -V > ${file_with_vcf_filenames_list_to_concat}
vcf_listed_one_line=$(echo "$(cat ${file_with_vcf_filenames_list_to_concat} | tr "\n" " ")")
echo ${vcf_listed_one_line}


# Bcftools command to merge properly vcf (headers are also handled)
# IMPORTANT!! VCF files must be sorted, indexed and given with proper chromosome vcf order if you want the final output sorted
bcftools concat -a -Oz -o ${vcf_path}/${out_filename}.sorted.vcf.gz ${vcf_listed_one_line}

# Index
tabix -p vcf ${vcf_path}/${out_filename}.sorted.vcf.gz
