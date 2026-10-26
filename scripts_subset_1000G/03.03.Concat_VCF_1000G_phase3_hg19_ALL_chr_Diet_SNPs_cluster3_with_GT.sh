#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G
#SBATCH --job-name=vcf-concat
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/err/03.03.Concat_VCF_1000G_phase3_hg19_ALL_chr_Diet_SNPs_cluster3_with_GT_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/out/03.03.Concat_VCF_1000G_phase3_hg19_ALL_chr_Diet_SNPs_cluster3_with_GT_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules


# CESGA
module load bcftools/1.19

# Paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
wd=$(pwd)
cluster=cluster3
vcf_path=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted
file_with_vcf_filenames_list_to_concat=${vcf_path}/vcf_filtering_name_list_sorted.txt # it is expected that each line depicts a vcf filename


# vcf filename list stored in one line 
cd ${vcf_path}
ls *.vcf.gz | sort -V > ${file_with_vcf_filenames_list_to_concat}
vcf_listed_one_line=$(echo "$(cat ${file_with_vcf_filenames_list_to_concat} | tr "\n" " ")")
echo ${vcf_listed_one_line}


# Bcftools command to merge properly vcf (headers are also handled)
# IMPORTANT!! VCF files must be sorted, indexed and given with proper chromosome vcf order if you want the final output sorted
bcftools concat -a -Oz -o ${vcf_path}/SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.vcf.gz ${vcf_listed_one_line}

# Index
tabix -p vcf ${vcf_path}/SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.vcf.gz
