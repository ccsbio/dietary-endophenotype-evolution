#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics
#SBATCH --job-name=vcf-concat
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics/err/03.02.Concat_VCF_Vindija_ALL_chr_Diet_SNPs_cluster2_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics/out/03.02.Concat_VCF_Vindija_ALL_chr_Diet_SNPs_cluster2_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules
# CESGA 
module load bcftools/1.19


# Paths
in_path=${LUSTRE}/ccs
cluster=cluster2
vcf_path=${in_path}/Diet_project/Vindija_VCFs_Diet_clusters/${cluster}/sorted
file_with_vcf_filenames_list_to_concat=${vcf_path}/vcf_filtering_name_list_sorted.txt # it is expected that each line depicts a vcf filename


# vcf filename list stored in one line 
cd ${vcf_path}
vcf_listed_one_line=$(printf "%s\n" *.vcf.gz | sort -V | paste -sd' ' -)
echo ${vcf_listed_one_line}


# Bcftools command to merge properly vcf (headers are also handled)
# IMPORTANT!! VCF files must be sorted, indexed and given with proper chromosome vcf order if you want the final output sorted
bcftools concat -a -Oz -o ${vcf_path}/SNPs_Diet_${cluster}_aDNA_Vindija_all_chr_sorted.vcf.gz ${vcf_listed_one_line}

# Index
tabix -p vcf ${vcf_path}/SNPs_Diet_${cluster}_aDNA_Vindija_all_chr_sorted.vcf.gz
