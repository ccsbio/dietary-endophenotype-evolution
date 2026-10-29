#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft
#SBATCH --job-name=vcf-concat
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/err/03.Concat_Allentoft_ALL_chr_filtered_VCF_selected_SNPs_linked_protein_coding_genes_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/out/03.Concat_Allentoft_ALL_chr_filtered_VCF_selected_SNPs_linked_protein_coding_genes_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules


# CESGA
module load bcftools/1.19


# Paths
in_path=${LUSTRE}/ccs
vcf_path=${in_path}/SNPs_linked_to_genes_Allentoft/Allentoft_data_final_set_protein_coding/sorted
#file_with_vcf_filenames_list_to_concat=$(ls ${vcf_path}) # it is expected that each line depicts a vcf filename


# vcf filename list stored in one line 
#vcf_listed_one_line=$(echo "$(cat ${file_with_vcf_filenames_list_to_concat} | tr "\n" " ")")
#echo ${vcf_listed_one_line}


# Bcftools command to merge properly vcf (headers are also handled)
# IMPORTANT!! VCF files must be sorted, indexed and given with proper chromosome vcf order if you want the final output sorted
cd ${vcf_path}
vcf_listed_one_line=$(printf "%s\n" *.vcf.gz | sort -V | paste -sd' ' -)
echo ${vcf_listed_one_line}

bcftools concat -a -Oz -o ${vcf_path}/SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_all_chr_EUR_samples_sorted.vcf.gz ${vcf_listed_one_line}

# Index
tabix -p vcf ${vcf_path}/SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_all_chr_EUR_samples_sorted.vcf.gz
