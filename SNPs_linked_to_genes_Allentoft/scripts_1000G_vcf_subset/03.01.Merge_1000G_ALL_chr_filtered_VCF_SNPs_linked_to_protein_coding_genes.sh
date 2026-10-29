#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset
#SBATCH --job-name=vcf-merge
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/err/03.01.Merge_1000G_ALL_chr_filtered_VCF_SNPs_linked_to_protein_coding_genes_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/out/03.01.Merge_1000G_ALL_chr_filtered_VCF_SNPs_linked_to_protein_coding_genes_%j.out
#SBATCH --mem=5000
#SBATCH --time=04:00:00

# Load modules



# Paths
in_path=${LUSTRE}/ccs
vcf_path=${in_path}/SNPs_linked_to_genes_Allentoft/1000G_data_final_set_protein_coding/VCF_filt_9_cols
vcf_with_header=SNPs_linked_to_protein_coding_genes_chr1_9fields.vcf.gz

# Create command for all VCFs to merge
# Double parenthesis so it becomes an array
files_to_iterate=($( ls ${vcf_path}/*.gz | grep -v ${vcf_path}/${vcf_with_header} ))

# Get the length of the arrays
num_elements=${#files_to_iterate[@]}
echo ${num_elements}

# Create the file with the first VCF
echo " <(zcat ${vcf_path}/${vcf_with_header})" > ${vcf_path}/command.txt
# Now add the rest
for ((i = 0; i < num_elements; i++))
do  
  f="${files_to_iterate[i]}"
  echo " <(zcat ${f})" >> ${vcf_path}/command.txt
done  

commands=$(echo "cat $(cat ${vcf_path}/command.txt | tr "\n" " ")")
echo ${commands}
eval ${commands} | bgzip > ${vcf_path}/SNPs_linked_to_protein_coding_genes_and_Allentoft_all_chr_1000G_ph3_hg19_9fields.vcf.gz

# Remove intermediate file
rm ${vcf_path}/command.txt
