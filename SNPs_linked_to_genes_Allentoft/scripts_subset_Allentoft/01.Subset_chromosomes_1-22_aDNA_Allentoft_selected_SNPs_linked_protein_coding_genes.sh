#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft
#SBATCH --job-name=vcf-filt
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/err/01.Subset_chromosomes_1-22_aDNA_Allentoft_selected_SNPs_linked_protein_coding_genes_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/out/01.Subset_chromosomes_1-22_aDNA_Allentoft_selected_SNPs_linked_protein_coding_genes_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules


# CESGA 
module load bedtools/2.30.0


# Paths
in_path=${LUSTRE}/ccs
wd=$(pwd)
vcf_path=${in_path}/DATA/ancient_DNA_Denmark/VCF_EURsamples_after_QC_diploid_imputed
bed_path=${wd}/SNP_files
out_path=${in_path}/SNPs_linked_to_genes_Allentoft/Allentoft_data_final_set_protein_coding
bedfile=SNPs_genome_wide_max50KB_linked_to_protein_coding_genes_GENCODEv47_hg19_levels2-3.bed

mkdir -p ${out_path}


# Subset by chromosome 
# If you want the vcf header add -header before -wa, which just prints the results
touch ${out_path}/vcf_filtering_name_list.txt
for i in {1..22}
do
    bedtools intersect -a ${vcf_path}/${i}.Allentoft_imputed_ph_EUR_967samples_diploid_INFO_0.8_MAF_0.01_GP_0.8_missing.vcf.gz -b ${bed_path}/${bedfile} -header -wa | bgzip > ${out_path}/SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_chr${i}.vcf.gz
    echo "SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_chr${i}.vcf.gz" >> ${out_path}/vcf_filtering_name_list.txt
done