#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset
#SBATCH --job-name=vcf-filt
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/err/01.01.Subset_SNPs_linked_to_protein_coding_genes_chr1-22_1000G_hg19_%j.err
#SBATCH --output=//home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/out/01.01.Subset_SNPs_linked_to_protein_coding_genes_chr1-22_1000G_hg19_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules


# CESGA 
module load bedtools/2.30.0


# Paths
in_path=${LUSTRE}/ccs
wd=$(pwd)
vcf_path=${in_path}/DATA/1000G_human_variation_phase3/VCFs_1000G_phase3_hg19
out_path=${in_path}/SNPs_linked_to_genes_Allentoft/1000G_data_final_set_protein_coding
bedfile=${wd}/SNP_files/SNPs_genome_wide_max50KB_linked_to_protein_coding_genes_GENCODEv47_hg19_levels2-3.bed

mkdir -p ${out_path}

# Subset by chromosome 
# If you want the vcf header add -header before -wa, which just prints the results
for i in {1..22}
do
   if [ $i -eq 1 ]
   then
      bedtools intersect -a ${vcf_path}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -b ${bedfile} -header -wa | bgzip > ${out_path}/SNPs_linked_to_protein_coding_genes_chr${i}.vcf.gz
   else
      bedtools intersect -a ${vcf_path}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -b ${bedfile} -wa | bgzip > ${out_path}/SNPs_linked_to_protein_coding_genes_chr${i}.vcf.gz
   fi
done
