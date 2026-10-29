#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes
#SBATCH --job-name=filt-gencode
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/err/00.01.Filter_complex_LD_regions_GENCODE_gene_annotation_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/out/00.01.Filter_complex_LD_regions_GENCODE_gene_annotation_%j.out
#SBATCH --mem=10000
#SBATCH --time=03:00:00


# Load modules


# CESGA 
module load bedtools/2.30.0


# Paths
wd=$(pwd)
in_path=${LUSTRE}/ccs/DATA/GENCODE/hg19_release_47
bedname=gencode.v47lift37.annotation_protein_coding_genes_level2-3
bedfile_gene_annotation=${in_path}/${bedname}.bed
bedfile_complex_LD_regions=${wd}/complex_regions.bed


# Intersect bedfiles with genecode annotation and complex regions (high LD, hypervariable...)
# -v returns non-overlapping regions
bedtools intersect -v -a ${bedfile_gene_annotation} -b ${bedfile_complex_LD_regions} > ${in_path}/${bedname}.no_complex_regions.bed