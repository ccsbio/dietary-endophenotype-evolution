#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes
#SBATCH --job-name=snps-closest
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/err/01.01.Subset_SNPs_from_Allentoft_aDNA_linked_to_genes_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/out/01.01.Subset_SNPs_from_Allentoft_aDNA_linked_to_genes_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=03:00:00
#SBATCH --array=1-22

# Load modules


# CESGA 
module load bedtools/2.30.0


# Array variables
in_path=${LUSTRE}/ccs
wd=$(pwd)
array_file=${wd}/array_chr.txt
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)


# Paths
bedfile_gene_annotation=${in_path}/DATA/GENCODE/hg19_release_47/gencode.v47lift37.annotation_protein_coding_genes_level2-3.no_complex_regions.bed
vcf_path=${in_path}/DATA/ancient_DNA_Denmark/VCF_EURsamples_after_QC_diploid_imputed
bed_path=${vcf_path}/bedfiles
out_path=${wd}/bedfiles_results_snps_Allentoft_closest_to_genes

mkdir -p ${out_path}
mkdir -p ${bed_path}


# Subset by chromosome 
# bedtools closest requires that all input files are presorted data by chromosome and then by start position (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files).
# -D option will use negative integers to report distances to “upstream” features.
# 1. First create a bedfile
zcat ${vcf_path}/${chr}.Allentoft_imputed_ph_EUR_967samples_diploid_INFO_0.8_MAF_0.01_GP_0.8_missing.vcf.gz | grep -v "#" | cut -f1-3 | awk 'BEGIN{FS=OFS="\t"}{print "chr"$1, $2-1, $2, $3}' | sort -k1,1n -k2,2n > ${bed_path}/chr${chr}_Allentoft_imputed_ph_afterQC.bed

genes=$(awk -v chr="chr"${chr} 'BEGIN{FS=OFS="\t"}{if ($1==chr) print}' ${bedfile_gene_annotation} | wc -l)
echo ${genes}

# 2. Now check out SNPs linked to genes
# -D ref
bedtools closest -a <(awk -v chr="chr"${chr} 'BEGIN{FS=OFS="\t"}{if ($1==chr) print}' ${bedfile_gene_annotation} | sort -k1,1 -k2,2n) -b ${bed_path}/chr${chr}_Allentoft_imputed_ph_afterQC.bed -d > ${out_path}/${chr}.SNPs_Allentoft_afterQC_linked_to_genes_GENCODEv47.bed
