#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft
#SBATCH --job-name=vcf-filt
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/err/01.03.Subset_chromosomes_1-22_aDNA_Allentoft_Diet_array_clusters_outer_merge_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/out/01.03.Subset_chromosomes_1-22_aDNA_Allentoft_Diet_array_clusters_outer_merge_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00
#SBATCH --array=1-4

# Load modules


# CESGA 
module load bedtools/2.30.0


# WARNING!
# This script only can wok with some conditions:
# - SNP file name is expected to have words connected with '_'
# - SNP file name contains the pval threshold without spaces near to 'pval' (ex: 1e-08pval)
# - 'cluster' variable in the array file and extracted 'pval threshold' from the filename are expected to match strings in the subdirectories with SNP files for that cluster


# Array and global variables
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg/ccs
wd=$(pwd)
array_file=${wd}/array_Allentoft_SNPs_4clusters_diet_pval1e-03_outer_merge.tsv
cluster=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
snps_bedfile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
pval_threshold=$(echo ${snps_bedfile} | tr "_" "\n" | grep "pval" | sed 's/pval\.bed//')


# Paths
vcf_path=${main}/DATA/ancient_DNA_Denmark/VCF_EURsamples_after_QC_diploid_imputed
bed_path=${wd}/SNP_files/outer_merge_${pval_threshold}
bedfile=${bed_path}/${snps_bedfile}
out_path=${main}/Diet_project/aDNA_trend_analysis/Allentoft_VCF_afterQC/outer_merge_pval_${pval_threshold}/${cluster}

mkdir -p ${out_path}


# Subset by chromosome 
# If you want the vcf header add -header before -wa, which just prints the results
touch ${out_path}/vcf_filtering_name_list.txt
for i in {1..22}
do
   bedtools intersect -a ${vcf_path}/${i}.Allentoft_imputed_ph_EUR_967samples_diploid_INFO_0.8_MAF_0.01_GP_0.8_missing.vcf.gz -b ${bedfile} -header -wa | \
   bgzip > ${out_path}/Diet_SNPs_aDNA_Allentoft_afterQC_chr${i}_${pval_threshold}_outer_merge.vcf.gz
   echo "Diet_SNPs_aDNA_Allentoft_afterQC_chr${i}_${pval_threshold}_outer_merge.vcf.gz" >> ${out_path}/vcf_filtering_name_list.txt
done
