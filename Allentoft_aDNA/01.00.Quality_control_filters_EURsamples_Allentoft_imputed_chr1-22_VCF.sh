#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=QC-diploid-Allentoft
#SBATCH -D /home/csic/gpo/olg/ccs/Allentoft_aDNA
#SBATCH --error=/home/csic/gpo/olg/ccs/Allentoft_aDNA/err/01.00.Quality_control_filters_EURsamples_Allentoft_imputed_chr1-22_VCF_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Allentoft_aDNA/out/01.00.Quality_control_filters_EURsamples_Allentoft_imputed_chr1-22_VCF_%A_%a.out
#SBATCH --mem=20000
#SBATCH --time=05:00:00
#SBATCH --array=1-22

# Load modules
# CESGA
module load bcftools/1.19




# Variables
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg/ccs
inpath=${main}/DATA/ancient_DNA_Denmark/VCF
outpath=${main}/DATA/ancient_DNA_Denmark/VCF_EURsamples_after_QC_diploid_imputed
sample_file=./sample_files/Sample_file_967_aDNA_EUR_ok.txt
array_file=./array_chromosomes_Allentoft.txt
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

mkdir -p ${outpath}

# Quality filters:
# - Just European samples in the Allentoft dataset
# - Remove variants with imputation INFO score < 0.8
# - Remove variants with MAF < 0.01
# - Just diploid data remains. Pseudohaploid calls were removed
# - Genotypes with genotype posterior probabilities < 0.8 are converted to missing data
# - SNPs from 1000G MASK and 1240K array remain 

bcftools view -S ${sample_file} -i 'MAF >= 0.01 & INFO >= 0.8' ${inpath}/${chr}.neo_impute_ph.vcf.gz | bcftools +setGT -Oz -- -t q -n . -e 'GP[*] >= 0.8' > ${outpath}/${chr}.Allentoft_imputed_ph_EUR_967samples_diploid_INFO_0.8_MAF_0.01_GP_0.8_missing.vcf.gz


#for i in {1..22}
#do
#done