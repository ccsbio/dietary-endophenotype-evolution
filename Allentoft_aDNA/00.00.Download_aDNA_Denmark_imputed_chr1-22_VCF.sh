#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=allentoft
#SBATCH -D /home/csic/gpo/olg/ccs/Allentoft_aDNA
#SBATCH --ntasks=6 
#SBATCH --error=/home/csic/gpo/olg/ccs/Allentoft_aDNA/err/00.00.Download_aDNA_Denmark_imputed_chr1-22_VCF_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Allentoft_aDNA/out/00.00.Download_aDNA_Denmark_imputed_chr1-22_VCF_%j.out
#SBATCH --mem=10000
#SBATCH --time=05:00:00

# Variables
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg
outpath=${main}/ccs/DATA/ancient_DNA_Denmark/VCF

mkdir -p ${outpath}

# Download VCF for all autosomal chromosomes and tabix
for i in {1..22}
do
    wget https://erda.ku.dk/archives/917f1ac64148c3800ab7baa29402d088/neo/genotypes-public/${i}.neo_impute_ph.vcf.gz -P ${outpath}
    wget https://erda.ku.dk/archives/917f1ac64148c3800ab7baa29402d088/neo/genotypes-public/${i}.neo_impute_ph.vcf.gz.tbi -P ${outpath}
done

wget https://erda.ku.dk/archives/917f1ac64148c3800ab7baa29402d088/neo/genotypes-public/README.md -P ${outpath}