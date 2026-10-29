#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Download_archaic_genomes
#SBATCH --job-name=vindija-download
#SBATCH --error=/home/csic/gpo/olg/ccs/Download_archaic_genomes/err/00.01.Download_Vindija_chr1-22_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Download_archaic_genomes/out/00.01.Download_Vindija_chr1-22_%j.out
#SBATCH --mem=10000
#SBATCH --time=05:00:00

# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
outpath=${main}/DATA/Archaic_genomes/Neanderthal/Vindija/

mkdir -p ${outpath}

# Download VCF for all autosomal chromosomes
wget http://ftp.eva.mpg.de/neandertal/Vindija/VCF/README -P ${outpath}
for i in {1..22}
do
    wget http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${i}_mq25_mapab100.vcf.gz -P ${outpath}
done