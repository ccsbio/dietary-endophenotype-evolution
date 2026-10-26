#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment
#SBATCH --job-name=BED-1000G
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/err/02.04.Make_BED_from_VCF_Diet_cluster4_SNPs_1e-03_1000G_EUR_pop_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/out/02.04.Make_BED_from_VCF_Diet_cluster4_SNPs_1e-03_1000G_EUR_pop_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Cesga
module load plink/1.9b5



# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
pval=1e-03
cluster=cluster4
inpath=${main}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted/European_populations
vcf=Diet_${cluster}_SNPs_${pval}_1000G_ph3_merge_GWAS_and_aDNA.EUR.recalc.sorted.vcf.gz
bed_path=${inpath}/BED_Diet_${cluster}_SNPs_${pval}_1000G_EUR

mkdir -p ${bed_path}


# Create bedfile for ASD-risk SNPs by pval threshold 
# the --const-fid avoids the problem of samples that can contain '-' within the name, bacause by default it would be 'family-sample'
plink --vcf ${inpath}/${vcf} --alleleACGT --allow-no-sex --const-fid --make-bed --out ${bed_path}/1000G_ph3_EUR_Diet_${cluster}_${pval}pval