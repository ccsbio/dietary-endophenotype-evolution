#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G
#SBATCH --job-name=1000G-EUR    
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/err/04.02.Filter_EUR_population_1000G_ALL_chr_filtered_VCF_SNPs_Diet_cluster2_with_GT_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/out/04.02.Filter_EUR_population_1000G_ALL_chr_filtered_VCF_SNPs_Diet_cluster2_with_GT_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules
# CESGA
module load bcftools/1.19




# Paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
cluster=cluster2
vcf_path=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted
vcf=SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.vcf.gz
sample_file=${in_path}/DATA/1000G_human_variation_phase3/VCFs_1000G_phase3_hg19/EUR.samples.txt
out_path=${vcf_path}/European_populations

mkdir -p ${out_path}

# Fix ALT column and INFO field to match genotype data only in EUR populations
#bcftools view -S ${sample_file} -Ou ${vcf_path}/${vcf} \
#| bcftools view --trim-alt-alleles -Ou \
#| bcftools +fill-tags -Ou -- -t AC,AN,AF,NS,DP \
#| bcftools annotate -x 'INFO/.+_AF' -Ou \
#| bcftools view -i 'AC>0' -Oz -o ${vcf_path}/ASD_ALL_risk_SNPs_${pval}_1000G_ph3.EUR.recalc.vcf.gz

#| bcftools annotate -x 'INFO/.+_AF,INFO/MULTI_ALLELIC' -Ou \
bcftools view -S ${sample_file} -Ou ${vcf_path}/${vcf} \
| bcftools view --trim-alt-alleles -Ou \
| bcftools +fill-tags -Ou -- -t AC,AN,AF,NS \
| bcftools view -i 'AC>0' -Ou \
| bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz \
-o ${out_path}/SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.EUR.recalc.vcf.gz