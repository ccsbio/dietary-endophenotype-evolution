#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G
#SBATCH --job-name=vcf-filt
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/err/01.01.Subset_chromosomes_1-22_1000G_phase3_hg19_Diet_cluster1_with_GT_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_1000G/out/01.01.Subset_chromosomes_1-22_1000G_phase3_hg19_Diet_cluster1_with_GT_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules


# CESGA 
module load bedtools/2.30.0


# Paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
wd=$(pwd)
pval=1e-03
cluster=cluster1
vcf_path=${in_path}/DATA/1000G_human_variation_phase3/VCFs_1000G_phase3_hg19
bed_path=${wd}/SNP_files/outer_merge_${pval}
out_path=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}
bedfile=SNPs_diet_${cluster}_outer_merge_${pval}pval.bed
#contigs_header=${in_path}/Diet_project/scripts_subset_vcf_archaics/contigs_header_VCF_human_genome_hg19.txt

mkdir -p ${out_path}


# Subset by chromosome 
# If you want the vcf header add -header before -wa, which just prints the results
touch ${out_path}/vcf_filtering_name_list.txt
#for i in 2 14 19
for i in {1..22}
do
    #cat <(zcat ${vcf_path}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | grep "##") \
    bedtools intersect -a ${vcf_path}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -b ${bed_path}/${bedfile} -header -wa | \
    bgzip > ${out_path}/SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_chr${i}.vcf.gz
    echo "SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_chr${i}.vcf.gz" >> ${out_path}/vcf_filtering_name_list.txt
done