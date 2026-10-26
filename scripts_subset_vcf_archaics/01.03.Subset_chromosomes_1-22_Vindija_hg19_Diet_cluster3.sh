#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics
#SBATCH --job-name=vcf-filt
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics/err/01.03.Subset_chromosomes_1-22_Vindija_hg19_Diet_cluster3_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_vcf_archaics/out/01.03.Subset_chromosomes_1-22_Vindija_hg19_Diet_cluster3_%j.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00

# Load modules
# CESGA 
module load bedtools/2.30.0
module load bcftools/1.19


# Paths
in_path=${LUSTRE}/ccs
wd=$(pwd)
cluster=cluster3
vcf_path=${in_path}/DATA/Archaic_genomes/Neanderthal/Vindija
bed_path=${wd}/bedfiles/outer_merge_1e-03
out_path=${in_path}/Diet_project/Vindija_VCFs_Diet_clusters/${cluster}
bedfile=SNPs_diet_${cluster}_outer_merge_1e-03pval.bed
contigs_header=${wd}/contigs_header_VCF_human_genome_hg19.txt

mkdir -p ${out_path}


# Subset by chromosome 
# If you want the vcf header add -header before -wa, which just prints the results
touch ${out_path}/vcf_filtering_name_list.txt
for i in {1..22}
do
    cat <(bcftools view -h ${vcf_path}/chr${i}_mq25_mapab100.vcf.gz | grep ^##) <(cat ${contigs_header}) <(bcftools view -h ${vcf_path}/chr${i}_mq25_mapab100.vcf.gz | grep -v ^##) > ${out_path}/header_chr${i}.txt
    bcftools reheader --header ${out_path}/header_chr${i}.txt ${vcf_path}/chr${i}_mq25_mapab100.vcf.gz | \
    bedtools intersect -a stdin -b ${bed_path}/${bedfile} -header -wa | bgzip > ${out_path}/Vindija_SNPs_Diet_${cluster}_chr${i}.vcf.gz
    echo "Vindija_SNPs_Diet_${cluster}_chr${i}.vcf.gz" >> ${out_path}/vcf_filtering_name_list.txt
    rm ${out_path}/header_chr${i}.txt
done