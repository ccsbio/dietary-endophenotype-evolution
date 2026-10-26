#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=add-AF
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs/err/00.02.Launch_add_EUR_AF_risk_allele_prototypes_traits_significant_Diet_cluster2_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_null_distribution_linked_SNPs/out/00.02.Launch_add_EUR_AF_risk_allele_prototypes_traits_significant_Diet_cluster2_%A_%a.out
#SBATCH --mem=20000
#SBATCH --time=05:00:00
#SBATCH --array=1-2

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Main variables
cluster=cluster2
wd=$(pwd)
in_path=${LUSTRE}/ccs
python_dir=${wd}/python_scripts

# Array variables
array_file=${wd}/arrays_clusters/array_Diet_${cluster}_unique_SNPs_trend_analysis_BH_significant_traits.tsv
SNP_dict=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
out_name="$(echo ${trait} | sed -r 's/_UKB.+//' | sed -r 's/_([0-9]+){2,}_.+//' | sed -r 's/_Nat.+//' | sed -r 's/_GIANT.+//' | sed -r 's/_GLGC.+//' | sed 's/_BioPsych//' | sed -r 's/indpMarkers.+//')" 
echo ${out_name}

# Variables
SNP_dict_path=${HOME}/ccs/Diet_project/scripts_SNP_ascertainment/SNP_dictionaries/${cluster}
vcf_path=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted
vcf_file=${vcf_path}/SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.vcf.gz
tabix=${vcf_file}.tbi
out_dir=${wd}/files_SNPs_sets_with_AF/${cluster}/

mkdir -p ${out_dir}


# Launch vegeta
python ${python_dir}/Add_EUR_AF_from_1000G_ph3_hg19_using_SNP_dict_target_trait.py \
--SNP_dict ${SNP_dict_path}/${SNP_dict} \
--vcf_file ${vcf_file} \
--tabix_file ${tabix} \
--out_dir ${out_dir} \
--info_outfile ${out_name}
