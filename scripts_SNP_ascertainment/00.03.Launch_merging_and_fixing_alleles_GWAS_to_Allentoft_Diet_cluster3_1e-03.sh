#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=gwas-aDNA-merge
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/err/00.03.Launch_merging_and_fixing_alleles_GWAS_to_Allentoft_Diet_cluster3_1e-03_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/out/00.03.Launch_merging_and_fixing_alleles_GWAS_to_Allentoft_Diet_cluster3_1e-03_%j.out
#SBATCH --mem=30000
#SBATCH --time=04:00:00

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last

# Main paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
path_scripts=$(pwd)

# Variables
pval=1e-03
cluster=cluster3
cluster_name_gwas=Cluster_3_TC_prot
gwas_file=${in_path}/Diet_project/merge_GWAS_diet_by_clusters/outer_${pval}/${cluster_name_gwas}/Merged_GWAS_meta_analyses_Diet_${cluster_name_gwas}_outer_pval_${pval}.tsv.gz
vcf_path=${in_path}/Diet_project/aDNA_trend_analysis/Allentoft_VCF_afterQC/outer_merge_pval_${pval}/${cluster}/sorted/
vcf=Diet_SNPs_${cluster}_${pval}_outer_merge_aDNA_Allentoft_afterQC_all_chr_EUR_samples.sorted.vcf.gz
python_dir=${path_scripts}/python_scripts
out_dir=${in_path}/Diet_project/GWAS_matching_aDNA_data_${pval}pval/
filename_out=Merged_GWAS_meta_analyses_Diet_${cluster_name_gwas}_outer_pval_${pval}_common_and_fixed_Allentoft_aDNA

mkdir -p ${out_dir}


# Launch vegeta
python ${python_dir}/Merging_and_fixing_alleles_GWAS_to_Allentoft.py --gwas_data_file ${gwas_file} --in_dir ${vcf_path} --vcf ${vcf} --out_dir ${out_dir} --info_outfile ${filename_out}
