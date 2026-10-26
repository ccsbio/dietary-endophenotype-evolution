#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment
#SBATCH --job-name=assoc-file
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/err/03.03.02.Launch_create_assoc_file_by_trait_Diet_SNPs_cluster3_1e-05_with_IDs_1000G_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/out/03.03.02.Launch_create_assoc_file_by_trait_Diet_SNPs_cluster3_1e-05_with_IDs_1000G_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=02:00:00
#SBATCH --array=1-7

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Array and global variables
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg/ccs
wd=$(pwd)
pval=1e-03
pval_threshold=1e-05
cluster=cluster3
array_file=${wd}/arrays/array_traits_${cluster}.txt
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# Variables
cluster_name_gwas=Cluster_3_TC_prot
gwas_path=${main}/Diet_project/GWAS_matching_aDNA_data_${pval}pval/
gwas_file=${gwas_path}Merged_GWAS_meta_analyses_Diet_${cluster_name_gwas}_outer_pval_${pval}_common_and_fixed_Allentoft_aDNA_update1000G.tsv.gz
out_path=${wd}/assoc_files/${cluster}/
python_dir=${wd}/python_scripts
filename_out=Diet_${cluster}_${trait}_SNPs_${pval_threshold}_1000G_ph3_merge_GWAS_and_aDNA
log_file=${out_path}/log_${cluster}_assoc_files.txt

mkdir -p ${out_path}


# Initializate log
if [ ! -f ${log_file} ]
then
    echo "trait n_snps_uniq_associated pval_threshold" | tr " " "\t" > ${log_file}
fi


# Launch vegeta
python ${python_dir}/Create_assoc_file_merged_GWAS_by_trait_for_clumping_with_plink_and_1000G_data.py --gwas_data_file ${gwas_file} --phenotype_name ${trait} --pval_threshold ${pval_threshold} --log_filename ${log_file} --out_dir ${out_path} --info_outfile ${filename_out} 
