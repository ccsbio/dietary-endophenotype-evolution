#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=vegeta-diet
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_merge
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/err/01.04.Launch_outer_merge_vegeta_cluster_4_Trunk_fat_free_prot_pval_1e-03_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/out/01.04.Launch_outer_merge_vegeta_cluster_4_Trunk_fat_free_prot_pval_1e-03_%j.out
#SBATCH --mem=50000
#SBATCH --time=09:00:00




# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Variables
main_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs/
scripts_path=$(pwd)
pval_threshold=1e-03
out_dir=${main_path}/Diet_project/merge_GWAS_diet_by_clusters/outer_${pval_threshold}/Cluster_4_Trunk_fat_free_prot/
python_dir=${scripts_path}/python_scripts
pheno=array_cluster4_Trunk_fat_free_prototype.tsv
array_path=${scripts_path}/arrays/

mkdir -p ${out_dir}
cp ${array_path}${pheno} ${out_dir}${pheno}

# Launch vegeta
python ${python_dir}/01.04.Cluster_4.Outer_merge_meta_analyses_Diet_pval_${pval_threshold}.py --array_file ${pheno} --out_dir ${out_dir} --pval_threshold ${pval_threshold}