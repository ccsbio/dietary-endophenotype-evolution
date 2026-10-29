#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft
#SBATCH --job-name=bedfile-clust
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/err/00.04.Create_bedfile_Diet_cluster4_outer_merge_1e-03pval_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_subset_aDNA/subset_Allentoft/out/00.04.Create_bedfile_Diet_cluster4_outer_merge_1e-03pval_%j.out
#SBATCH --mem=20000
#SBATCH --time=04:00:00




# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Paths
pval=1e-03
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
scripts_path=$(pwd)
python_dir=${scripts_path}/python_scripts
merged_gwas=${main}/Diet_project/merge_GWAS_diet_by_clusters/outer_${pval}/Cluster_4_Trunk_fat_free_prot/Merged_GWAS_meta_analyses_Diet_Cluster_4_Trunk_fat_free_prot_outer_pval_${pval}.tsv.gz
out_path=${scripts_path}/SNP_files/outer_merge_${pval}/
filename_out=SNPs_diet_cluster4_outer_merge_${pval}pval

mkdir -p ${out_path}


# Run create bedfile with genomic positions from the whole cluster
python ${python_dir}/01.Create_bedfile_outer_merge_meta_analyses_Diet.py --merged_gwas ${merged_gwas} --out_dir ${out_path} --info_outfile ${filename_out}
