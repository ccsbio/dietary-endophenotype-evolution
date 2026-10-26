#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_GA_clustering
#SBATCH --job-name=GA-clust
#SBATCH --cpus-per-task=20
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_GA_clustering/err/02.03.Elitism.launch_GA_clustering_KMEDOIDS_traits_Diet_100iter_400pop_0_2MR_K6_4dim_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_GA_clustering/out/02.03.Elitism.launch_GA_clustering_KMEDOIDS_traits_Diet_100iter_400pop_0_2MR_K6_4dim_%j.out
#SBATCH --mem=10000
#SBATCH --time=03:00:00

# Load R


# CESGA
module load cesga/2022 gcc/system
module load R/4.3.3


# Main paths
in_path=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
wd=$(pwd)
# Variables
file_MDS_from_LDSC=${wd}/data/LDSC_MDS_matrix_Diet_58_GWAS_interval_4dim_BMIUKB.tsv
out_dir=${in_path}/Diet_project/GA_clustering_traits_Diet/using_kmedoids/dimensions_4/Elitism/

mkdir -p ${out_dir}

# GA args
pop_size=400
mutation_rate=0.2
iter=100
rep_type="Elitism"
cores=16
k_clusters=6

# Launch GA
Rscript ${wd}/R_scripts/clustering_kmedoids/Run_GA_clustering_traits_Diet_by_kmedoids.R ${pop_size} ${mutation_rate} ${iter} ${rep_type} ${k_clusters} ${file_MDS_from_LDSC} ${out_dir} ${cores}
