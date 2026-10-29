#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_GA_geo_trend
#SBATCH --job-name=GA-geo
#SBATCH --cpus-per-task=20
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_GA_geo_trend/err/01.01.Tournament.launch_GA_geographic_trend_split_Diet_cluster1_100iter_200pop_0_2MR_5tournament_size_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_GA_geo_trend/out/01.01.Tournament.launch_GA_geographic_trend_split_Diet_cluster1_100iter_200pop_0_2MR_5tournament_size_%A_%a.out
#SBATCH --mem=20000
#SBATCH --time=05:00:00
#SBATCH --array=1-3

# Load R


# CESGA
module load cesga/2022 gcc/system
module load R/4.3.3


# Main paths
in_path=${LUSTRE}/ccs
wd=$(pwd)

# Array variables
cluster=cluster1
array_file=${wd}/arrays_clusters/array_Diet_${cluster}_unique_SNPs_trend_analysis_BH_significant_traits.tsv
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
trait_clean="$(echo ${trait} | sed -r 's/_UKB.+//' | sed -r 's/_([0-9]+){2,}_.+//' | sed -r 's/_Nat.+//' | sed -r 's/_GIANT.+//' | sed -r 's/_GLGC.+//' | sed 's/_BioPsych//' | sed -r 's/indpMarkers.+//')"
n_prot="$(echo "${trait}" | tr "_" "\n" | grep "indpMarkers" | sed -r 's/indpMarkers//')"
echo ${trait_clean}

# Variables
samples_with_metadata=${wd}/data/Diet_all_traits_${cluster}_classified_Allentoft_temporal_trends_400-10yearsBP_results.tsv
out_dir=${in_path}/Diet_project/GA_geo_trend_results_Allentoft/${cluster}/

mkdir -p ${out_dir}

# GA args
pop_size=200
mutation_rate=0.2
iter=100
rep_type="Tournament"
cores=16
x_variable_name="Age.average"
y_variable_name="Proportion_risk_allele"
long="Longitude.ok"
lat="Latitude.ok"
tournament_size=5

# Launch GA
Rscript ${wd}/R_scripts/Run_GA_geographic_trends_split_Diet_by_trait.R \
${pop_size} \
${mutation_rate} \
${iter} \
${rep_type} \
${samples_with_metadata} \
${x_variable_name} \
${y_variable_name} \
${long} \
${lat} \
${out_dir} \
${trait_clean} \
${cores} \
${tournament_size}
