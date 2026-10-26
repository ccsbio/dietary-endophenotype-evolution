#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression
#SBATCH --job-name=liftOver
#SBATCH --cpus-per-task=6
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/Build_conversion_from_hg38_to_hg19_unprocessed_GWAS_Diet_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/Build_conversion_from_hg38_to_hg19_unprocessed_GWAS_Diet_%A_%a.out
#SBATCH --mem=50000
#SBATCH --time=04:00:00
#SBATCH --array=1-2


# Libraries
# Museo
#module load userspace
#module load python/conda
#source activate vegeta_last # Activate vegeta environment

# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Variables
packages_path=/home/csic/gpo/olg/ccs/packages
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/DATA/GWAS_data_diet


# Slurm array variables
wd=$(pwd)
#echo ${wd}
array_file=${wd}/Example_arrays/array_meta-analyses_Diet_build_conversion_hg38_to_hg19.tsv
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
subdirectory=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
chr_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f3)
pos_col=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f4)
python_dir=${wd}/python_scripts
liftOver_dir=${packages_path}/vegeta_last_version/liftOver_linux/

# Run liftOver conversion
cd ${input_dir}/${subdirectory}
python ${python_dir}/Fast_build_conversion_GWAS_sumstats.py --file ${file} --dir ${input_dir}/${subdirectory} --chr_col ${chr_col} --pos_col ${pos_col} --original_build hg38 --directory_liftOver ${liftOver_dir}



