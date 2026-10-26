#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=coords-parser
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_merge/parsers_GWAS_cols
#SBATCH --ntasks=6 
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/parsers_GWAS_cols/err/parse_coords_column_GWAS_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/parsers_GWAS_cols/out/parse_coords_column_GWAS_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=01:30:00
#SBATCH --array=1-3


# WARNINGS!!
# Remember to do following modifications each time it is going to be used:
# 1. Check if the file is compressed to change it by cat or zcat when reading
# 2. Check the file extension to have the correct file name 


# Extract the second column (path) using SLURM_ARRAY_TASK_ID
#path=$(awk -v row=$SLURM_ARRAY_TASK_ID -F '\t' 'NR==row {print $2}' tabulated_file.tsv)

# Paths
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg/ccs
input_dir=${main}/DATA/GWAS_data_diet/Diet_composition/Diet.GWAS.Merino.NatHumBehav.2022/GWAS.SumStats.BMI
out_dir=${input_dir}/parsed_diet_composition

mkdir -p ${out_dir}


# Slurm array variables
array_file=${input_dir}/array_parser_diet_composition.tsv
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file})
filename=$(basename ${file} .txt.gz) # Change file extensions when needed

# Variables intermediate file
header_file=header_${filename}.txt
header_file_ok=header_${filename}_ok.txt


# RUNNING PARSER
# Create header file first for the given file
zcat ${input_dir}/${file} | head -n +1 > ${input_dir}/${header_file}
# Expand the column that is being parsed with the expected output columns
awk 'NR==1 {sub(/MarkerName/, "chromosome position no_important1 no_important2")} 1' ${input_dir}/${header_file} | tr ' ' '\t' | awk 'BEGIN{FS=OFS="\t"}{print "MarkerName", $0}' > ${input_dir}/${header_file_ok}

# One liner parse: chromosome, position, ancestral and effect allele columns from UKB format
cat ${input_dir}/${header_file_ok} <(zcat ${input_dir}/${file} | tail -n +2 | awk 'BEGIN{FS=OFS="\t"}{print $1, $0}' | awk 'BEGIN{FS=OFS="\t"}{gsub(/:/, "\t", $2)} 1') | gzip > ${out_dir}/${filename}_parsed_coords.gz

# Remove header
rm ${input_dir}/${header_file}
rm ${input_dir}/${header_file_ok}