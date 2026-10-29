#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression
#SBATCH --job-name=add-rsID
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/00.05.merge_ldsc_IDs_to_sumstats_UKB_Ant_by_GenomicCoord_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/00.05.merge_ldsc_IDs_to_sumstats_UKB_Ant_by_GenomicCoord_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=02:00:00
#SBATCH --array=1-26


# Libraries


# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/DATA/GWAS_data_diet/UKB/UKB_Anthropometric_pheno
out_dir=${main}/Diet_project/ldsc_diet/ldsc_input/UKB_parsed_and_with_ldsc_rsID/UKB_Anthropometric

mkdir -p ${out_dir}


# Slurm array variables
wd=$(pwd)
array_file=${input_dir}/array_UKB_Anthropometric.txt
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file})
filename=$(basename ${file} .tsv.bgz)
rsID_file=${wd}/files_rsID_ldsc/ldsc_rsID_with_coords_list.tsv


# Variables intermediate file
header_file=header_${filename}.txt
new_header_file=new_header_${filename}.txt
intermediate_file=intermediate_${filename}.tsv
#tmp_sort=${out_dir}/tmp_sort_${filename}
# ${VAR:-DEFAULT} is parameter expansion: use $VAR if it’s set and non-empty, otherwise use DEFAULT.
tmp_sort=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/sort_${SLURM_ARRAY_TASK_ID}.XXXX")
# Conditional statements in the trap command protects you from the rare case where $tmp_sort is unset or accidentally empty (which could otherwise turn rm -rf "$tmp_sort" into rm -rf /).
# With the trap in place you don’t need any further manual rm -rf "$tmp_sort"—and you’re guaranteed that the scratch directory disappears no matter how the script exits.
trap 'd=${tmp_sort:-}; [[ -n $d && -d $d ]] && { rm -rf "$d" && echo "Removed $d"; }' EXIT
echo ${tmp_sort}
if [[ -n $tmp_sort && -d $tmp_sort ]]
then
    echo "Exists and non-empty"
fi
#mkdir -p ${tmp_sort}


# RUNNING MERGING
# 1. Create header file first for the given file
zcat ${input_dir}/${file} | head -n +1 > ${input_dir}/${header_file}
# 2. Expand column named variant which will transform in 4: chromosome, position, reference and alternative. Finally Let a copy in the first column with the original ID from the UKB and the last column will be for the rsID
awk 'NR==1 {sub(/variant/, "chromosome position ref alt")} 1' ${input_dir}/${header_file} | tr " " "\t" | awk 'BEGIN{FS=OFS="\t"}{print "chr:pos", $0, "alleles_merged_ldsc", "rsID"}' > ${input_dir}/${new_header_file} # add colname for rsID


# 3. Merging files by UKB ID to get rsID that match the ldsc package
join <(zcat ${input_dir}/${file} | awk 'BEGIN{FS=OFS="\t"}{gsub(/:/, "\t", $1)} 1' | awk 'BEGIN{FS=OFS="\t"}{print $1":"$2, $0}' | sort -T ${tmp_sort} -k1) <(cat ${rsID_file} | sort -T ${tmp_sort} -k1) -1 1 -2 1 | tr " " "\t" > ${tmp_sort}/${intermediate_file}
cat ${input_dir}/${new_header_file} ${tmp_sort}/${intermediate_file} | bgzip > ${out_dir}/${filename}_rsID_merged.tsv.gz # Add header


# Remove intermediate files
rm ${input_dir}/${header_file}
rm ${input_dir}/${new_header_file}
#rm ${out_dir}/${intermediate_file}
#rm -rf ${tmp_sort}