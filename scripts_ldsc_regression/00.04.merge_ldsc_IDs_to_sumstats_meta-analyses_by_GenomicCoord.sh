#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression
#SBATCH --job-name=add-rsID
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/00.04.merge_ldsc_IDs_to_sumstats_meta-analyses_by_GenomicCoord_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/00.04.merge_ldsc_IDs_to_sumstats_meta-analyses_by_GenomicCoord_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=01:00:00
#SBATCH --array=1-39


# Libraries


# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/DATA/GWAS_data_diet
out_dir=${main}/Diet_project/ldsc_diet/ldsc_input/gwas_meta-analyses_parsed_and_with_ldsc_rsID

mkdir -p ${out_dir}


# Slurm array variables
wd=$(pwd)
array_file=${wd}/Example_arrays/array_meta-analyses_Diet_by_GenomicCoords_ok.tsv
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
subdirectory=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f3)
chr_col_number=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f5)
pos_col_number=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f6)
# Check if the file is compressed
compressed_ext=$(echo "${file}"  | grep ".gz")


# rsID from ldsc package 
rsID_file=${wd}/files_rsID_ldsc/ldsc_rsID_with_coords_list.tsv


# Variables intermediate file
header_file=header_${filename}.txt
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


# Checks
echo "${compressed_ext}"
echo "${filename}"
echo "${subdirectory}"
echo "${chr_col_number}"
echo "${pos_col_number}"


# Check if the input file is separated by whitespaces instead of tabulations
# Some GWAS can have windows end of line. It is considered in the pipeline
# If so, change it before merge
if [[ ${compressed_ext} == ${file} ]]
then
    total_lines=$(zcat ${input_dir}/${subdirectory}/${file} | wc -l)
    n_fields_when_whistespace=$(zcat ${input_dir}/${subdirectory}/${file} | tr -d '\r' | grep ' ' | tr ' ' 'e' | awk '{print NF}' | sort | uniq -c | cut -d " " -f2)
    n_lines_whitespace=$(zcat ${input_dir}/${subdirectory}/${file} | tr -d '\r' | grep ' ' | tr ' ' 'e' | awk '{print NF}' | sort | uniq -c | cut -d " " -f1)
    if [[ ${n_lines_whitespace} -eq ${total_lines} ]] && [[ ${n_fields_when_whistespace} -eq 1 ]]
    then
        echo "Is not separeted by tabs. Conversion from whitespace to tab separated is saved as "${input_dir}/${subdirectory}/${filename}.tsv.gz
        zcat ${input_dir}/${subdirectory}/${file} | tr -d '\r' | tr " " "\t" | bgzip > ${input_dir}/${subdirectory}/${filename}.tsv.gz
        file_parsed=${input_dir}/${subdirectory}/${filename}.tsv.gz
        compressed_ext=${file} # we compressed the file, so this will make the checkpoint to read it as if compressed
    else
        file_parsed=no
    fi
else
    total_lines=$(cat ${input_dir}/${subdirectory}/${file} | wc -l)
    n_fields_when_whistespace=$(cat ${input_dir}/${subdirectory}/${file} | tr -d '\r' | grep ' ' | tr ' ' 'e' | awk '{print NF}' | sort | uniq -c | cut -d " " -f2)
    n_lines_whitespace=$(cat ${input_dir}/${subdirectory}/${file} | tr -d '\r' | grep ' ' | tr ' ' 'e' | awk '{print NF}' | sort | uniq -c | cut -d " " -f1)
    if [[ ${n_lines_whitespace} -eq ${total_lines} ]] && [[ ${n_fields_when_whistespace} -eq 1 ]]
    then
        echo "Is not separeted by tabs. Conversion from whitespace to tab separated is saved as "${input_dir}/${subdirectory}/${filename}.tsv.gz
        cat ${input_dir}/${subdirectory}/${file} | tr -d '\r' | tr " " "\t" | bgzip > ${input_dir}/${subdirectory}/${filename}.tsv.gz
        file_parsed=${input_dir}/${subdirectory}/${filename}.tsv.gz
        compressed_ext=${file} # we compressed the file, so this will make the checkpoint to read it as if compressed
    else
        file_parsed=no
    fi
fi


# Define if the input is the original file or the parsed file
if [[ ${file_parsed} == "no" ]]
then
    file_to_read=${input_dir}/${subdirectory}/${file}
else
    file_to_read=${input_dir}/${subdirectory}/${filename}.tsv.gz
fi



# RUNNING MERGING
# 1. Create header file first for the given file
if [[ ${compressed_ext} == ${file} ]]
then
    zcat ${file_to_read} | head -n +1 | tr -d '\r' | awk 'BEGIN{FS=OFS="\t"}{print "chr:pos", $0, "ref:alt_merging_file", "ldsc_rsID"}' > ${input_dir}/${subdirectory}/${header_file} # add colname for rsID 
else
    cat ${file_to_read} | head -n +1 | tr -d '\r' | awk 'BEGIN{FS=OFS="\t"}{print "chr:pos", $0, "ref:alt_merging_file", "ldsc_rsID"}' > ${input_dir}/${subdirectory}/${header_file} # add colname for rsID 
fi

# 2. Merging files by a new ID (which is chr:pos) to get rsID that match the ldsc package
# Some GWAS can have windows end of line. It is considered in the pipeline
# It is very important to fill missing data of every cell with NA. Otherwise when join is performed, some data can be located in the wrong columns
if [[ ${compressed_ext} == ${file} ]]
then
    echo "File compressed"
    join <(zcat ${file_to_read} | tr -d '\r' | awk -v chr=${chr_col_number} -v pos=${pos_col_number} 'BEGIN{FS=OFS="\t"}{print $chr":"$pos, $0}' |  awk 'BEGIN{FS=OFS="\t"}{for(i=1;i<=NF;i++) if($i=="") $i="NA"; print}' | sort -T ${tmp_sort} -k1) <(cat ${rsID_file} | sort -T ${tmp_sort} -k1) -1 1 -2 1 | tr " " "\t" > ${tmp_sort}/${intermediate_file}
else
    echo "File not compressed"
    join <(cat ${file_to_read} | tr -d '\r' | awk -v chr=${chr_col_number} -v pos=${pos_col_number} 'BEGIN{FS=OFS="\t"}{print $chr":"$pos, $0}' |  awk 'BEGIN{FS=OFS="\t"}{for(i=1;i<=NF;i++) if($i=="") $i="NA"; print}' | sort -T ${tmp_sort} -k1) <(cat ${rsID_file} | sort -T ${tmp_sort} -k1) -1 1 -2 1 | tr " " "\t" > ${tmp_sort}/${intermediate_file}
fi
# Add the header
cat ${input_dir}/${subdirectory}/${header_file} ${tmp_sort}/${intermediate_file} | bgzip > ${out_dir}/${filename}_rsID_merged.tsv.gz # Add header


# Remove intermediate files
rm ${input_dir}/${subdirectory}/${header_file}
#rm ${tm}/${intermediate_file}
#rm -rf ${tmp_sort}
