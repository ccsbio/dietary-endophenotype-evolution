#!/bin/bash

# Store arguments
file_array=$1
list_traits=$2
file_name=$(basename ${file_array} .tsv)


# Filter
awk 'BEGIN{FS=OFS="\t"} NR==FNR {traits[$1]; next} ($2 in traits)' ${list_traits} ${file_array} > ${file_name}_BH_significant_traits.tsv