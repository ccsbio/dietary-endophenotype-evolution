#!/bin/bash
log_file=$1
cluster=$2
awk 'BEGIN{FS=OFS="\t"} {if ($4 > 20) print}' ${log_file} | cut -f1 | tail -n +2 > list_traits_Diet_${cluster}_unique_SNPs_trend_analysis.tsv
