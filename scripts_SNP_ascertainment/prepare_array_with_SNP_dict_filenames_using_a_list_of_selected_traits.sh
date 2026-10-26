#!/bin/bash
directory=$1
list_substrings_file=$2
cluster=$3
pval=$4
find ${directory} -type f | grep -Ff ${list_substrings_file} | tr "/" "\n" | grep ".pkl" | awk 'BEGIN{FS=OFS="\t"}{print $1, $1}' | \
awk -v pval="$pval" 'BEGIN{FS=OFS="\t"}{gsub("dict_Allentoft_SNP_IDs_" pval "pval_", "", $2)}1' | \
awk 'BEGIN{FS=OFS="\t"}{gsub(/_clumping_r0.2_500kb_1000G_EUR.pkl/, "", $2)}1' > array_Diet_${cluster}_unique_SNPs_trend_analysis.tsv
