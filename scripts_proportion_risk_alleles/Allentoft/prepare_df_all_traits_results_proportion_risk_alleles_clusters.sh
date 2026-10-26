#!/bin/bash
# Stores .tsv files to process
# Double parenthesis so it becomes an array
cd ${1}
files_to_iterate=($( ls Proportion_* ))

# Get the length of the arrays
num_elements=${#files_to_iterate[@]}
echo ${num_elements}

# Create outpath for intermediate files and output
outpath="./processed_df"
mkdir ${outpath}

# Create the file with the first line (the header)
echo "Sample_ID Proportion_risk_allele N_non_missing N_missing N_markers Prop_missing Trait Trait_label" | tr " " "\t" > ${outpath}/header.txt
# Prepare command file
echo "cat ${outpath}/header.txt" > ${outpath}/command.txt


# Now iterate over files
for ((i = 0; i < num_elements; i++))
do  
  f="${files_to_iterate[i]}"
  filename="$(basename ${f} .tsv)"
  trait="$(echo ${f} | sed 's/^Proportion_risk_allele_//' | sed -r 's/_([0-9]+){2,}_.+//' |  sed -r 's/_Nat.+//' | sed 's/_BioPsych//' | sed -r 's/_GIANT.+//')"
  echo ${trait}
  tail -n+2 ${f} | awk 'BEGIN{FS=OFS="\t"}{print $1, $2, $3, $4}' | awk 'BEGIN{FS=OFS="\t"}{print $0, $3+$4}' |\
   awk -v t=${trait} 'BEGIN{FS=OFS="\t"}{print $0, $4/$5, t, t" n="$5}' > ${outpath}/${filename}_processed.tsv
  echo "${outpath}/${filename}_processed.tsv" >> ${outpath}/command.txt
done  


# Concat processed files
commands=$(echo "$(cat ${outpath}/command.txt | tr "\n" " ")")
echo ${commands}
eval ${commands} > ${outpath}/Proportion_risk_allele_Diet_cluster${2}.tsv

# Remove
rm ${outpath}/command.txt
rm ${outpath}/header.txt
