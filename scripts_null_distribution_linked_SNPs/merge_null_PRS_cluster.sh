#!/bin/bash

# first change to input path
cd $1

# Stores .tsv files to process
# Double parenthesis so it becomes an array
files_to_iterate=($( ls Null_* ))
# store also cluster name
cluster=$2 

# Get the length of the arrays
num_elements=${#files_to_iterate[@]}
echo ${num_elements}

# Create outpath for intermediate files and output
outpath="./processed_df"
mkdir ${outpath}

# Create the file with the first line (the header)
head -n1 "${files_to_iterate[0]}" | awk 'BEGIN{FS=OFS="\t"}{print "Trait", $0}' > ${outpath}/header.txt
# Prepare command file
echo "cat ${outpath}/header.txt" > ${outpath}/command.txt


# Now iterate over files
for ((i = 0; i < num_elements; i++))
do  
  f="${files_to_iterate[i]}"
  filename="$(basename ${f} .tsv)"
  trait="$(echo ${f} | tr "_" "\n" | tail -n +14 | tr "\n" "_" | cut -d "." -f1)"
  echo ${trait}
  tail -n+2 ${f} | awk -v t=${trait} 'BEGIN{FS=OFS="\t"}{print t, $0}' > ${outpath}/${filename}_processed.tsv
  echo "${outpath}/${filename}_processed.tsv" >> ${outpath}/command.txt
done  

# Concat processed files
commands=$(echo "$(cat ${outpath}/command.txt | tr "\n" " ")")
echo ${commands}
eval ${commands} > ${outpath}/Null_PRS_1000replicates_allentoft_SNPs_50Kb_togenes_${cluster}.tsv

# Remove
rm ${outpath}/command.txt
rm ${outpath}/header.txt
