#!/bin/bash
out_dir=${LUSTRE}/ccs/DATA/GENCODE/hg19_release_47
mkdir -p ${out_dir}

# Download comprehensive gene annotation in reference chromosomes of gh38 mapped to gh37
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gff3.gz -P ${out_dir}
