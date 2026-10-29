#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset
#SBATCH --job-name=vcf-filt-samples
#SBATCH --error=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/err/02.01.Filtering_sample_genos_cols_per_chr_1000G_SNPs_linked_to_protein_coding_genes_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/out/02.01.Filtering_sample_genos_cols_per_chr_1000G_SNPs_linked_to_protein_coding_genes_%j.out
#SBATCH --mem=5000
#SBATCH --time=04:00:00

# Load modules



# Paths
in_path=${LUSTRE}/ccs
vcf_path=${in_path}/SNPs_linked_to_genes_Allentoft/1000G_data_final_set_protein_coding
out_path=${vcf_path}/VCF_filt_9_cols

mkdir -p ${out_path}

# Subset by chromosome 
# If you want the vcf header add -header before -wa, which just prints the results
for i in {1..22}
do
   vcf_name=SNPs_linked_to_protein_coding_genes_chr${i}
   if [[ $i -eq 1 ]]
   then
      # VCF of chr1 has the header, it must be also adjuster to 9 columns
      zcat ${vcf_path}/${vcf_name}.vcf.gz | grep "#" | tail -n1 | cut -f1-9 > ${out_path}/header_vcf_9fields.txt
      zcat ${vcf_path}/${vcf_name}.vcf.gz | grep -v "#" | cut -f1-9 > ${out_path}/${vcf_name}_9fields.vcf
      cat ${out_path}/header_vcf_9fields.txt ${out_path}/${vcf_name}_9fields.vcf | bgzip > ${out_path}/${vcf_name}_9fields.vcf.gz
      rm ${out_path}/${vcf_name}_9fields.vcf # Remove intermediate file
   else
      n_lines=$(zcat ${vcf_path}/${vcf_name}.vcf.gz | grep -v "#" | wc -l | cut -d " " -f1)
      echo ${n_lines}
      if [[ ${n_lines} -gt 0 ]]
      then
         zcat ${vcf_path}/${vcf_name}.vcf.gz | grep -v "#" | cut -f1-9 | bgzip > ${out_path}/${vcf_name}_9fields.vcf.gz
      fi
   fi
done
