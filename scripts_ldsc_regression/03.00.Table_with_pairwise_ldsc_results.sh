#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression
#SBATCH --job-name=parse-ldsc
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/err/03.00.Table_with_pairwise_ldsc_result_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_ldsc_regression/out/03.00.Table_with_pairwise_ldsc_result_%j.out
#SBATCH --mem=10000
#SBATCH --time=01:00:00


# Libraries


# Variables
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
input_dir=${main}/Diet_project/ldsc_diet/ldsc_regression_Diet
out_dir=${main}/Diet_project/ldsc_diet
outfile=Diet_GWAS_pairwise_ldsc_regression_94_traits.tsv

# List of files with results
cd ${input_dir}
list_res=$(ls *.log)

# Output file
echo "Trait1 Trait2 N_SNPs genetic_corr genetic_corr_se genetic_corr_z genetic_corr_p h2_trait1 h2_trait2 h2_trait1_intercept h2_trait1_intercept_se h2_trait1_ratio mean_chi_trait1 lambda_GC_trait1 h2_trait2_intercept h2_trait2_intercept_se h2_trait2_ratio mean_chi_trait2 lambda_GC_trait2 covariance_intercept covariance_intercept_se" | tr " " "\t" > ${out_dir}/${outfile}


# Iterate over all files except the one we're using
for f in ${list_res}
do  
  echo ${f}

  # Get genetic correlation stats in the following order: corr, se, z and p
  genetic_corr=$(cat ${f} | grep -A1 "gcov_int_se" |  awk '{gsub(/^  */,"")}1' | tr -s " " | tr " " "\t" | cut -f3-6 | tail -n +2)
  # Get heritability estimates, intercept and ratio
  # Phenotype 1
  lambda_gc_1=$(cat ${f} | grep "Lambda" | head -n1 | cut -d " " -f3)
  h2_1=$(cat ${f} | grep "Total Observed scale h2" | cut -d " " -f5 | head -n1)
  intercept1=$(cat ${f} | grep "Intercept" | head -n1 | cut -d " " -f2)
  intercept1_se=$(cat ${f} | grep "Intercept" | cut -d " " -f3 | awk '{gsub(/\(/, "")}1' | awk '{gsub(/\)/, "")}1' | head -n1)
  ratio1=$(cat ${f} | grep "Ratio" | head -n1 | cut -d " " -f2)
  if [[ "$ratio1" =~ ^[[:space:]]*\< ]]
  then
    ratio1=0
  fi

  # Phenotype 2
  lambda_gc_2=$(cat ${f} | grep "Lambda" | tail -n1 | cut -d " " -f3)
  h2_2=$(cat ${f} | grep "Total Observed scale h2" | cut -d " " -f5 | tail -n1)
  intercept2=$(cat ${f} | grep "Intercept" | tail -n2 | head -n1 | cut -d " " -f2)
  intercept2_se=$(cat ${f} | grep "Intercept" | cut -d " " -f3 | awk '{gsub(/\(/, "")}1' | awk '{gsub(/\)/, "")}1' | tail -n2 | head -n1)
  ratio2=$(cat ${f} | grep "Ratio" | tail -n1 | cut -d " " -f2)
  if [[ "$ratio2" =~ ^[[:space:]]*\< ]]
  then
    ratio2=0
  fi

  # Get the names of the genotypes
  trait1=$(echo ${f} | cut -d "." -f1)
  trait2=$(echo ${f} | cut -d "." -f2)

  # Get Number of SNPs used
  N_SNPs=$(cat ${f} | grep "valid" | cut -d " " -f1)

  # Get mean X^2
  mean_chi_1=$(cat ${f} | grep "Chi^2" | cut -d " " -f3 | head -n1)
  mean_chi_2=$(cat ${f} | grep "Chi^2" | cut -d " " -f3 | tail -n1)

  # Get Bivariate relevant statistics
  intercept_cov_stats=$(cat ${f} | grep -A1 "gcov_int_se" |  awk '{gsub(/^  */,"")}1' | tr -s " " | tr " " "\t" | cut -f11-12 | tail -n +2)

  # Append extracted fields in output file
  echo ${trait1} ${trait2} ${N_SNPs} ${genetic_corr} ${h2_1} ${h2_2} ${intercept1} ${intercept1_se} ${ratio1} ${mean_chi_1} ${lambda_gc_1} ${intercept2} ${intercept2_se} ${ratio2} ${mean_chi_2} ${lambda_gc_2} ${intercept_cov_stats} | tr " " "\t" >> ${out_dir}/${outfile}

  #echo ${intercept1}
  #echo ${intercept2}
done

echo "done"