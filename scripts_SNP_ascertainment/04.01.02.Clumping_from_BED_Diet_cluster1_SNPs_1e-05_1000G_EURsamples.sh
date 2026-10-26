#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=clumping-1000G
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/err/04.01.02.Clumping_from_BED_Diet_cluster1_SNPs_1e-05_1000G_EURsamples_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/out/04.01.02.Clumping_from_BED_Diet_cluster1_SNPs_1e-05_1000G_EURsamples_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=04:00:00
#SBATCH --array=1-9


# Load modules
# Cesga
module load plink/1.9b5



# WARNING!
# The p.value of the .assoc file must be the same one as the clump, since we look at most significant SNPs UNIQUELY associated to a trait within a genetic cluster
# If you use a higher p.value in the .assoc file rather than the one selected for clumping, then you will loss a lot of relevant SNPs because they'll be more shared across the cluster

# Array and global variables
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg/ccs
wd=$(pwd)
pval=1e-03
cluster=cluster1
array_file=${wd}/arrays/array_traits_${cluster}.txt
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# Variables
in_path=${main}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted/European_populations
BED=${in_path}/BED_Diet_${cluster}_SNPs_${pval}_1000G_EUR/1000G_ph3_EUR_Diet_${cluster}_${pval}pval
ld_threshold=0.2
pval_clump=1e-05
assoc_file=${wd}/assoc_files/${cluster}/Diet_${cluster}_${trait}_SNPs_${pval_clump}_1000G_ph3_merge_GWAS_and_aDNA.assoc
outpath=${wd}/clumping/pval${pval_clump}_LD${ld_threshold}_1000G_EUR_pop/${cluster}
filename_out=1000G_EUR_Diet_${cluster}_${trait}_${pval_clump}pval_${ld_threshold}r2_500kb

mkdir -p ${outpath}


# Clump index SNPs (smaller pvals) and prune ld SNPs within 500kb
#      --keep ${sample_file_plink} \
plink --bfile ${BED} \
      --geno 0.8 \
      --clump ${assoc_file} \
      --clump-p1 ${pval_clump} \
      --clump-p2 ${pval_clump} \
      --clump-r2 ${ld_threshold} \
      --clump-kb 500 \
      --out ${outpath}/${filename_out}


# Parse clumped output
cat ${outpath}/${filename_out}.clumped | tr -s " " "\t " > ${outpath}/${filename_out}_parsed.clumped
rm ${outpath}/${filename_out}.clumped