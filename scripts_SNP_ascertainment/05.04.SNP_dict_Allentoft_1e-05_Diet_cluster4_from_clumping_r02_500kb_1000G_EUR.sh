#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=SNP-dict
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/err/05.04.SNP_dict_Allentoft_1e-05_Diet_cluster4_from_clumping_r02_500kb_1000G_EUR_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_SNP_ascertainment/out/05.04.SNP_dict_Allentoft_1e-05_Diet_cluster4_from_clumping_r02_500kb_1000G_EUR_%A_%a.out
#SBATCH --mem=10000
#SBATCH --time=02:00:00
#SBATCH --array=1-21

# Load conda
# CESGA
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last





# Array and global variables
main=/mnt/lustre/scratch/nlsas/home/csic/gpo/olg/ccs
path_scripts=$(pwd)
pval=1e-03
ld_threshold=0.2
cluster=cluster4
array_file=${path_scripts}/arrays/array_traits_${cluster}.txt
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)

# This is an optional pval threshold, smaller than the one used to create intermediate files.
# Nevertheless, it should be the same one as the one used for the clump. Otherwise, maybe SNPs are not found. 
# It can be an input with the argument --pval_threshold
pval_threshold=1e-05 


# Variables
cluster_name_gwas=Cluster_4_Trunk_dat_free_prot
gwas_path=${main}/Diet_project/GWAS_matching_aDNA_data_${pval}pval/
gwas_file=${gwas_path}Merged_GWAS_meta_analyses_Diet_${cluster_name_gwas}_outer_pval_${pval}_common_and_fixed_Allentoft_aDNA_update1000G.tsv.gz
plink_clumped=${path_scripts}/clumping/pval${pval_threshold}_LD${ld_threshold}_1000G_EUR_pop/${cluster}/1000G_EUR_Diet_${cluster}_${trait}_${pval_threshold}pval_${ld_threshold}r2_500kb_parsed.clumped
python_dir=${path_scripts}/python_scripts
out_dir=${path_scripts}/SNP_dictionaries/${cluster}/
log_file=${out_dir}/log_${cluster}_ascertainment_plink_clumped.txt
filename_out=clumping_r${ld_threshold}_500kb_1000G_EUR
qc_file=${path_scripts}/complex_LD_regions_hg19_chr.tsv

mkdir -p ${out_dir}


# Initializate log
if [ ! -f ${log_file} ]
then
    echo "trait n_snps_uniq_associated n_items_ascertained n_indep_ascertained pval_threshold" | tr " " "\t" > ${log_file}
fi


# Launch vegeta
# --pval_threshold ${pval_threshold}
python ${python_dir}/SNP_ascertainment_from_PLINK_clumping_and_fixed_GWAS_with_Allentoft_aDNA_by_trait_from_merge.py --gwas_data_file ${gwas_file} --plink_clumped_file ${plink_clumped} --pval_threshold ${pval_threshold} --phenotype_name ${trait} --log_filename ${log_file} --out_dir ${out_dir} --info_outfile ${filename_out} --regions_qc ${qc_file}
