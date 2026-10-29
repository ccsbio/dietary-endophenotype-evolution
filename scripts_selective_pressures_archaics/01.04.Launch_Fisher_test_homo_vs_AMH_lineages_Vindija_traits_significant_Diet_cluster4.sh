#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Fisher-archaics
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_selective_pressures_archaics
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_selective_pressures_archaics/err/01.04.Launch_Fisher_test_homo_vs_AMH_lineages_Vindija_traits_significant_Diet_cluster4_%A_%a.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_selective_pressures_archaics/out/01.04.Launch_Fisher_test_homo_vs_AMH_lineages_Vindija_traits_significant_Diet_cluster4_%A_%a.out
#SBATCH --mem=15000
#SBATCH --time=04:00:00
#SBATCH --array=1-20

# Load conda



# Cesga
module --force purge
module load cesga/system miniconda3/22.11.1-1
conda activate vegeta_last


# Main paths
in_path=${LUSTRE}/ccs
wd=$(pwd)
python_dir=${wd}/python_scripts


# Array variables
cluster=cluster4
array_file=${wd}/arrays_clusters/array_Diet_${cluster}_unique_SNPs_trend_analysis.tsv
SNP_dict=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f1)
trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${array_file} | cut -f2)
trait_clean="$(echo ${trait} | sed -r 's/_UKB.+//' | sed -r 's/_([0-9]+){2,}_.+//' | sed -r 's/_Nat.+//' | sed -r 's/_GIANT.+//' | sed -r 's/_GLGC.+//' | sed 's/_BioPsych//' | sed -r 's/indpMarkers.+//')"
n_prot="$(echo "${trait}" | tr "_" "\n" | grep "indpMarkers" | sed -r 's/indpMarkers//')"
trait_label_merged_GWAS="$(echo "${trait}" | sed -r 's/_[0-9]+indpMarkers//')"
echo ${trait_label_merged_GWAS}


# Variables
# Trait data
SNP_dict_path=${HOME}/ccs/Diet_project/scripts_SNP_ascertainment/SNP_dictionaries/${cluster}
gwas=${in_path}/Diet_project/merge_GWAS_diet_by_clusters/outer_1e-03/Cluster_4_Trunk_fat_free_prot/Merged_GWAS_meta_analyses_Diet_Cluster_4_Trunk_fat_free_prot_outer_pval_1e-03.tsv.gz
pval_thres=1e-05

# 1000G
vcf_path_1000G=${in_path}/Diet_project/1000G_ph3_hg19_VCFs_Diet_clusters/${cluster}/sorted
vcf_1000G=${vcf_path_1000G}/SNPs_1000G_ph3_hg19_with_GT_Diet_${cluster}_all_chr_sorted.vcf.gz
tabix_1000G=${vcf_1000G}.tbi

# Vindija
vcf_path_Vindija=${in_path}/Diet_project/Vindija_VCFs_Diet_clusters/${cluster}/sorted
vcf_Vindija=${vcf_path_Vindija}/SNPs_Diet_${cluster}_aDNA_Vindija_all_chr_sorted.vcf.gz
tabix_Vindija=${vcf_Vindija}.tbi

# Outpath
#alternative="two-sided"
out_dir=${in_path}/Diet_project/Homo_lineages_tests/Vindija/${cluster}/
log_file=${out_dir}/log_fisher_test_archaics_${cluster}.txt
mkdir -p ${out_dir}


# Initializate log
if [ ! -f ${log_file} ]
then
    echo "trait oddsratio pval_two_sided pval_greater pval_less archaic n_prototypes n_prot_original" | tr " " "\t" > ${log_file}
fi



# Launch vegeta
#--alternative ${alternative}
python ${python_dir}/Assess_selective_pressures_associated_SNPs_homo_vs_AMH_lineages.py \
--SNP_dict ${SNP_dict_path}/${SNP_dict} \
--vcf_archaic ${vcf_Vindija} \
--vcf_1000G ${vcf_1000G} \
--tabix_file_archaic ${tabix_Vindija} \
--tabix_file_1000G ${tabix_1000G} \
--out_dir ${out_dir} \
--info_outfile ${trait_clean} \
--gwas ${gwas} \
--pval_threshold ${pval_thres} \
--trait_label_merged_GWAS ${trait_label_merged_GWAS} \
--log_filename ${log_file}
