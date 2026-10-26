#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet
#SBATCH --job-name=download-UKB
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet/err/download_GWAS_UKB_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet/out/download_GWAS_UKB_%j.out
#SBATCH --mem=10000
#SBATCH --time=05:00:00


# Subdirectories
main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
path_data=${main}/DATA/GWAS_data_diet/UKB
ukb_PSY_pheno=${path_data}/UKB_Psychiatric_pheno
ukb_anthropometric_pheno=${path_data}/UKB_Anthropometric_pheno
ukb_metabolic_pheno=${path_data}/UKB_Metabolic_pheno

# Create subdirectories
mkdir -p ${ukb_PSY_pheno}
mkdir -p ${ukb_anthropometric_pheno}
mkdir -p ${ukb_metabolic_pheno}



# Antropometric
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23111_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23111_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23112_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23112_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23113_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23113_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23114_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23114_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23115_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23115_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23116_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23116_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23117_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23117_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23118_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23118_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23119_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23119_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23120_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23120_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23121_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23121_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23122_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23122_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23123_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23123_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23124_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23124_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23125_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23125_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23126_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23126_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23127_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23127_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23128_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23128_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23129_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23129_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23130_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23130_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23099_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23099_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23101_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23101_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23102_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23102_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23104_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_anthropometric_pheno}/23104_irnt.gwas.imputed_v3.both_sexes.tsv.bgz



# Metabolic
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23105_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/23105_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100002_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100002_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100003_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100003_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100004_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100004_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100005_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100005_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100006_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100006_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100007_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100007_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100008_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100008_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100009_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100009_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100011_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100011_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100012_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100012_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100013_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100013_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100014_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100014_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100015_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100015_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100016_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100016_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100017_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100017_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100018_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100018_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100019_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100019_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100021_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100021_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100023_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100023_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100024_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100024_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/100025_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_metabolic_pheno}/100025_irnt.gwas.imputed_v3.both_sexes.tsv.bgz



# Psychiatric
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20536_1.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_PSY_pheno}/20536_1.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20536_2.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_PSY_pheno}/20536_2.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20536_3.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_PSY_pheno}/20536_3.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_16.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_PSY_pheno}/20544_16.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_13.gwas.imputed_v3.both_sexes.v2.tsv.bgz -O ${ukb_PSY_pheno}/20544_13.gwas.imputed_v3.both_sexes.v2.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_12.gwas.imputed_v3.both_sexes.tsv.bgz -O ${ukb_PSY_pheno}/20544_12.gwas.imputed_v3.both_sexes.tsv.bgz