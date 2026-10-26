#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet
#SBATCH --job-name=download-GWAS
#SBATCH --cpus-per-task=6
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet/err/download_GWAS_diet_last_candidates_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet/out/download_GWAS_diet_last_candidates_%j.out
#SBATCH --mem=10000
#SBATCH --time=05:00:00

main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
path_data=${main}/DATA/GWAS_data_diet/Last_candidates_diet

# Create directories
mkdir -p ${path_data}/Metabolic
mkdir -p ${path_data}/Anthropometric
mkdir -p ${path_data}/Food_liking
mkdir -p ${path_data}/Diet_habits

# Download
# Vitamin D meta-analysis
wget https://cnsgenomics.com/data/revez_20/Revezetal2020_25OHD.gz -P ${path_data}/Metabolic
wget https://cnsgenomics.com/data/revez_20/Revezetal2020_25OHD_Readme.pdf -P ${path_data}/Metabolic

# RG and FG meta-analysis
wget https://cnsgenomics.com/data/qiao_et_al_2023_nc/Qiao_et_al_2023_metaGLU.txt -P ${path_data}/Metabolic
wget https://cnsgenomics.com/data/qiao_et_al_2023_nc/Description_RG_GWAS_summary_statistics.pdf -P ${path_data}/Metabolic

# Estimated bone mineral density
wget http://www.gefos.org/sites/default/files/BEurope-Bmd-As-C-Gwas-SumStats.txt_0.gz -P ${path_data}/Antropometric

# BMI newer and larger GWAS (Veteran program US)
# BMI mean value across time points
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475156/GCST90475156.tsv.gz -P ${path_data}/Antropometric
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475156/GCST90475156.tsv.gz-meta.yaml -P ${path_data}/Antropometric
# BMI maximum value across time points
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475153/GCST90475153.tsv.gz -P ${path_data}/Antropometric
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475153/GCST90475153.tsv.gz-meta.yaml -P ${path_data}/Antropometric

# Liking food
# F-Acquired
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094791/GCST90094791_buildGRCh37.tsv -P ${path_data}/Food_liking
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094791/GCST90094791_buildGRCh37.tsv-meta.yaml -P ${path_data}/Food_liking
# F-highly palatable and savory 
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094833/GCST90094833_buildGRCh37.tsv -P ${path_data}/Food_liking
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094833/GCST90094833_buildGRCh37.tsv-meta.yaml -P ${path_data}/Food_liking
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094785/GCST90094785_buildGRCh37.tsv -P ${path_data}/Food_liking
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094785/GCST90094785_buildGRCh37.tsv-meta.yaml -P ${path_data}/Food_liking
# F-vegetarian liking
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094863/GCST90094863_buildGRCh37.tsv -P ${path_data}/Food_liking
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094863/GCST90094863_buildGRCh37.tsv-meta.yaml -P ${path_data}/Food_liking

# Diet habit
# fresh fruit consumption
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132979/GCST90132979_buildGRCh37.tsv.gz -P ${path_data}/Diet_habits
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132979/GCST90132979_buildGRCh37.tsv.gz-meta.yaml -P ${path_data}/Diet_habits
# Alcohol consumption
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132992/GCST90132992_buildGRCh37.tsv.gz -P ${path_data}/Diet_habits
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132992/GCST90132992_buildGRCh37.tsv.gz-meta.yaml -P ${path_data}/Diet_habits
# Oily fish consumption
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132993/GCST90132993_buildGRCh37.tsv.gz -P ${path_data}/Diet_habits
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132993/GCST90132993_buildGRCh37.tsv.gz-meta.yaml -P ${path_data}/Diet_habits
# Added salt
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90133001-GCST90134000/GCST90133001/GCST90133001_buildGRCh37.tsv.gz -P ${path_data}/Diet_habits
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90133001-GCST90134000/GCST90133001/GCST90133001_buildGRCh37.tsv.gz-meta.yaml -P ${path_data}/Diet_habits
# Bread consumption
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132981/GCST90132981_buildGRCh37.tsv.gz -P ${path_data}/Diet_habits
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132981/GCST90132981_buildGRCh37.tsv.gz-meta.yaml -P ${path_data}/Diet_habits