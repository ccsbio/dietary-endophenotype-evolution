#!/bin/bash
#SBATCH --nodes=1
#SBATCH -D /home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet
#SBATCH --job-name=download-GWAS
#SBATCH --cpus-per-task=6
#SBATCH --error=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet/err/download_GWAS_diet_meta_%j.err
#SBATCH --output=/home/csic/gpo/olg/ccs/Diet_project/scripts_merge/download_GWAS_diet/out/download_GWAS_diet_meta_%j.out
#SBATCH --mem=20000
#SBATCH --time=05:00:00

main=/mnt/lustre/scratch/nlsas//home/csic/gpo/olg/ccs
path_data=${main}/DATA/GWAS_data_diet


# Create directories
mkdir -p ${path_data}/Diet_composition  
mkdir -p ${path_data}/Glucemic_traits  
mkdir -p ${path_data}/Anthropometric
mkdir -p ${path_data}/Lipids  
mkdir -p ${path_data}/Liver  
mkdir -p ${path_data}/Metabolites  
mkdir -p ${path_data}/Physical_activity  

# Download
wget https://zenodo.org/records/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1 -O ${path_data}/Anthropometric/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz
wget https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz -P ${path_data}/Lipids/
wget https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz -P ${path_data}/Lipids/
wget https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz -P ${path_data}/Lipids/
wget https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/nonHDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz -P ${path_data}/Lipids/
wget https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz -P ${path_data}/Lipids/
wget https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/README.txt -P ${path_data}/Lipids/
#wget https://cnsgenomics.com/data/qiao_et_al_2023_nc/Qiao_et_al_2023_metaGLU.txt
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013405/GCST90013405_buildGRCh37.tsv -P ${path_data}/Liver/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013405/GCST90013405_buildGRCh37.tsv-meta.yaml -P ${path_data}/Liver/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013406/GCST90013406_buildGRCh37.tsv -P ${path_data}/Liver/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013406/GCST90013406_buildGRCh37.tsv-meta.yaml -P ${path_data}/Liver/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013407/GCST90013407_buildGRCh37.tsv -P ${path_data}/Liver/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013407/GCST90013407_buildGRCh37.tsv-meta.yaml -P ${path_data}/Liver/
wget https://magicinvestigators.org/downloads/MAGIC1000G_2hGlu_EUR.tsv.gz -P ${path_data}/Glucemic_traits/ --no-check-certificate
wget https://magicinvestigators.org/downloads/MAGIC1000G_FG_EUR.tsv.gz -P ${path_data}/Glucemic_traits/ --no-check-certificate
wget https://magicinvestigators.org/downloads/MAGIC1000G_HbA1c_EUR.tsv.gz -P ${path_data}/Glucemic_traits/ --no-check-certificate
wget https://magicinvestigators.org/downloads/MAGIC1000G_FI_EUR.tsv.gz -P ${path_data}/Glucemic_traits/ --no-check-certificate
wget https://magicinvestigators.org/downloads/Trans_readme.pdf -P ${path_data}/Glucemic_traits/ --no-check-certificate
wget https://csg.sph.umich.edu/willer/public/TSH2020/2020_Zhou_et_al_TSH_meta_HUNT-MGI-ThyroidOmics.txt -P ${path_data}/Metabolites/
#wget https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
#wget https://figshare.com/ndownloader/files/28169271 -O pgcAN2.2019-07.vcf.tsv.gz
wget https://personal.broadinstitute.org/ryank/Merino_2023_dietary_intake_PMID34426670.zip -P ${path_data}/Diet_composition/
wget https://kp4cd.org/sites/default/files/READMEs/README_Merino_PMID34426670.txt -P ${path_data}/Diet_composition/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104341/GCST90104341_buildGRCh37.tsv -P ${path_data}/Physical_activity/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104341/README.txt -P ${path_data}/Physical_activity/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104345/GCST90104345_buildGRCh37.tsv -P ${path_data}/Physical_activity/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104339/GCST90104339_buildGRCh37.tsv -P ${path_data}/Physical_activity/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269542/GCST90269542.tsv.gz -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269542/GCST90269542.tsv.gz-meta.yaml -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269541/GCST90269541.tsv.gz -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269541/GCST90269541.tsv.gz-meta.yaml -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269547/GCST90269547.tsv.gz -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269547/GCST90269547.tsv.gz-meta.yaml -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269572/GCST90269572.tsv.gz -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269572/GCST90269572.tsv.gz-meta.yaml -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269573/GCST90269573.tsv.gz -P ${path_data}/Metabolites/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269001-GCST90270000/GCST90269573/GCST90269573.tsv.gz-meta.yaml -P ${path_data}/Metabolites/
#wget https://figshare.com/ndownloader/files/30731276 -O pgc-panic2019.vcf.tsv.gz
#wget https://figshare.com/ndownloader/files/31424374 -O jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz
