# Dietary endophenotype evolution: provisional code repository

Code supporting the manuscript:

**Beyond single-trait analysis: ancient genomes reveal coordinated evolution of dietary endophenotypes in Holocene Europe**

Carla Casanova-Suárez, Karen Hardy and Oscar Lao.

Prepared for scientific review.

Preprint (bioRxiv): [10.64898/2026.10.27.762183](https://doi.org/10.64898/2026.10.27.762183).

## Release status and scope

This independent publication repository contains selected scripts from the computational pipeline executed on the CESGA high-performance computing cluster. It documents the cluster analyses, their configured parameters, and their input/output relationships. Some intermediate files and results were subsequently processed locally; those local workflows are outside the present release.

**This is a preliminary, partial code release. It has not been demonstrated to execute independently or reproduce the manuscript results from the contents supplied here.** Internal analytical functions, GWAS/genotype datasets, intermediate products and generated results are omitted. Selected job arrays, supporting annotations, sample lists and temporal-analysis metadata are included; an array listing a file does not mean that the listed file is distributed.

The descriptions below are based on the included source and author clarification. Explicit script calls and file references establish configured relationships, but do not establish that every launcher was submitted with its present array range or that all file handoffs have been validated. The manuscript remains the reference for the scientific interpretation; this README does not assign unverified figure numbers to scripts.

## Directory guide

The existing organization is retained. Most folders correspond to one analytical component; R and Python entry points are generally called by shell launchers in the same component.

```text
Allentoft_aDNA/                         Allentoft download, sample lists and QC
Download_archaic_genomes/               Vindija download
Download_human_variation_data_1000G_phase3_hg19/
                                       1000G download, population panel and EUR list
SNPs_linked_to_genes_Allentoft/          Gene-linked null-set preparation
  scripts_subset_SNPs_linked_to_genes/   GENCODE annotation and gene/SNP proximity
  scripts_subset_Allentoft/             Allentoft null-set VCF preparation
  scripts_1000G_vcf_subset/              1000G annotations for the gene-linked set
scripts_ldsc_regression/                 GWAS preparation and pairwise LDSC
  Example_arrays/                       GWAS and LDSC job specifications
  python_scripts/                       Genome-build conversion entry point
scripts_GA_clustering/                   Genetic-algorithm trait clustering with k-medoids
  R_scripts/clustering_kmedoids/
scripts_merge/                          Cluster-specific GWAS merging
  arrays/                               GWAS filenames, labels, builds and paths
  download_GWAS_diet/                    Original acquisition commands
  parsers_GWAS_cols/                     Coordinate-column preparation
  python_scripts/
scripts_subset_aDNA/subset_Allentoft/    Allentoft VCF subsetting and processing
  python_scripts/                       BED-file creation
scripts_subset_1000G/                   1000 Genomes Phase 3 VCF processing
scripts_SNP_ascertainment/              Allele matching, PLINK clumping, dictionaries
  arrays/                               Trait lists used by job arrays
  python_scripts/
scripts_proportion_risk_alleles/Allentoft/
  python_scripts/                       Proportions and temporal statistics
  metadata_samples_trend_analysis/       Metadata for 946 temporal-analysis samples
scripts_null_distribution_linked_SNPs/  Frequency-matched, gene-linked null analyses
  arrays_clusters/                      Configured trait selections
  python_scripts/
scripts_GA_geo_trend/                   Geographic differences in temporal trends
  arrays_clusters/                      Configured significant-trait selections
  R_scripts/
scripts_subset_vcf_archaics/            Vindija VCF preparation
scripts_selective_pressures_archaics/   Archaic-versus-AMH lineage comparisons
  arrays_clusters/
  python_scripts/
```

`AGENTS.md` contains repository-maintenance instructions. `.codex/` contains local tooling configuration; neither is part of the scientific pipeline.

## Workflow overview

Solid arrows summarize configured analytical dependencies. Dashed arrows mark preparation, selection or staging steps that are not fully represented in this release. The diagram is a reading guide, not an automated execution graph. File handoff qualifications are listed under [Current limitations](#current-limitations).

```mermaid
flowchart TD
    GWAS[GWAS acquisition and coordinate preparation] --> LDSC[LDSC annotation, munging and pairwise regression]
    LDSC -. MDS preparation outside this release .-> GA[GA trait clustering using k-medoids]
    GA -. Cluster post-refinement and specification .-> MERGE[Cluster-specific GWAS merging]
    GWAS --> MERGE
    MERGE --> BED[Cluster SNP BED files]
    RAWA[Allentoft download] --> QC[QC using the 967-sample EUR list]
    QC --> ADNA
    RAWK[1000 Genomes download and EUR sample list] --> KG
    RAWV[Vindija download] --> VIND
    BED --> ADNA[Allentoft VCF subsets]
    BED --> KG[1000 Genomes VCF subsets and EUR processing]
    BED --> VIND[Vindija VCF subsets]
    ADNA --> ASC[Allele matching, within-cluster trait-specific SNP selection, clumping and dictionaries]
    KG --> ASC
    MERGE --> ASC
    ASC --> TEMP[Risk-allele proportions and temporal permutation tests]
    TEMP -. Multiple-testing processing outside this release .-> SEL[Included BH-labelled trait selections]
    ASC --> AF[Index SNP risk-allele frequency distribution]
    KG --> AF
    AF --> NULL[Gene-linked null PRS matched to the frequency distribution]
    ANNO[GENCODE annotation and complex-region exclusion] --> NEAR[Gene and Allentoft SNP proximity]
    QC --> NEAR
    NEAR -. Final distance filtering outside this release .-> NSET[Gene-linked SNP BED]
    NSET --> NVCF[Allentoft null-set VCF preparation]
    QC --> NVCF
    NSET --> KANNO[1000 Genomes null-set annotation preparation]
    RAWK --> KANNO
    KANNO -. Local addition of 1000G allele frequencies outside this release .-> NAF[Gene-linked SNP frequency table]
    NVCF --> NULL
    NAF --> NULL
    SEL --> NULL
    NULL --> NTAU[Aggregation and weighted tau from null PRS]
    TEMP --> NTAU
    TEMP -. Metadata combination and input preparation .-> GEO[Geographic GA]
    SEL --> GEO
    ASC --> ARCH[Archaic-versus-AMH Fisher tests]
    KG --> ARCH
    VIND --> ARCH
```

## Analytical components and execution order

Numbering orders steps within a component. It does not define a single repository-wide submission sequence. Cluster identifiers in filenames refer to analytical trait groups, whereas CESGA refers to the computing environment.

### 1. GWAS preparation and genetic correlations

The [download commands](scripts_merge/download_GWAS_diet/) record original GWAS acquisition locations. The [coordinate parsers](scripts_merge/parsers_GWAS_cols/) prepare coordinate fields in selected summary-statistic formats. Raw summary statistics are not included, and current URL accessibility and redistribution permissions have not been verified.

The LDSC sequence is documented in [scripts_ldsc_regression](scripts_ldsc_regression/):

| Step | Included launcher or entry point | Required inputs | Declared products |
| --- | --- | --- | --- |
| Build conversion | [Build_conversion_from_hg38_to_hg19_unprocessed_GWAS_Diet.sh](scripts_ldsc_regression/Build_conversion_from_hg38_to_hg19_unprocessed_GWAS_Diet.sh) calls [Fast_build_conversion_GWAS_sumstats.py](scripts_ldsc_regression/python_scripts/Fast_build_conversion_GWAS_sumstats.py) | GWAS files, chromosome/position columns, internal methods, liftOver executable and chain files | Converted coordinates and mapping products; naming is partly delegated to internal methods |
| rsID annotation | `00.04` through `00.07` shell scripts | GWAS files, column specifications, coordinate-to-rsID table | Summary statistics with LDSC rsIDs |
| Munging | `01.*Munge_sumstats*` launchers | Annotated statistics, sample-size and allele fields, LDSC, HapMap3 allele list | `.sumstats.gz` files and munging logs |
| Pairwise analysis | [02.00.Ldsc_regression_Diet.sh](scripts_ldsc_regression/02.00.Ldsc_regression_Diet.sh) | Munged statistics and European reference LD scores | Pairwise LDSC logs and an analysis-issue log |
| Tabulation | [03.00.Table_with_pairwise_ldsc_results.sh](scripts_ldsc_regression/03.00.Table_with_pairwise_ldsc_results.sh) | LDSC result logs | `Diet_GWAS_pairwise_ldsc_regression_94_traits.tsv` |

The analysis uses hg19/GRCh37 coordinates. The author confirms that two GWAS were converted from hg38 for the LDSC analysis. The [conversion array](scripts_ldsc_regression/Example_arrays/array_meta-analyses_Diet_build_conversion_hg38_to_hg19.tsv) identifies `GCST90475156.tsv.gz` and `GCST90475153.tsv.gz`. The pairwise LDSC array contains 94 entries; this is distinct from the 58-trait input named by the clustering launchers.

[IMPORTANT_munge_sumstats_bug.txt](scripts_ldsc_regression/IMPORTANT_munge_sumstats_bug.txt) describes a manual LDSC modification for extremely small p-values. The installed, modified LDSC source is not distributed here. Its exact version and behavior therefore remain an external dependency.

### 2. Genetic-algorithm clustering of traits with k-medoids

[scripts_GA_clustering](scripts_GA_clustering/) contains K2 through K10 launchers named `02.03.Elitism.launch_GA_clustering_KMEDOIDS_traits_Diet_100iter_400pop_0_2MR_K*_4dim.sh`, calling [Run_GA_clustering_traits_Diet_by_kmedoids.R](scripts_GA_clustering/R_scripts/clustering_kmedoids/Run_GA_clustering_traits_Diet_by_kmedoids.R).

- **Input:** `LDSC_MDS_matrix_Diet_58_GWAS_interval_4dim_BMIUKB.tsv`, with a `Codes` column and MDS coordinate columns. This table and its preparation workflow are not included.
- **Visible operation:** the R entry point computes Euclidean distances from `D1` and `D2` and passes the distance matrix to `evolve.or.die.clusters.parallel.kmedoids`. The `4dim` filename does not establish that all four dimensions are used by this entry point. The function is expected from the external `Functions_GA_clustering_traits_by_kmedoids.R` file sourced under `R_scripts/clustering_kmedoids/`; that file and its implementation are not included.
- **Output:** `Best_clustering_GA_Diet_*.rda`, under the configured `GA_clustering_traits_Diet/using_kmedoids/dimensions_4/Elitism/` results directory.

The current k-medoids/interval configuration has not been executed or validated against manuscript outputs as part of repository preparation.

The four-cluster GWAS specifications supplied in the merging component reflect the post-refinement assignments confirmed by the author. In particular, `DHA_BioPsych_2023` was reassigned from cluster 4 to cluster 3 during post-refinement; the supplied merging and downstream ascertainment arrays use cluster 3. The conversion from initial clustering results to these refined specifications is not implemented by an included script.

### 3. Cluster-specific GWAS merging

The four `01.*Launch_outer_merge*` launchers in [scripts_merge](scripts_merge/) call their corresponding Python scripts in `python_scripts/`. Each uses an included specification table in `arrays/`, containing GWAS filenames, phenotype labels, genome builds and original input paths.

The merge threshold is `1e-03`. This broader threshold was used to retain more data in the initial intermediate files, allowing the downstream analysis threshold to be changed without rebuilding those files. It is not the association threshold used for the final SNP ascertainment analysis, which is `1e-05`. The custom merge call enables INDEL filtering and restricts variants to autosomal, biallelic sites. Outputs are `Merged_GWAS_meta_analyses_Diet_*_outer_pval_1e-03.tsv.gz`. Duplicate genomic positions (`chr`, `pos`) are checked as a merge QC step: when one to nine repeated rows are detected, an additional `_dup_filt.tsv.gz` table is written, excluding all rows at the repeated positions. The QC count measures repeated rows beyond the first occurrence; the filtering mask uses `keep=False` to remove every occurrence at those positions. At ten or more repeated rows, the script issues a warning for human review rather than automatically filtering them. After each merge job, inspect the log and output tables before proceeding. The cluster 1 and 2 BED-creation launchers use the duplicate-filtered tables.

The configured groups use labels based on logTG, trunk fat percentage, total cholesterol and trunk fat-free mass. The historical string `Trunk_dat_free_prot` is retained in some paths and output names; these literal names are part of the file interfaces.

### 4. Ancient and present-day VCF processing

#### Downloads and Allentoft QC

| Component | Included entry point | Inputs and products |
| --- | --- | --- |
| Allentoft download | [00.00.Download_aDNA_Denmark_imputed_chr1-22_VCF.sh](Allentoft_aDNA/00.00.Download_aDNA_Denmark_imputed_chr1-22_VCF.sh) | Downloads the 22 autosomal imputed VCFs, tabix indexes and source README from the configured ERDA archive to `DATA/ancient_DNA_Denmark/VCF/` |
| Allentoft QC | [01.00.Quality_control_filters_EURsamples_Allentoft_imputed_chr1-22_VCF.sh](Allentoft_aDNA/01.00.Quality_control_filters_EURsamples_Allentoft_imputed_chr1-22_VCF.sh) | Uses the included chromosome array and 967-sample EUR list; writes chromosome VCFs under `VCF_EURsamples_after_QC_diploid_imputed/`, as consumed by subsequent subsetting |
| 1000 Genomes download | [Download_hg37_vcf_phase3_chr1-22.sh](Download_human_variation_data_1000G_phase3_hg19/Download_hg37_vcf_phase3_chr1-22.sh) | Downloads the 22 Phase 3 autosomal VCFs from release `20130502` and the population panel to `DATA/1000G_human_variation_phase3/VCFs_1000G_phase3_hg19/` |
| Vindija download | [00.01.Download_Vindija_chr1-22_VCF.sh](Download_archaic_genomes/00.01.Download_Vindija_chr1-22_VCF.sh) | Downloads the 22 `Vindija33.19/chr*_mq25_mapab100.vcf.gz` files and source README to `DATA/Archaic_genomes/Neanderthal/Vindija/` |

The Allentoft QC command selects samples with `bcftools view -S`, applies `MAF >= 0.01 & INFO >= 0.8`, and passes the result through `bcftools +setGT` with `-t q -n . -e 'GP[*] >= 0.8'`. SNPs from both the 1000G MASK and the 1240K array are retained, subject to these QC criteria; this step does not restrict the input to one of the two sets. The script describes the data as diploid but contains no separate ploidy-selection command. The original imputation workflow is not included.

The initial European subset used for Allentoft QC contains 967 samples, including five from the Canary Islands, and is specified by [Sample_file_967_aDNA_EUR_ok.txt](Allentoft_aDNA/sample_files/Sample_file_967_aDNA_EUR_ok.txt). **Canary Islands samples are excluded from downstream temporal analyses through the 946-sample metadata table described below.** This later selection also restricts the analysis to the 400-15,000-year BP interval. The QC sample list and the downstream analysis metadata therefore serve different stages.

The included [EUR.samples.txt](Download_human_variation_data_1000G_phase3_hg19/EUR.samples.txt) contains exactly the 503 samples labelled `EUR` in the included [1000G population panel](Download_human_variation_data_1000G_phase3_hg19/integrated_call_samples_v3.20130502.ALL.panel). The author generated the list from that panel; the generation command is not included. The EUR-processing launchers expect it in the downloaded VCF directory, so it must be staged there. Download URL accessibility has not been checked as part of this release.

#### Cluster-specific VCF processing

The branches operate on autosomal chromosome files, 1 through 22:

| Component | Local sequence | Required inputs | Declared products |
| --- | --- | --- | --- |
| [Allentoft](scripts_subset_aDNA/subset_Allentoft/) | `00.*` BED creation; `01.03` VCF subsetting; `02.*` sorting/indexing; `03.*` concatenation | Merged GWAS, preprocessed Allentoft VCFs, BED specifications and chromosome list | Cluster BED files, chromosome VCF subsets, concatenated sorted VCFs and tabix indexes |
| [1000 Genomes](scripts_subset_1000G/) | `01.*` subsetting; `02.*` sorting/indexing; `03.*` concatenation; `04.*` EUR processing | Phase 3 hg19 VCFs, cluster BED files, chromosome list and `EUR.samples.txt` | Cluster VCFs and EUR subsets with recalculated tags and assigned variant IDs |
| [Vindija](scripts_subset_vcf_archaics/) | `01.*` header preparation/subsetting; `02.*` sorting/indexing; `03.*` concatenation | Vindija VCFs, cluster BED files, included contig header and chromosome list | Cluster-specific sorted/indexed Vindija VCFs |

The cluster-specific Allentoft scripts consume the QC outputs described above, whose filenames include `EUR_967samples_diploid_INFO_0.8_MAF_0.01_GP_0.8_missing`.

Allentoft subsetting is performed first. The cluster BED files, `SNPs_diet_cluster1_outer_merge_1e-03pval.bed` through `SNPs_diet_cluster4_outer_merge_1e-03pval.bed`, are generated in the Allentoft component under `scripts_subset_aDNA/subset_Allentoft/SNP_files/outer_merge_1e-03/`. These same BED files are then shared with the 1000G and Vindija subsetting steps to use a common genomic-position selection across datasets; separate SNP selections are not generated for those branches. Their configured BED locations are `scripts_subset_1000G/SNP_files/outer_merge_1e-03/` and `scripts_subset_vcf_archaics/bedfiles/outer_merge_1e-03/`, respectively.

The initial GWAS-to-Allentoft matching launchers use the concatenated VCFs under each cluster's `sorted/` directory, with the suffix `all_chr_EUR_samples.sorted.vcf.gz`, as do the risk-allele proportion analyses.

### 5. Allele matching and SNP ascertainment

The sequence in [scripts_SNP_ascertainment](scripts_SNP_ascertainment/) is:

| Prefix | Operation and analytical entry point | Products consumed downstream |
| --- | --- | --- |
| `00.*` | GWAS-to-Allentoft matching: [Merging_and_fixing_alleles_GWAS_to_Allentoft.py](scripts_SNP_ascertainment/python_scripts/Merging_and_fixing_alleles_GWAS_to_Allentoft.py) | Matched GWAS tables |
| `01.*` | 1000G matching: [Merging_alleles_GWAS_1000G_data_fix_VCF.py](scripts_SNP_ascertainment/python_scripts/Merging_alleles_GWAS_1000G_data_fix_VCF.py) | Updated GWAS tables with 1000G identifiers and filtered/sorted VCFs |
| `02.*` | PLINK binary-file preparation | `.bed`, `.bim` and `.fam` files |
| `03.*.02` | Within-cluster trait-specific SNP selection and association-file creation: [Create_assoc_file_merged_GWAS_by_trait_for_clumping_with_plink_and_1000G_data.py](scripts_SNP_ascertainment/python_scripts/Create_assoc_file_merged_GWAS_by_trait_for_clumping_with_plink_and_1000G_data.py) | `.assoc` files containing SNPs uniquely associated with the target trait within its cluster at `1e-05` |
| `04.*.02` | PLINK clumping | Parsed `.clumped` files |
| `05.*` | [SNP_ascertainment_from_PLINK_clumping_and_fixed_GWAS_with_Allentoft_aDNA_by_trait_from_merge.py](scripts_SNP_ascertainment/python_scripts/SNP_ascertainment_from_PLINK_clumping_and_fixed_GWAS_with_Allentoft_aDNA_by_trait_from_merge.py) | `.pkl` dictionaries of index SNPs and linked candidates, plus ascertainment logs |

The unique-SNP selection occurs in [Create_assoc_file_merged_GWAS_by_trait_for_clumping_with_plink_and_1000G_data.py](scripts_SNP_ascertainment/python_scripts/Create_assoc_file_merged_GWAS_by_trait_for_clumping_with_plink_and_1000G_data.py), before LD clumping. It extracts the p-value columns from the cluster-specific merged GWAS dataframe and calls the internal `number_of_unique_SNPs_associated_by_trait_and_pval_threshold` function. As clarified by the author, SNPs are selected for association with the target trait at `p < 1e-05` without association at that threshold with another trait in the same cluster dataframe. These represent **trait-specific within-cluster components**; uniqueness is assessed relative to the other traits in that cluster, not across all traits in the study. The internal selection function is not included in this release.

Uniqueness is assessed at the analysis threshold, rather than inherited from the broader `1e-03` intermediate selection. Clumping and dictionary creation also use `1e-05`. References to `1e-03` in these launchers identify upstream merged GWAS/VCF products retained for threshold flexibility; the downstream association threshold remains `1e-05`. Here, index SNPs refers to the ascertained independent markers represented by the SNP-dictionary keys, distinct from the traits representing clusters in trait clustering.

Clumping uses `--clump-p1 1e-05`, `--clump-p2 1e-05`, `--clump-r2 0.2`, `--clump-kb 500` and `--geno 0.8`. Dictionary creation supplies [complex_LD_regions_hg19_chr.tsv](scripts_SNP_ascertainment/complex_LD_regions_hg19_chr.tsv) to the internal QC function. The region table includes `chr`, `start`, `end` and `annotation`; its filtering implementation is in the omitted module.

The four trait arrays are in `arrays/`. [select_traits_by_number_of_independent_positions.sh](scripts_SNP_ascertainment/select_traits_by_number_of_independent_positions.sh) selects entries with more than 20 independent positions from an ascertainment log; [prepare_array_with_SNP_dict_filenames_using_a_list_of_selected_traits.sh](scripts_SNP_ascertainment/prepare_array_with_SNP_dict_filenames_using_a_list_of_selected_traits.sh) creates dictionary/trait arrays from a supplied selection. These helpers require explicit command-line inputs.

### 6. Risk-allele proportions and temporal trends

In [scripts_proportion_risk_alleles/Allentoft](scripts_proportion_risk_alleles/Allentoft/):

- `01.*Launch_proportion_risk_alleles*` calls [calculate_proportion_risk_alleles_Allentoft_and_weighted_or_regular_tau_using_ancient_DNA_and_SNP_dict.py](scripts_proportion_risk_alleles/Allentoft/python_scripts/calculate_proportion_risk_alleles_Allentoft_and_weighted_or_regular_tau_using_ancient_DNA_and_SNP_dict.py).
- `02.*Launch_weighted_Tau_permutation_test*` calls [calculate_observed_weighted_or_regular_tau_and_null_distribution_proportion_risk_alleles_using_ancient_DNA_and_SNP_dict.py](scripts_proportion_risk_alleles/Allentoft/python_scripts/calculate_observed_weighted_or_regular_tau_and_null_distribution_proportion_risk_alleles_using_ancient_DNA_and_SNP_dict.py).
- [prepare_df_all_traits_results_proportion_risk_alleles_clusters.sh](scripts_proportion_risk_alleles/Allentoft/prepare_df_all_traits_results_proportion_risk_alleles_clusters.sh) aggregates supplied proportion outputs.

**Inputs:** sorted/indexed Allentoft VCFs, SNP dictionaries, included trait arrays, and the included [temporal-analysis metadata](scripts_proportion_risk_alleles/Allentoft/metadata_samples_trend_analysis/Metadata_Allentoft_EUR_samples_400-15KyearsBP_946indv_diploid_imputed_noCanary.tsv). The table contains `Sample_ID`, `Age.average`, `Latitude.ok`, `Longitude.ok` and `Cluster.group`. Its 946 unique sample IDs form a subset of the 962-sample no-Canary list; none of the five Canary Islands IDs occur. Recorded ages range from 465 to 14,800 years BP, within the configured 400-15,000-year interval. Metadata membership controls which samples enter the temporal analysis, while the input VCFs retain the earlier 967-sample selection. Preparation of this metadata table is not scripted in the release.

**Products:** `Proportion_risk_allele_*.tsv`, `Presence_risk_allele_indep_*.tsv`, observed-tau logs, `Null_distribution_Kendall_tau_prop_risk_*.tsv`, permutation-test logs and aggregated proportion tables.

Although Canary Islands samples are present in the initial European VCF subset, they are excluded when selecting samples for these temporal analyses using the metadata `Sample_ID` column. The null weighted-tau comparisons use the same no-Canary temporal metadata; retaining these samples in an upstream VCF does not include them in those downstream tests.

Launchers enable weighted tau, using bin arguments start `400`, end `15000`, size `500`, and `--filt_outlier_bins_by_y 10000`. The exact weighting and outlier-bin behavior is delegated to internal methods; this argument alone should not be interpreted as a complete description of sample filtering. Permutations use 10,000 replicates and seed `3895725384`.

The supplied BH-labelled selections in subsequent components identify configured downstream traits. The multiple-testing calculation and its test family are not implemented by an included script.

### 7. Gene-linked, allele-frequency-matched null analyses

#### Gene-linked SNP set preparation

[SNPs_linked_to_genes_Allentoft](SNPs_linked_to_genes_Allentoft/) contains the upstream preparation for a set of Allentoft SNPs near protein-coding regions and their 1000G frequency annotations:

1. [download_GENCODE_release47_lift37_gene_anno.sh](SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/download_GENCODE_release47_lift37_gene_anno.sh) obtains `gencode.v47lift37.annotation.gff3.gz`. [create_bedfile_genes_from_gff3_using_filters.sh](SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/create_bedfile_genes_from_gff3_using_filters.sh) selects `gene` records with `gene_type=protein_coding` and `level` 2 or 3, converting GFF3 starts to BED coordinates.
2. [00.01.Filter_complex_LD_regions_GENCODE_gene_annotation.sh](SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/00.01.Filter_complex_LD_regions_GENCODE_gene_annotation.sh) removes gene intervals overlapping the included [complex_regions.bed](SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/complex_regions.bed), using `bedtools intersect -v`.
3. [01.01.Subset_SNPs_from_Allentoft_aDNA_linked_to_genes.sh](SNPs_linked_to_genes_Allentoft/scripts_subset_SNPs_linked_to_genes/01.01.Subset_SNPs_from_Allentoft_aDNA_linked_to_genes.sh) converts QC Allentoft variant positions to BED and runs `bedtools closest`, with genes as `-a` and SNPs as `-b`, writing per-chromosome proximity tables. This command contains no 50-kb cutoff. The final `SNPs_genome_wide_max50KB_linked_to_protein_coding_genes_GENCODEv47_hg19_levels2-3.bed` and its generation step are not included; the author indicates that distance filtering may have been performed locally in Python.
4. The [Allentoft branch](SNPs_linked_to_genes_Allentoft/scripts_subset_Allentoft/) uses that final BED to subset the QC VCFs, sort/index chromosome products and concatenate/index `SNPs_linked_to_protein_coding_genes_aDNA_Allentoft_afterQC_all_chr_EUR_samples_sorted.vcf.gz`, at the path used by the null-simulation launchers.
5. The [1000G branch](SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/) intersects Phase 3 VCFs with the same BED and retains the first nine fields for annotation processing. The final null-set dataframe, `SNPs_linked_to_protein_coding_genes_Allentoft_AF_EUR.tsv`, was prepared on the author's local computer by adding 1000G allele frequencies to the selected gene-linked Allentoft SNPs. This local preparation code and the resulting table are not included in this release.

#### Null simulations and temporal comparison

[scripts_null_distribution_linked_SNPs](scripts_null_distribution_linked_SNPs/) contains:

1. `00.*Launch_add_EUR_AF*` calls [Add_EUR_AF_from_1000G_ph3_hg19_using_SNP_dict_target_trait.py](scripts_null_distribution_linked_SNPs/python_scripts/Add_EUR_AF_from_1000G_ph3_hg19_using_SNP_dict_target_trait.py), using SNP dictionaries and 1000G annotations to produce `SNPs_*_indpMarkers_with_EUR_AF_risk_alleles.tsv`. A separate cluster-1 negative-control annotation launcher is included; a complete downstream negative-control simulation sequence is not supplied.
2. `01.*Launch_null_PRS*` calls [Null_distribution_PRS_SNPs_linked_to_genes_matching_AF.py](scripts_null_distribution_linked_SNPs/python_scripts/Null_distribution_PRS_SNPs_linked_to_genes_matching_AF.py). Inputs include the index-SNP risk-allele frequency table, an external gene-linked null SNP table with `chr`, `pos` and `ALT_AF`, an indexed Allentoft VCF for that set, and GWAS data for excluding associated SNPs. Outputs are per-trait `Null_PRS_*.tsv` tables.
3. [merge_null_PRS_cluster.sh](scripts_null_distribution_linked_SNPs/merge_null_PRS_cluster.sh) combines per-trait products into `processed_df/Null_PRS_1000replicates_allentoft_SNPs_50Kb_togenes_<cluster>.tsv`.
4. `02.*Launch_weighted_tau_distribution*` calls [From_null_PRS_to_weighted_tau_distribution.py](scripts_null_distribution_linked_SNPs/python_scripts/From_null_PRS_to_weighted_tau_distribution.py), using the combined null PRS, sample metadata and observed permutation-test tau logs. Outputs are weighted-tau distributions and comparison logs.

For each trait, matching targets the empirical distribution of European risk-allele frequencies across its ascertained index SNPs. The entry point supplies the full table of individual `Risk_allele_AF` values to the internal matching function, rather than a single mean frequency. The matching implementation is omitted; exact matching behavior cannot be verified from this release alone. The historical argument name `--prototype_SNPs_file` refers to this index-SNP frequency table.

The simulation launchers configure 1,000 replicates, seed `1234`, AF span `0.05`, and exclusion of target-trait-associated SNPs at `1e-05`. The target null set is described by filenames and interfaces as SNPs within 50 kb of protein-coding genes. Its preparation is only partially represented, as detailed above; final distance filtering and EUR frequency-table generation remain outside this release. The weighted-tau step uses the temporal bin and outlier arguments described above. The null weighted-tau launchers expect the temporal metadata under their own configured component path; staging the supplied table there is not automated.

### 8. Geographic variation in temporal trends

The three launchers in [scripts_GA_geo_trend](scripts_GA_geo_trend/) call [Run_GA_geographic_trends_split_Diet_by_trait.R](scripts_GA_geo_trend/R_scripts/Run_GA_geographic_trends_split_Diet_by_trait.R).

**Inputs:** included significant-trait arrays and external per-cluster tables named `Diet_all_traits_<cluster>_classified_Allentoft_temporal_trends_400-10yearsBP_results.tsv`. The entry point requires `Trait` and the configured columns `Age.average`, `Proportion_risk_allele`, `Longitude.ok` and `Latitude.ok`.

**Output:** `Best_geographic_trend_differences_GA_*.rda`. Preparation of the combined input tables and the custom geographic GA implementation are not included. The literal `400-10yearsBP` filename is retained without inferring its meaning from the name alone.

### 9. Archaic-versus-AMH lineage comparisons

Four launchers in [scripts_selective_pressures_archaics](scripts_selective_pressures_archaics/) call [Assess_selective_pressures_associated_SNPs_homo_vs_AMH_lineages.py](scripts_selective_pressures_archaics/python_scripts/Assess_selective_pressures_associated_SNPs_homo_vs_AMH_lineages.py).

**Inputs:** included trait arrays, SNP dictionaries, merged GWAS, indexed Vindija and 1000G VCFs. The entry point combines GWAS effects, ancestral-allele/frequency annotations and archaic calls, recovers independent positions using the dictionaries, and delegates Fisher tests to an internal method.

**Products:** annotated independent-marker tables and logs containing odds ratios and two-sided, greater and less p-values. The configured GWAS association threshold is `1e-05`. The internal test implementation is omitted, so its contingency-table construction cannot be assessed fully from this entry point.

## Computational environment and interfaces

The launchers were configured for Bash, CESGA modules and SLURM. Many use Bash process substitution, UNIX text-processing utilities, `bgzip` and `tabix`. Most jobs are intended to be submitted with `sbatch`, with their original `#SBATCH -D` working directories and log paths. They are not portable local-run instructions as supplied.

| Use | Active configuration visible in launchers |
| --- | --- |
| Python/internal methods | `cesga/system miniconda3/22.11.1-1`; conda environment `vegeta_last` |
| LDSC | The same Miniconda module; conda environment `ldsc` |
| R analyses | `cesga/2022 gcc/system`; `R/4.3.3` |
| VCF/BED operations | `bedtools/2.30.0`; `bcftools/1.19`; some steps use `cesga/2020 bcftools/1.19` |
| LD clumping | `plink/1.9b5` |

Direct Python imports include pandas, NumPy and pysam, together with the internal modules listed below. R entry points load doParallel, doRNG, dplyr, tidyr and, for clustering, stringr; they also call foreach interfaces. Python and R package versions, LDSC revision, liftOver version and complete environment definitions are not pinned in this release. Dependencies of the omitted functions cannot be enumerated completely.

In the gene-linked 1000G preparation branch, [02.01.Filtering_sample_genos_cols_per_chr_1000G_SNPs_linked_to_protein_coding_genes.sh](SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/02.01.Filtering_sample_genos_cols_per_chr_1000G_SNPs_linked_to_protein_coding_genes.sh) and [03.01.Merge_1000G_ALL_chr_filtered_VCF_SNPs_linked_to_protein_coding_genes.sh](SNPs_linked_to_genes_Allentoft/scripts_1000G_vcf_subset/03.01.Merge_1000G_ALL_chr_filtered_VCF_SNPs_linked_to_protein_coding_genes.sh) contain no active module-loading commands. They require `bgzip` and the shell/text-processing utilities to be available in the execution environment.

### Command-line arguments

Shell launchers provide the configured argument combinations and are the primary reference for a particular analysis. Python argument declarations are in each entry point's `argparse` section. Common groups are:

| Entry point/group | Main argument interface |
| --- | --- |
| Build conversion | `--file`, `--dir`, `--chr_col`, `--pos_col`, `--original_build`, `--directory_liftOver` |
| GWAS merging | `--array_file`, `--out_dir`, `--pval_threshold` |
| Allele matching | `--gwas_data_file`, `--in_dir`, `--vcf`, `--out_dir`, `--info_outfile` |
| Association files/dictionaries | `--gwas_data_file`, `--phenotype_name`, `--pval_threshold`, output/log arguments; dictionary creation additionally uses `--plink_clumped_file` and `--regions_qc` |
| Temporal statistics | VCF/index, `--SNP_dict`, metadata, y-variable, bin/filter and output/log arguments; permutations additionally require `--number_replicates`, `--cpus` and `--seed` |
| Null PRS | `--null_set_file`, `--prototype_SNPs_file`, VCF/index, replicate/CPU/seed/AF-span arguments, and configured GWAS-exclusion options |
| Null weighted tau | `--null_PRS_cluster_results`, `--tau_log_observed`, metadata, bin/filter, CPU and output/log arguments |
| Archaic comparisons | `--SNP_dict`, archaic/1000G VCFs and indexes, GWAS, trait label, threshold and output/log arguments |

The R runners use positional arguments in this order:

```text
Trait clustering:
population_size mutation_rate iterations reproduction_type K
MDS_file output_directory number_of_cores [tournament_size]

Geographic analysis:
population_size mutation_rate iterations reproduction_type
sample_table x_variable y_variable longitude_column latitude_column
output_directory trait number_of_cores [tournament_size]
```

Both R runners require an output directory ending in `/`; several Python interfaces also require trailing directory separators. Preserve the exact launcher arguments when interpreting file paths and outputs.

### GA settings and resource requests

Trait-clustering launchers configure k-medoids, population 400, mutation rate 0.2, 100 iterations, elitist reproduction, K2-K10 and 16 workers. Geographic launchers configure population 200, mutation rate 0.2, 100 iterations, tournament reproduction, tournament size 5 and 16 workers. The visible R runners do not set an explicit random seed; random-state handling inside omitted functions is not established here.

SLURM requests range from `--mem=5000` to `--mem=150000`, with walltimes from one to nine hours. Parallel launchers commonly request 20 CPUs and pass 16 or 17 workers to the analysis. These are configured allocations, not measured minimum requirements or execution benchmarks. Input data and memory use vary by stage and trait.

## Included specifications and omitted material

Included specifications comprise GWAS merging/LDSC arrays, SNP-ascertainment trait lists, temporal-analysis arrays, downstream significant-trait arrays, the Allentoft cluster-to-BED array, complex-LD-region annotations, and chromosome/contig support files. Included sample information comprises Allentoft sample lists, the 946-sample temporal metadata table, the 1000G population panel and its EUR sample list.

The following components are not supplied:

- Internal Python modules `merge_GWAS_phenos`, `work_vegeta` and `ancient_vegeta`, imported from the configured `vegeta_last_version` installation.
- `Functions_GA_clustering_traits_by_kmedoids.R`, including the called `evolve.or.die.clusters.parallel.kmedoids` implementation.
- `Functions_GA_geographic_trend.R`, including `evolve_or_die_geo_trends_across_time` and the geographic optimization details.
- GWAS summary statistics; ancient, present-day and archaic VCFs/indexes; gene-linked null SNP/frequency datasets. Sample information is limited to the lists, panel and temporal metadata identified above.
- Intermediate merged statistics, BED/PLINK files, clumping results, SNP dictionaries, MDS/combined metadata tables, simulation products, analysis logs and final results.
- Local workflows for processing intermediate files/results; the MDS preparation, final cluster selection/export, multiple-testing processing, geographic input preparation, final gene-distance filtering and null-set EUR frequency-table generation are not reconstructed in this release.
- LDSC software and its installed modification, reference LD scores, HapMap3 allele list, liftOver executable/chains and complete software environments.

These omissions must not be interpreted as instructions to recreate internal methods or substitute data. Sharing status must be assessed separately for each item.

## Current limitations

1. **External arrays:** six UKB lists referenced by LDSC are absent: `array_UKB_Anthropometric.txt`, `array_UKB_Metabolic.txt`, `array_UKB_Psychiatric.txt`, and their corresponding `_with_rsID.txt` lists. Four coordinate-parser lists are also absent: `array_parse_UKB_Ant.txt`, `array_parse_UKB_Met.txt`, `array_parse_UKB_Psyc.txt` and `array_parser_diet_composition.tsv`. Launchers point to their original CESGA data locations.
2. **Support files and staging:** `ldsc_rsID_with_coords_list.tsv` is referenced by LDSC annotation scripts but not distributed. `chromosomes.txt` is included under `scripts_subset_vcf_archaics/` and the gene-linked Allentoft branch, but is absent at the separate paths expected by the cluster-specific Allentoft and 1000G sorting launchers. The included `EUR.samples.txt` must be staged in the configured 1000G data directory; the null weighted-tau launchers also expect a separate component-local copy of the supplied temporal metadata. Final gene-linked SNP BED and EUR frequency tables remain unavailable.
3. **Original path assumptions:** absolute CESGA code, input, scratch, output and log locations remain in the scripts, along with `$HOME`, `$LUSTRE` and working-directory assumptions. Folder placement in this repository does not recreate those locations.
4. **Job arrays and orchestration:** launchers whose array files are included cover all rows, from 1 to the number of entries. Maximum indices for absent external arrays cannot be verified from this release. Array row order controls trait/job selection and is part of the interface. There is no master dependency-aware submission script.
5. **Incomplete methodological visibility:** core calculations depend on omitted functions. The visible entry points and parameters do not establish the complete algorithms, multiple-testing family, allele-handling behavior or geographic GA implementation.
6. **Validation status:** repository preparation has used static inspection and file/reference checks. No cluster analysis or comparison against manuscript outputs has been performed as part of preparing this release. Passing documentation checks would not demonstrate scientific reproduction.

## Editorial review and additional material

Requests for consideration of additional material should identify the relevant analytical component and files through the manuscript's editorial correspondence. This repository does not establish authorization to release the omitted internal modules, R functions, datasets or result files. Any provision would require confirmation by the authors of ownership, applicable sharing restrictions and specific release authorization. No availability-on-request commitment is made by this provisional release.

Please cite the [bioRxiv preprint](https://doi.org/10.64898/2026.10.27.762183) when referring to these analyses.
