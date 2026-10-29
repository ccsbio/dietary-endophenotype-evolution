#!/bin/bash
zcat gencode.v47lift37.annotation.gff3.gz | awk -F'\t' 'BEGIN{OFS="\t"}
$3=="gene"{
  gt="."; lv="."; gid=".";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i]~/^gene_type=/){sub(/^gene_type=/,"",a[i]); gt=a[i]}
    if(a[i]~/^level=/){sub(/^level=/,"",a[i]); lv=a[i]}
    if(a[i]~/^ID=/){sub(/^ID=/,"",a[i]); gid=a[i]}
  }

  # filter: protein_coding and level 2 or 3
  if(gt=="protein_coding" && (lv=="2" || lv=="3")){
    chrom=$1;
    start=$4-1;   # GFF3 → BED conversion
    end=$5;

    name=gid"|gene_type="gt"|level="lv;

    print chrom, start, end, name
  }
}' > gencode.v47lift37.annotation_protein_coding_genes_level2-3.bed
