featureCounts -T 4 -p -M -t exon -g gene_name \
  -a /media/yunus/TOSHIBA1TB/reference_genomes/sorgum/Sbicolor_454_v3.1.1.gene.gtf \
  -o GSE80699_raw_counts_hisat-sorted.txt \
  *.bam \
  --verbose


