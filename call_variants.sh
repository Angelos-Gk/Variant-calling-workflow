#!/bin/bash

mkdir variants

# generate a pileup file
ls BAMS/*.final.bam >bamlist.txt
samtools mpileup -Q 28 -q 1 -f data/ref/hg19.chr5_12_17.fa -b bamlist.txt -o variants/normal-tumor.mpileup

# call variants
varscan somatic variants/normal-tumor.mpileup variants/variants.vcf --mpileup 1 --normal-purity 1 --tumor-purity 0.5 --output-vcf 1 
#--tumor-purity estimated “Estimated purity (tumor content) of normal sample”
#--normal-purity  “Estimated purity (non-tumor content) of normal sample”

# rename files
mv variants/variants.vcf.snp variants/variants.snp.vcf
mv variants/variants.vcf.indel variants/variants.indel.vcf

# compress vcfs
bgzip variants/variants.snp.vcf 
bgzip variants/variants.indel.vcf

# index vcfs
bcftools index variants/variants.snp.vcf.gz
bcftools index variants/variants.indel.vcf.gz


# merge snp and indels in the same vcf
echo "merging snps and indels"
bcftools concat --allow-overlaps variants/variants.snp.vcf.gz variants/variants.indel.vcf.gz -Oz -o variants/variants.all.vcf.gz
bcftools index variants/variants.all.vcf.gz

