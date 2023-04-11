#!/bin/bash
# dont forget to set chmod 755 for the script (and for the scripts that are called) in order to be able to run it

SECONDS=0

# ==== OPTIONAL SETUP ======
# set up conda environment with the necessary packages 
conda env create -f environment.yml

# activate the environment
conda activate variants

# ==== START OF WORKFLOW ======
# navigate to working directory
cd /mnt/c/Users/angel/Desktop/Bash_workflows/variant_calling

# download the data
./download_data.sh

# run fastqc
fastqc data/fastq/raw_reads/*.fastq.gz # you can specify number of threads using -t [number of threads]

# run multiqc
multiqc data/fastq/raw_reads/*_fastqc*

# run trimmomatic
./trim_reads/sh

# run fastqc again 
fastqc data/fastq/trimmed_reads/*.fastq.gz

# run multiqc again
multiqc data/fastq/trimmed_reads/*_fastqc*

# map the reads
./maps_reads.sh

# process the bams further
./process_bams.sh

# call variants
./call_variants.sh


# ====BASIC STATISTICS====
cd variants

# list the sample groups
bcftools query -l variants.all.vcf.gz

# filter the variants
# first check the column names
zgrep '^#' variants.all.vcf.gz | tail -n1

# perform filtering keeping only the variants that have "PASS" on the FILTER column
bcftools filter -i 'FILTER="PASS"' variants.all.vcf.gz -Oz -o variants.all.pass.vcf.gz

# indexing
bcftools index variants.all.pass.vcf.gz

# check the variants that have been identified in the TUMOR sample
bcftools view -s TUMOR variants.all.pass.vcf.gz -c1 -Oz -o tumor.vcf.gz
bcftools index tumor.vcf.gz

# identify number of variants in the tumor vcf file
zgrep -v -c '^#' tumor.vcf.gz

# identify number of snps
bcftools view -v snps tumor.vcf.gz | grep -v -c '^#'

# identify number of indels
bcftools view -v indels tumor.vcf.gz | grep -v -c '^#'

# check where the variants are (which chromosome) -- use uniq so each chromosome appears once
bcftools query -f '%CHROM\n' tumor.vcf.gz | uniq

# check number of variants in chromosome 5
bcftools view --regions chr5 tumor.vcf.gz | grep -v -c '^#'

# check variants type
zgrep '^##INFO' tumor.vcf.gz #check the types of variants
bcftools query -f '%INFO/SS\n' tumor.vcf.gz | sort | uniq

# check the number of each type of variant
bcftools query -f '%INFO/SS\n' tumor.vcf.gz | sort | uniq -c


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"