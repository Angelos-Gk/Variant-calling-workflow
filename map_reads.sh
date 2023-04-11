#!/bin/bash

ref=data/ref/hg19.chr5_12_17.fa 

# specify threads to be used
threads=8

# generate index for the reference genome
bwa index data/ref/hg19.chr5_12_17.fa

mkdir BAMS

# perform mapping
read1=data/fastq/trimmed_reads/normal_trimmed_paired_r1.fastq.gz
read2=data/fastq/trimmed_reads/normal_trimmed_paired_r2.fastq.gz
RGID="231335"
RGSN="Normal"
bwa mem -t $threads -R "@RG\\tID:${RGID}\\tPL:Illumina\\tPU:None\\tLB:None\\tSM:${RGSN}" $ref $read1 $read2 | samtools view -h -b -o BAMS/normal.raw.bam


read1=data/fastq/trimmed_reads/tumor_trimmed_paired_r1.fastq.gz
read2=data/fastq/trimmed_reads/tumor_trimmed_paired_r2.fastq.gz
RGID="231336"
RGSN="Tumor"
bwa mem -t $threads -R "@RG\\tID:${RGID}\\tPL:Illumina\\tPU:None\\tLB:None\\tSM:${RGSN}" $ref $read1 $read2 | samtools view -h -b -o BAMS/tumor.raw.bam

