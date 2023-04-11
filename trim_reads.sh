#!/bin/bash


threads=8

cdir=$(pwd)

fastqdir=${cdir}/data/fastq/raw_reads
trimdir=${cdir}/data/fastq/trimmed_reads
mkdir $trimdir
cd $trimdir

for sample in normal tumor
do 
trimmomatic PE -threads $threads ${fastqdir}/${sample}_r1.fastq.gz ${fastqdir}/${sample}_r2.fastq.gz \
${sample}_trimmed_paired_r1.fastq.gz ${sample}_trimmed_unpaired_r1.fastq.gz \
${sample}_trimmed_paired_r2.fastq.gz ${sample}_trimmed_unpaired_r2.fastq.gz \
ILLUMINACLIP:$cdir/data/TruSeq3-PE.fa:2:30:10:8:true HEADCROP:3 TRAILING:10 MINLEN:25
done

cd $cdir

