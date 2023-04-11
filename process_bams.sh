#!/bin/bash

ref=data/ref/hg19.chr5_12_17.fa
threads=8

for sample in normal tumor 
do

echo "filtering bam file"
bamtools filter -in BAMS/${sample}.raw.bam -mapQuality ">=1" -isMapped true -isMateMapped true -out BAMS/${sample}.filtered.bam

echo "sorting by name"
samtools sort -@ $threads -n BAMS/${sample}.filtered.bam -o BAMS/${sample}.sorted.n.bam 

echo "applying fixmate"
samtools fixmate -m BAMS/${sample}.sorted.n.bam BAMS/${sample}.fixmate.bam

echo "sorting by coordinate"
samtools sort -@ $threads BAMS/${sample}.fixmate.bam -o BAMS/${sample}.sorted.p.bam 

echo "marking and removing duplicates"
samtools markdup -r -@ $threads BAMS/${sample}.sorted.p.bam BAMS/${sample}.dedup.bam 

echo "duplicates have been removed"

# left align reads

samtools faidx data/ref/hg19.chr5_12_17.fa

echo "left align reads"

cat BAMS/${sample}.dedup.bam | bamleftalign -f $ref  > BAMS/${sample}.leftalign.bam

# Recalibrate read mapping qualities

echo "recalibrating read mapping qualities"
samtools calmd -@ 8 -bAr BAMS/${sample}.leftalign.bam $ref  > BAMS/${sample}.calmd.bam


# final filtering
echo "filtering bam file"

bamtools filter -in BAMS/${sample}.calmd.bam -mapQuality "<=254" -out BAMS/${sample}.final.bam

bwa index BAMS/${sample}.final.bam

echo "bam ready file is ${sample}.final.bam"
done

