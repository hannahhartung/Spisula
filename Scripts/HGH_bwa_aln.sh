#!/bin/bash

echo Alignment started
date

ACC=Sol7
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo Alignment started
date

ACC=Sol8
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo Alignment started
date

ACC=Sol10
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo Alignment started
date

ACC=Sim4
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo Alignment started
date

ACC=Sim5
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo Alignment started
date

ACC=Sim6
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo Alignment started
date

ACC=Sol6
SAM=Solidissima
REFFASTA=Sol_ref.fa

bwa mem -M -t 22 \
 -R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
 ${ACC}_clean_R1.fastq.gz  ${ACC}_clean_R2.fastq.gz \
 | samtools view -Sb - > $ACC.bam
 

echo "$ACC" Alignment finished
date

echo All Alignments finished
date


#REFFASTA=Sol_cdhit_ref.fa

# Un-comment one of the accession/sample pair below to run alignment for this sample

#ACC=$1
#SAM=$2

#ACC=Sol6
#SAM=Solidissima

#ACC=SRR1663609
#SAM=ZW177

#ACC=SRR1663610
#SAM=ZW184

#ACC=SRR1663611
#SAM=ZW185

# bwa will run on 20 CPUs (-t 20)

#echo Alignment started
#date

#bwa mem -M -t 7 \
#-R "@RG\tID:${ACC}\tSM:${SAM}\tLB:${SAM}\tPL:ILLUMINA" $REFFASTA \
#${ACC}_R1.fastq.gz  ${ACC}_R2.fastq.gz \
#| samtools view -Sb - > $ACC.bam

#echo Alignment finished
#date
