#!/bin/bash

# A subset of the commands used in the pipeline described in Therkildsen and Palumbi: "Practical low-coverage genome-wide sequencing of hundreds of individually barcoded samples for population and evolutionary genomics in non-model species" published in Molecular Ecology Resources.
# usage: nohup ./Trimmomatic_mph.sh SAMPLELIST OUTPUTDIR >& NYC_trimmomatic.log &

SAMPLELIST=$1 # Filename and path to list of fastq file base names (remove _R1.fastq.gz) e.g. /workdir/mphwork/NYC_WGS_samples.txt
OUTPUTDIR=$2 # path to directory where output files are to be written, e.g. /workdir/mphwork/trimmed
##ADAPTERS=$3 # /workdir/Cod/ReferenceSeqs/NexteraPE_NT.fa is a list of adapter/index sequences for Nextera DNA and Nextera XT PE libraries

##### RUN EACH SAMPLE THROUGH PIPELINE #######
# Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST | cut -f1`; do
SAMPLE=`grep "$SAMPLEFILE" $SAMPLELIST`
SAMPLE=$OUTPUTDIR$SAMPLE'_1'

#### CLEANING THE READS ####
# Remove adapter sequence with Trimmomatic, minlength=80. 
java -jar /programs/trimmomatic/trimmomatic-0.36.jar PE -threads 22 -phred33 $SAMPLEFILE'_R1.fastq.gz' $SAMPLEFILE'_R2.fastq.gz' $SAMPLEFILE'_AdapterClipped_F_paired.fastq' $SAMPLEFILE'_AdapterClipped_F_unpaired.fastq' $SAMPLEFILE'_AdapterClipped_R_paired.fastq' $SAMPLEFILE'_AdapterClipped_R_unpaired.fastq' ILLUMINACLIP:/programs/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 MINLEN:80
done