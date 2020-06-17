#!/bin/bash

# index reference genome for bwa, create fasta indexes (fai and dict)

TMP=/workdir/$USER/tmp
GATKDIR=/programs/gatk-4.1.4.0
export PATH=$GATKDIR:$PATH

# cd genome

# Genome summary files needed and by GATK tools
gatk CreateSequenceDictionary -R Sol_cdhit_ref.fa  -O Sol_cdhit_ref.dict
samtools faidx Sol_cdhit_ref.fa

# index for BWA alignment
bwa index Sol_cdhit_ref.fa

# index image file needed by some Spark-based tools (if used)
#gatk --java-options "-Djava.io.tmpdir=$TMP" BwaMemIndexImageCreator \
#     -I Sol_cdhit_ref.fa \
#     -O Sol_cdhit_ref.fa.img
