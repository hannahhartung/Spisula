#!/bin/bash

REFFASTA=./Solidis_ref.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp

export PATH=$GATKDIR:$PATH

# run this script as follows
#
#    nohup ./sort_dedup_index.sh SRR1663609 >& log &
#

ACC=$1

# Note: GATK will create its temporary files in $TMP which is on large local disk.
# This is safer than putting them in default /tmp, which is usually small

echo Dedup/sorting started
date
gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

# Separate indexing not needed if CREATE_INDEX true in MarkDuplicates

#echo Indexing started
#date
#
#gatk --java-options "-Djava.io.tmpdir=$TMP" BuildBamIndex \
#    --INPUT=$ACC.sorted.dedup.bam

echo Run completed
date

