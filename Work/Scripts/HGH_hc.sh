#!/bin/bash

REFFASTA=Solidis_ref.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp

export PATH=$GATKDIR:$PATH

# We will run the genotyping on one chromosome only.
# Other chromosomes clould be handlen in separate runs, 
# possibly in parallel..

#REGION=chr2R
#No defined chromosome in mine

ACC=$1

# multi-threading does not work well here - they recommend using Queue to
# parallelize over regions, or just manually run a few in parallel...

echo Genotyping for $ACC started
date

#What is Xmz4g?
gatk --java-options "-Xmx4g"  HaplotypeCaller \
     --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf

echo Run ended
date

