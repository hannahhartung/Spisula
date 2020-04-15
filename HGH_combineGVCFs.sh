#!/bin/bash

REFFASTA=./Solidis_ref.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp

export PATH=$GATKDIR:$PATH

# the gVCFG files obyained before need to be combined to be used with GenotypeGVCFs tool

# REGION=chr2R


echo Combining GVCFs started 
date

gatk CombineGVCFs \
     --tmp-dir $TMP \
     -R $REFFASTA \
     # -L $REGION \
     --variant Sim.g.vcf \
     --variant Sol.g.vcf \
     -O all.g.vcf

echo Run ended
date

