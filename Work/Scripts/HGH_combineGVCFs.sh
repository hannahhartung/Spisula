#!/bin/bash

REFFASTA=./Sol_ref.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp

export PATH=$GATKDIR:$PATH

# the gVCFG files obyained before need to be combined to be used with GenotypeGVCFs tool

# REGION=chr2R


echo Combining GVCFs started 
date

gatk CombineGVCFs --tmp-dir $TMP -R $REFFASTA \
--variant Sim4.g.vcf.gz \
--variant Sim5.g.vcf.gz \
--variant Sim6.g.vcf.gz \
--variant Sol6.g.vcf.gz \
--variant Sol7.g.vcf.gz \
--variant Sol8.g.vcf.gz \
--variant Sol10.g.vcf.gz \
-O all.g.vcf.gz &

echo Run ended
date

