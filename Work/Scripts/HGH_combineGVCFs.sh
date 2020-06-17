#!/bin/bash

REFFASTA=./Solidis_ref.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp

export PATH=$GATKDIR:$PATH

# the gVCFG files obyained before need to be combined to be used with GenotypeGVCFs tool

# REGION=chr2R


echo Combining GVCFs started 
date

gatk CombineGVCFs --tmp-dir $TMP -R $REFFASTA \
		--variant Sim4.g.vcf --variant Sim5.g.vcf --variant Sim6.g.vcf\
		--variant Sol6.g.vcf --variant Sol7.g.vcf --variant Sol8.g.vcf --variant Sol10.g.vcf\
		-O all.g.vcf &

echo Run ended
date

