#!/bin/bash

REFFASTA= Solidis_ref.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp

export PATH=$GATKDIR:$PATH

VCF=$1     # this should be the name *without* .vcf extension!

echo Run started
date

# Separate snps 
gatk SelectVariants \
    --tmp-dir $TMP \
    -R $REFFASTA \
    -V $VCF.vcf \
    --select-type-to-include SNP \
    -O $VCF.snps.vcf

# Filter SNPs
gatk VariantFiltration \
    --tmp-dir $TMP \
    -R $REFFASTA \
    -V $VCF.snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter" \
    -O $VCF.snps.filtered.vcf

# Separate indels
gatk SelectVariants \
    --tmp-dir $TMP \
    -R $REFFASTA \
    -V $VCF.vcf \
    --select-type-to-include INDEL \
    -O $VCF.indels.vcf

# Filter indels
gatk VariantFiltration \
    --tmp-dir $TMP \
    -R $REFFASTA \
    -V $VCF.indels.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0" \
    --filter-name "my_indel_filter" \
    -O $VCF.indels.filtered.vcf

# Merge filtered files
gatk --java-options "-Djava.io.tmpdir=$TMP" MergeVcfs \
   -R $REFFASTA \
   -I $VCF.snps.filtered.vcf \
   -I $VCF.indels.filtered.vcf \
   -O $VCF.filtered.vcf \

echo Run ended
date
