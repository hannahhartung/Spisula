#!/bin/bash

#FB=/programs/freebayes/bin/freebayes
FB=/programs/freebayes-v1.3.1/freebayes

REFFASTA = Spisula_ref.fa

echo Started at
date

$FB -f $REFFASTA \
Sol.sorted.dedup.bam \
Sim.sorted.dedup.bam \
--min-mapping-quality 30 \
> fb.vcf

echo Completed at
date
