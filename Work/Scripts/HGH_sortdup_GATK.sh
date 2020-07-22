#!/bin/bash

#INITIALIZE JAVA AND GATK
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH
REFFASTA=./Sol_trim_ref_NoPf.fa
GATKDIR=/programs/gatk-4.1.4.0
TMP=/workdir/$USER/tmp
export PATH=$GATKDIR:$PATH

#BEGIN INDIVIDUALS

#SOL 6
echo Sorting Started
date

ACC=Sol6
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8' 
            
            
            >& ${ACC}_sortdup.log &

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date

#SOL 7
echo Sorting Started
date

ACC=Sol7
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date

#SOL 8
echo Sorting Started
date

ACC=Sol8
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date


#SOL 10
echo Sorting Started
date

ACC=Sol10
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date


#SIM 4
echo Sorting Started
date

ACC=Sim4
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date

#SIM 5
echo Sorting Started
date

ACC=Sim5
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date


#SIM 6
echo Sorting Started
date

ACC=Sim6
SAM=Solidissima
REFFASTA=Sol_ref.fa

gatk MarkDuplicatesSpark \
            -I ${ACC}.bam \
            -O ${ACC}.sorted.dedup.bam \
            -M ${ACC}.sorted.dedup.txt \
            --tmp-dir $TMP \
            --conf 'spark.executor.cores=8'

echo "$ACC" Sorting Finished
date

echo "$ACC" Haplotype Calling Started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
    --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC.sorted.dedup.bam  \
     -ERC GVCF \
     --native-pair-hmm-threads 10 \
     --minimum-mapping-quality 30 \
     -O $ACC.g.vcf
     
echo "$ACC" Haplotype Calling Finished
date



#ORIGINAL SCRIPT
#REFFASTA=./Solidis_ref.fa
#GATKDIR=/programs/gatk-4.1.4.0
#TMP=/workdir/$USER/tmp

#export PATH=$GATKDIR:$PATH

# run this script as follows
#
#    nohup ./sort_dedup_index.sh SRR1663609 >& log &
#

#ACC=$1

# Note: GATK will create its temporary files in $TMP which is on large local disk.
# This is safer than putting them in default /tmp, which is usually small

#echo Dedup/sorting started
#date
#gatk MarkDuplicatesSpark \
#            -I ${ACC}.bam \
#            -O ${ACC}.sorted.dedup.bam \
#            -M ${ACC}.sorted.dedup.txt \
#            --tmp-dir $TMP \
#            --conf 'spark.executor.cores=8'

# Separate indexing not needed if CREATE_INDEX true in MarkDuplicates

#echo Indexing started
#date
#
#gatk --java-options "-Djava.io.tmpdir=$TMP" BuildBamIndex \
#    --INPUT=$ACC.sorted.dedup.bam

#echo Run completed
#date

