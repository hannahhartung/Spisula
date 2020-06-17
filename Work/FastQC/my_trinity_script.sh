#!/bin/bash

TRINITY_HOME=/programs/trinityrnaseq-v2.8.6
TRINITY_OUT=./trinity_out

# Make sure we can find Jellyfish (for K-mer counting)
export PATH=/programs/jellyfish-2.2.3/bin:$PATH
# Make sure we can find salmon (for expression filtering)
export PATH=/programs/salmon-0.11.3/bin:$PATH
# Make sure we use the latest Bowtie2
 export PATH=/programs/bowtie2-2.3.4.3:$PATH

$TRINITY_HOME/Trinity --seqType fq \
--left Sp_ds.left.fq.gz,Sp_hs.left.fq.gz,Sp_log.left.fq.gz,Sp_plat.left.fq.gz \
--right Sp_ds.right.fq.gz,Sp_hs.right.fq.gz,Sp_log.right.fq.gz,Sp_plat.right.fq.gz \
--SS_lib_type RF \
--max_memory 2G  \
--trimmomatic \
--CPU 2 \
--output $TRINITY_OUT


# To change default Trimmomatic options, add the folllowing to the Trinity command line, 
# after changing parameters to your liking. (Do NOT include the "#", but remeber "\" at the end of the line )
#
# --quality_trimming_params "ILLUMINACLIP:$TRINITY_HOME/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"

# To skip read normalization, add the following option 
#
# --no_normalize_reads
#
# By fdefault, Trinity normalizes reads using option
# --normalize_max_read_cov 200
