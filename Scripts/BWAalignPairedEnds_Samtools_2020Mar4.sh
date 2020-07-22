#bwa sampe to create .sam
#paired ends

bwa sampe Solidis_Trinity.fasta SimR1.sai SimR2.sai Sim_reads.fq Sim_mates.fq > Sim_12.sam
bwa sampe Solidis_Trinity.fasta SolR1.sai SolR2.sai Sol_reads.fq Sol_mates.fq > Sol_12.sam

#assuming that worked, convert to .bam
#unknown whether it will have a header or not so over night try both

#assume no header

samtools view -bT Solidis_Trinity.fasta Sim_12.sam > Sim_nohead.bam 
samtools view -bT Solidis_Trinity.fasta Sol_12.sam > Sol_nohead.bam
# when no header
#view these (and original .sam next time to check which was right)

samtools view -bS Sim12.sam > Sim_whead.bam
samtools view -bS Sol12.sam > Sol_whead.bam
 # when SAM header present
 
 
# sort by coordinate to streamline data processing 
samtools sort -O bam -o Sim_whead.sorted.bam -T temp Sim_whead.bam 
samtools sort -O bam -o Sol_whead.sorted.bam -T temp Sol_whead.bam
samtools sort -O bam -o Sim_nohead.sorted.bam -T temp Sim_nohead.bam  
samtools sort -O bam -o Sol_nohead.sorted.bam -T temp Sol_nohead.bam 


# a position-sorted BAM file can also be indexed
samtools index Sim_whead.sorted.bam 
samtools index Sol_whead.sorted.bam 
samtools index Sim_nohead.sorted.bam 
samtools index Sol_nohead.sorted.bam 