#!/bin/bash -l

#SBATCH --partition=regular
#SBATCH --job-name=test_stacks1
#SBATCH --output=test_stacks1.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hh693@cornell.edu

for M in 1 2 3 4 5 6 7 8 9
do
popmap=../info/popmap.test_samples.tsv
reads_dir=../cleaned
out_dir=stacks.M$M
log_file=$out_dir/denovo_map.oe
/programs/stacks/bin/denovo_map.pl --samples $reads_dir -T 10 -O $popmap -o $out_dir -M $M -n $M -m 3 -b 1 -S &> $log_file
done