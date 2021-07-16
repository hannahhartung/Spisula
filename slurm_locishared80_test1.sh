#!/bin/bash -l

#SBATCH --partition=regular
#SBATCH --job-name=sharedloci80_test
#SBATCH --output=test_sharedloci80_1.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hh693@cornell.edu

for M in 1 2 3 4 5 6 7 8 9
do
stacks_dir=stacks.M$M
out_dir=$stacks_dir/populations.r80
log_file=$out_dir/populations.oe
/programs/stacks/bin/populations -P $stacks_dir -O $out_dir -t 24 -r 0.80 &> $log_file
done
