Filtering Steps
================

## The Steps Being Used to Filter the .vcf Files.

#### Part 1: Subsetting to work with a smaller data set

Printing columns convienent for grep and awk

    grep "#" Sim.g.vcf > Sim.test.vcf
    
    grep -v "#" Sim.subset.vcf | awk '{ {print $10,$1,$2,$3,$4,$5,$6,$7,$8,$9} }' | tail

Counting the number pre-sorting

    wc -l Sim.g.vcf
    7,046,542 Sim.g.vcf
    
    wc -l Sim.testplus.vcf 
    6868878 Sim.testplus.vcf
    
    grep "0/0" Sim.testplus.vcf | wc -l
    5900698

Filtering Similis

    grep "1/1" Sim.testplus.vcf > Sim.filtered2.vcf
    awk '{if($3>100) {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1} }' Sim.filtered2.vcf | head
    awk -F "," '{if($3>100) {print} }' Sim.filtered2.vcf | awk '{{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}}' > Sim.filtered3.vcf

Attempts to Filter Sol by the same locations Alternative recommended by
Dr. Hare: bedtools from quinlab

    awk '{print $1, $2}' Sim.filtered3.vcf | head
    grep -v "#" Sim.g.vcf | awk '{{print}}' > Sol.testplus.vcf
    
    awk '{print $1, $2}' Sim.filtered3.vcf > Sim.filter.index.vcf
    
    
    
    Sub2 Test
    head Sol.testplus.vcf > Sol.sub2.vcf
    awk '{print $1, $2}' Sim.filtered3.vcf | head > Sim.sub2.vcf
    
    
    Sim.sub2.vcf | for i in $(cat); do
        echo "$i"
    done
    
    while IFS='' read -r LINE || [ -n "${LINE}" ]; do
        echo "${LINE}"
        grep "${LINE}" Sol.sub2.vcf | awk '{print}'
    done < Sim.sub2.vcf 
    
    while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    echo ${LINE}
    echo ${LINE} | cut -d " " --f 1
    echo ${LINE} | cut -d " " --f 2
    echo a
    echo b
    awk '{if($1==a && $2==b) {print} }' Sol.sub2.vcf
    done < Sim.sub2.vcf 
    
    
    
    
    
    
    while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    echo ${LINE}
    echo "${LINE}" | cut -d " " --f 1 > a.txt
    echo "${LINE}" | cut -d " " --f 2 > b.txt
    awk '{if($1== && $2==b.txt) {print} }' Sol.sub2.vcf
    done < Sim.sub2.vcf 
    
    
    
    while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    echo ${LINE}
    cut -d " " --f 1 "${LINE}" 
    cut -d " " --f 2 "${LINE}"
    awk '{if($1== cut -d " " --f 1 "${LINE}"  && $2==cut -d " " --f 2 "${LINE}") {print} }' Sol.sub2.vcf
    done < Sim.sub2.vcf
