Marker Development: Allele Frequency Calculations
================

## For Primer Pairs L116, L148, L203, L430, and L447

From my first two rounds of marker creation using alignments to the
Solidissima Transcriptome.

#### Set Up

Including notes to myself.

``` r
library(tidyverse)
library(knitr)

#Import
samples<-read_csv("allelefreq_spreadsheet.csv")
#kable(head(samples))

#Remove collumns to the right filled with NAs left over from excel
samples <- samples[colSums(!is.na(samples)) > 0]
#Filter rows at bottom with same thing (using tidyverse)
samples <- samples %>% filter_all(any_vars(!is.na(.)))
```

``` r
#going to use mutate
#Similis is 0, Solidissima is 1
#if L116 = hybrid, L116_allele1 = 0, L116_allele2 = 1
#samples %>% mutate(allele = 
#                     ifelse(L116 == "hybrid", 0, 
#                            ifelse(L116 == "similis", 0, 
#                                  ifelse(L116 == "solidissima", 1, NA))))
#cannot easily/effectively turn this into a loop because it is difficult to iterate over collumn names and even more so difficult to create new collumn name within mutate.

#Do over normal for loop
i <- c("L116", "L148","L203","L430","L447")
for (j in i){
  a1 <- paste(j,'allele','1', sep='_')
  samples[[a1]] <- ifelse(samples[[j]] == "hybrid", 0, 
                   ifelse(samples[[j]] == "similis", 0, 
                   ifelse(samples[[j]] == "solidissima", 1, NA)))
  a2 <- paste(j,'allele','2', sep='_')
  samples[[a2]] <- ifelse(samples[[j]] == "hybrid", 1, 
                   ifelse(samples[[j]] == "similis", 0, 
                   ifelse(samples[[j]] == "solidissima", 1, NA)))
}

#Remove full row NAs (from excel document - not all were used for tests)
samples <- samples[which(rowMeans(is.na(samples)) < 0.70), ] #row with 70% or more NA. 15/21 of NA variables (the number of loci*3), is 71%

kable(head(samples))
```

| Sample ID    | Site    | Location             | Site Expected Species | Year | Mitochondrial Marker | L116    | L148        | L203        | L430        | L447        | L116\_allele\_1 | L116\_allele\_2 | L148\_allele\_1 | L148\_allele\_2 | L203\_allele\_1 | L203\_allele\_2 | L430\_allele\_1 | L430\_allele\_2 | L447\_allele\_1 | L447\_allele\_2 |
| :----------- | :------ | :------------------- | :-------------------- | ---: | :------------------- | :------ | :---------- | :---------- | :---------- | :---------- | --------------: | --------------: | --------------: | --------------: | --------------: | --------------: | --------------: | --------------: | --------------: | --------------: |
| A - ARC 006  | A - ARC | Reitsman Hatchery    | Solidissima           | 2020 | solidissima          | NA      | NA          | solidissima | solidissima | solidissima |              NA |              NA |              NA |              NA |               1 |               1 |               1 |               1 |               1 |               1 |
| A - ARC 007  | A - ARC | Reitsman Hatchery    | Solidissima           | 2020 | solidissima          | NA      | NA          | solidissima | solidissima | solidissima |              NA |              NA |              NA |              NA |               1 |               1 |               1 |               1 |               1 |               1 |
| B - ICO 008  | B - ICO | Reitsman Hatchery    | Solidissima           | 2020 | solidissima          | hybrid  | solidissima | NA          | NA          | NA          |               0 |               1 |               1 |               1 |              NA |              NA |              NA |              NA |              NA |              NA |
| B - ICO 015  | B - ICO | Reitsman Hatchery    | Solidissima           | 2020 | solidissima          | NA      | NA          | solidissima | solidissima | solidissima |              NA |              NA |              NA |              NA |               1 |               1 |               1 |               1 |               1 |               1 |
| BLP0819\_181 | BLP     | Southern Long Island | Solidissima           | 2019 | solidissima          | hybrid  | solidissima | similis     | hybrid      | solidissima |               0 |               1 |               1 |               1 |               0 |               0 |               0 |               1 |               1 |               1 |
| BLP0819\_182 | BLP     | Southern Long Island | Solidissima           | 2019 | solidissima          | similis | solidissima | solidissima | solidissima | solidissima |               0 |               0 |               1 |               1 |               1 |               1 |               1 |               1 |               1 |               1 |

#### Calculations

##### Comparing Summation and Mean Methods

``` r
usingmean <- samples %>%
  group_by(`Mitochondrial Marker`) %>%
  summarize(L116 = (mean(L116_allele_1[!is.na(L116_allele_1)]) + mean(L116_allele_2[!is.na(L116_allele_2)]))/2, 
            L148 = (mean(L148_allele_1[!is.na(L148_allele_1)]) + mean(L148_allele_2[!is.na(L148_allele_2)]))/2, 
            L203 = (mean(L203_allele_1[!is.na(L203_allele_1)]) + mean(L203_allele_2[!is.na(L203_allele_2)]))/2,
            L430 = (mean(L430_allele_1[!is.na(L430_allele_1)]) + mean(L430_allele_2[!is.na(L430_allele_2)]))/2,
            L447 = (mean(L447_allele_1[!is.na(L447_allele_1)]) + mean(L447_allele_2[!is.na(L447_allele_2)]))/2)

#     of L116 allele 1,     in rows that aren't NA,        how many are exactly 1 (solidissima allele)
sum(samples$L116_allele_1[!is.na(samples$L116_allele_1)] == 1)
```

    ## [1] 28

``` r
#     the same for allele 2
sum(samples$L116_allele_2[!is.na(samples$L116_allele_2)] == 1) 
```

    ## [1] 53

``` r
# this is more because hybrids have a 0 for allele 1 and a 1 for allele 2

#the total number of non-na's from the two alleles
sum(!is.na(samples$L116_allele_1)) + sum(!is.na(samples$L116_allele_2))     #119 from each is 238 total
```

    ## [1] 238

``` r
# double check by adding all of the 1 and 0s from both alleles together
sum(samples$L116_allele_1[!is.na(samples$L116_allele_1)] == 1) +
  sum(samples$L116_allele_2[!is.na(samples$L116_allele_2)] == 1) +
  sum(samples$L116_allele_1[!is.na(samples$L116_allele_1)] == 0) +
  sum(samples$L116_allele_2[!is.na(samples$L116_allele_2)] == 0)
```

    ## [1] 238

``` r
# also 238

# Ratio of allele = 1 / allele = not NA
( sum(samples$L116_allele_1[!is.na(samples$L116_allele_1)] == 1) +
  sum(samples$L116_allele_2[!is.na(samples$L116_allele_2)] == 1) ) /
  ( sum(!is.na(samples$L116_allele_1)) + sum(!is.na(samples$L116_allele_2)) )
```

    ## [1] 0.3403361

``` r
# Repeat using tidyverse group_by and summarize
usingbasic <- samples %>%
  group_by(`Mitochondrial Marker`) %>%
  summarize(L116 = ( sum(L116_allele_1[!is.na(L116_allele_1)] == 1) +                  # allele 1 = 1
                      sum(L116_allele_2[!is.na(L116_allele_2)] == 1) ) /               # allele 2 = 1
                      ( sum(!is.na(L116_allele_1)) + sum(!is.na(L116_allele_2)) ),     # allele 1 & 2 is 1 or 0
            L148 = ( sum(L148_allele_1[!is.na(L148_allele_1)] == 1) +  
                      sum(L148_allele_2[!is.na(L148_allele_2)] == 1) ) /       
                      ( sum(!is.na(L148_allele_1)) + sum(!is.na(L148_allele_2)) ),
            L203 = ( sum(L203_allele_1[!is.na(L203_allele_1)] == 1) +  
                      sum(L203_allele_2[!is.na(L203_allele_2)] == 1) ) /       
                      ( sum(!is.na(L203_allele_1)) + sum(!is.na(L203_allele_2)) ),
            L430 = ( sum(L430_allele_1[!is.na(L430_allele_1)] == 1) +  
                      sum(L430_allele_2[!is.na(L430_allele_2)] == 1) ) /       
                      ( sum(!is.na(L430_allele_1)) + sum(!is.na(L430_allele_2)) ),
            L447 = ( sum(L447_allele_1[!is.na(L447_allele_1)] == 1) +  
                      sum(L447_allele_2[!is.na(L447_allele_2)] == 1) ) /       
                      ( sum(!is.na(L447_allele_1)) + sum(!is.na(L447_allele_2)) ))
       
kable(usingmean, digits = 4)
```

| Mitochondrial Marker |   L116 |   L148 |   L203 |   L430 |   L447 |
| :------------------- | -----: | -----: | -----: | -----: | -----: |
| similis              | 0.1356 | 0.0410 | 0.0846 | 0.4237 | 0.0985 |
| solidissima          | 0.5417 | 0.7097 | 0.9452 | 0.9133 | 0.9803 |

``` r
kable(usingbasic, digits = 4)
```

| Mitochondrial Marker |   L116 |   L148 |   L203 |   L430 |   L447 |
| :------------------- | -----: | -----: | -----: | -----: | -----: |
| similis              | 0.1356 | 0.0410 | 0.0846 | 0.4237 | 0.0985 |
| solidissima          | 0.5417 | 0.7097 | 0.9452 | 0.9133 | 0.9803 |

##### Using Mean to Calculate Frequency

  - Mitochondrial and each locus allele frequencies over all samples.
  - Allele frequencies of each locus splitting into two populations:
    mitochondrial similis and solidissima.
  - Mitochondrial and locus frequencies at sites and locations.

<!-- end list -->

``` r
#To calculate Solidissima Frequency
# ( Mean(allele1) + Mean(allele2) )
# Because 0,1 0,0 as alleles = 1 solidissima / 4 total = mean = allele freq

#Frequencies of Solidissimas in each populations (ideally 1 in solidissima, and 0 in similis)

byspecies <- samples %>%
  group_by(`Mitochondrial Marker`) %>%
  summarize(L116 = (mean(L116_allele_1[!is.na(L116_allele_1)]) + mean(L116_allele_2[!is.na(L116_allele_2)]))/2, 
            L148 = (mean(L148_allele_1[!is.na(L148_allele_1)]) + mean(L148_allele_2[!is.na(L148_allele_2)]))/2, 
            L203 = (mean(L203_allele_1[!is.na(L203_allele_1)]) + mean(L203_allele_2[!is.na(L203_allele_2)]))/2,
            L430 = (mean(L430_allele_1[!is.na(L430_allele_1)]) + mean(L430_allele_2[!is.na(L430_allele_2)]))/2,
            L447 = (mean(L447_allele_1[!is.na(L447_allele_1)]) + mean(L447_allele_2[!is.na(L447_allele_2)]))/2)

#Formatting
colnames(byspecies)[1] <- ""
byspecies <- byspecies[-1]
rownames(byspecies) <- c("Similis", "Solidissima")
```

**Solidissima Allele Frequency** in each population (as dictated by
mitochondrial marker). A 1.0 for Solidissima means that all Solidissima
samples according to the mitochondrial marker are homozygous
Solidissima, and a 0.0 under Similis means the same for Similis.

|             |  L116 |  L148 |  L203 |  L430 |  L447 |
| ----------- | ----: | ----: | ----: | ----: | ----: |
| Similis     | 0.136 | 0.041 | 0.085 | 0.424 | 0.098 |
| Solidissima | 0.542 | 0.710 | 0.945 | 0.913 | 0.980 |

These values are **significantly** different than what I calculated
before. **Please check my math.** These confirm similar relative
patterns in effectiveness between the primers (L116 is still bad and
L447 is still the best), but all have far closer values for similis
detection to expected.

For comparison:

|             |  L116 |  L148 |  L203 |  L430 |  L447 |
| ----------- | ----: | ----: | ----: | ----: | ----: |
| Similis     | 0.575 | 0.528 | 0.569 | 0.750 | 0.581 |
| Solidissima | 0.771 | 0.855 | 0.971 | 0.952 | 0.989 |

These were the values that lead us to question if the transcriptome used
was important. For clarity, below I have shown difference from ideal.
absolute(expected frequency - observed frequency)

``` r
diverge <- byspecies
diverge[2,] <- 1 - diverge[2,]
diverge <- diverge %>%
  rbind(c(diverge$L116[1]-diverge$L116[2],diverge$L148[1]-diverge$L148[2],diverge$L203[1]-diverge$L203[2],diverge$L430[1]-diverge$L430[2],diverge$L447[1]-diverge$L447[2]))
#Similis - Solidissima (positive numbers mean similis is being detected less - higher divergence from expected)
diverge <- full_join(diverge, t_freq)
diverge[5,] <- 1 - diverge[5,]
diverge <- diverge %>%
  rbind(c(diverge$L116[4]-diverge$L116[5],diverge$L148[4]-diverge$L148[5],diverge$L203[4]-diverge$L203[5],diverge$L430[4]-diverge$L430[5],diverge$L447[4]-diverge$L447[5]))

rownames(diverge) <- c("New Sim (O-E)", " New Sol (O-E)", "Difference in Detection", "Old Sim (O-E)", "Old Sol (O-E)", "Old Difference in Detection")

kable(diverge, digits = 2)
```

|                             |   L116 |   L148 | L203 | L430 | L447 |
| --------------------------- | -----: | -----: | ---: | ---: | ---: |
| New Sim (O-E)               |   0.14 |   0.04 | 0.08 | 0.42 | 0.10 |
| New Sol (O-E)               |   0.46 |   0.29 | 0.05 | 0.09 | 0.02 |
| Difference in Detection     | \-0.32 | \-0.25 | 0.03 | 0.34 | 0.08 |
| Old Sim (O-E)               |   0.57 |   0.53 | 0.57 | 0.75 | 0.58 |
| Old Sol (O-E)               |   0.23 |   0.15 | 0.03 | 0.05 | 0.01 |
| Old Difference in Detection |   0.35 |   0.38 | 0.54 | 0.70 | 0.57 |

Note here than in rows 3 and 6, *postive* values for difference in
detection mean that similis is being detected less well at this locus
than solidissima. In the newly calculated values, there is a fairly even
split between positive and negative values (to me indicating that there
is not a bias against similis), whereas in the older differences in
detection values all are postive and very large (indicating the
anti-similis bias that we previously concluded).

I have looked through my excel calculations from earlier and have still
not found the error yet but clearly this was a better way to do the
calculation if you agree with my math.

### By population

``` r
bysite <- samples %>%
  group_by(`Site`) %>%
  summarize(L116 = (mean(L116_allele_1[!is.na(L116_allele_1)]) + mean(L116_allele_2[!is.na(L116_allele_2)]))/2, 
            L148 = (mean(L148_allele_1[!is.na(L148_allele_1)]) + mean(L148_allele_2[!is.na(L148_allele_2)]))/2, 
            L203 = (mean(L203_allele_1[!is.na(L203_allele_1)]) + mean(L203_allele_2[!is.na(L203_allele_2)]))/2,
            L430 = (mean(L430_allele_1[!is.na(L430_allele_1)]) + mean(L430_allele_2[!is.na(L430_allele_2)]))/2,
            L447 = (mean(L447_allele_1[!is.na(L447_allele_1)]) + mean(L447_allele_2[!is.na(L447_allele_2)]))/2)
#Set NAs to be a - to reduce clutter
options(knitr.kable.NA = '-')
kable(bysite, digits=1)
```

| Site     | L116 | L148 | L203 | L430 | L447 |
| :------- | ---: | ---: | ---: | ---: | ---: |
| A - ARC  |   \- |   \- |  1.0 |  1.0 |  1.0 |
| B - ICO  |  0.5 |  1.0 |  1.0 |  1.0 |  1.0 |
| BLP      |  0.3 |  1.0 |  0.9 |  0.9 |  1.0 |
| C - MVSG |  0.0 |  0.0 |  0.0 |  0.3 |  0.0 |
| CUP      |  0.3 |  1.0 |  1.0 |  0.9 |  1.0 |
| GA       |  0.1 |  0.0 |  0.0 |  0.3 |  0.1 |
| GBE      |  0.9 |  0.2 |  1.0 |  1.0 |  1.0 |
| GLD      |  0.2 |  0.0 |  0.3 |  0.2 |  0.0 |
| MCX      |  0.8 |  1.0 |  1.0 |  1.0 |  1.0 |
| MW       |  0.2 |  0.1 |  0.1 |  0.6 |  0.1 |
| PEC      |  0.0 |  0.0 |  0.2 |  0.0 |  0.2 |
| PT       |  0.9 |  1.0 |  1.0 |  0.8 |  1.0 |
| SHN      |  0.7 |  0.0 |  0.9 |  0.9 |  1.0 |
| TEST     |   \- |   \- |  0.8 |   \- |  0.5 |

Showing Solidissima allele frequency (a 1.0 means that all those were
solidissima, a 0.0 means all similis)

These are condensed into broader areas below.

| Location             | L116 | L148 | L203 | L430 | L447 |
| :------------------- | ---: | ---: | ---: | ---: | ---: |
| George’s Bank        |  0.9 |  0.2 |  1.0 |  1.0 |  1.0 |
| Georgia              |  0.1 |  0.0 |  0.0 |  0.3 |  0.1 |
| Massachusetts        |  0.3 |  0.3 |  0.2 |  0.6 |  0.3 |
| Northern Long Island |  0.1 |  0.0 |  0.4 |  0.1 |  0.2 |
| Reitsman Hatchery    |  0.1 |  0.2 |  0.4 |  0.6 |  0.4 |
| Southern Long Island |  0.4 |  0.7 |  1.0 |  0.9 |  1.0 |

  - Where George’s Bank is just GBE - expected solidissima (1.0).
  - Georgia is just GA - expected similis (0.0).
  - Massachusetts, PT and MW - mixed.
  - Northern Long Islang - mostly similis (0.0).
  - Southern Long Island - mostly solidissima (1.0).
  - Reisman Hatery - mixed, more solidissima.
