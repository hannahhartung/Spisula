Allele Frequency Calculations
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

rownames(diverge) <- c("New Sim |O-E|", " New Sol |O-E|", "Difference in Detection", "Old Sim |O-E|", "Old Sol |O-E|", "Old Difference in Detection")

kable(diverge)
```

|                             |        L116 |        L148 |      L203 |      L430 |      L447 |
| --------------------------- | ----------: | ----------: | --------: | --------: | --------: |
| New Sim |O-E|               |   0.1355932 |   0.0409836 | 0.0846154 | 0.4237288 | 0.0984848 |
| New Sol |O-E|               |   0.4583333 |   0.2903226 | 0.0547945 | 0.0866667 | 0.0197368 |
| Difference in Detection     | \-0.3227401 | \-0.2493390 | 0.0298209 | 0.3370621 | 0.0787480 |
| Old Sim |O-E|               |   0.5750000 |   0.5282258 | 0.5688406 | 0.7500000 | 0.5809859 |
| Old Sol |O-E|               |   0.2291667 |   0.1451613 | 0.0289855 | 0.0477941 | 0.0105634 |
| Old Difference in Detection |   0.3458333 |   0.3830645 | 0.5398551 | 0.7022059 | 0.5704225 |

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
