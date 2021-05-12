Quantifying Mollusc Kit Digestion Gels
================

## Tidy Data

Jan 14 2021

``` r
library(dplyr)
library(ggplot2)
```

``` r
dat <- read.table("/Users/hannah/gitHub/Spisula/Work/Quantification/Ss_MollsucKitDigest_1_12_20_r2_a.txt",header=F)
#make the second row the header after removing the first row
names(dat) <- as.matrix(dat[2, ])
dat <- dat[-2, ]
dat[] <- lapply(dat, function(x) type.convert(as.character(x)))
#rename column headers
colnames(dat)[which(names(dat) == "Offset(%)")] <- "Offset_percent"
colnames(dat)[which(names(dat) == "Offset(pix)")] <- "Offset_pix"
#remove the lines with names, replace with a collumn that has the lane name
dat <- filter(dat, Background != "----------")
#remove other headers
dat <- filter(dat, Background != "Background")
#add lane name
Lanes <- c(rep("MW_uncut",271),rep("MW_5",271),rep("MW_35",271),rep("MW_uncut2",271),rep("MW_1",271),rep("MW_05",271))
dat <- cbind(dat, Lanes)


head(dat)
```

    ##   Offset_percent Offset_pix Profile Background    Lanes
    ## 1         0.0000          0 63.0000    63.0000 MW_uncut
    ## 2         0.3704          1 65.0000    63.0116 MW_uncut
    ## 3         0.7407          2 64.0000    63.0233 MW_uncut
    ## 4         1.1111          3 65.0000    63.0349 MW_uncut
    ## 5         1.4815          4 66.0000    63.0465 MW_uncut
    ## 6         1.8519          5 66.0000    63.0581 MW_uncut

``` r
#create a for loop to do this across all of them,
#plot x=offset_percent, y= profile, and y2= background, grouped in plots by lanes 
```
