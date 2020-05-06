Figure Creation for Quarterly Report
================

``` r
library(tidyverse)
library(knitr)
```

## GitHub Documents

Length Histogram

``` r
Spislengths <- read.csv("./SpisulaR_2020May6.csv")
colnames(Spislengths) <- c("Old_ID","DNAPlate", "PlateWell", "ID", "Label", "Add","ShellStatus","CollectionTeamShell","HingeHeightmm","ShellLengthmm", "EstimateShellLength","Extracted","Species")

Spislengths$Species <- recode(Spislengths$Species, "failed to amplify; needs retest" = 'Not Yet Identified', "Failed to amplify; needs retest" = 'Not Yet Identified', 'S. similis'='S. similis', 'S. solidisima'= 'S. solidissima')

Spislengths %>%
  select(ID, ShellLengthmm, HingeHeightmm, EstimateShellLength, Species) %>%
  filter(EstimateShellLength > 0) %>%
  ggplot() +
  geom_histogram(aes(x=EstimateShellLength, fill=Species), bins=20) +
  labs(y=' ', x = "Estimate for Shell Length (mm)") +
  theme_minimal()
```

![](figureCreationQR_files/figure-gfm/cars-1.png)<!-- -->