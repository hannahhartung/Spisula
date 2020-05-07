library(tidyverse)
library(knitr)

Sim <- tribble(
  ~c, ~p, ~stuff, ~pf,
  "chrom1", 2, "apple","P",
  "chrom1", 3, "apple","P",
  "chrom4", 3, "apple","P",
  "chrom2", 3, "orrange","P"
) 


Sol <- tribble(
  ~c, ~p, ~stuff,
  "chrom1", 1, "apple",
  "chrom1", 2, "apple",
  "chrom4", 3, "apple",
  "chrom2", 1, "orrange"
) 

#Make them dataframes
Sol <- as.data.frame(Sol)
Sim <- as.data.frame(Sim) x

for (i in 1:nrow(Sol)) {
  chrom <- Sol[i,1]
  pos <- Sol[i,2]
  for (j in 1:nrow(Sim)) {
    if (Sim[j,1]==chrom & Sim[j,2]==pos) {
      Sim[j,4] <- "F"
    }
  }
}

Sim %>%
  filter(pf == "P")
