library(ggplot2)
library(tidyverse)
library(knitr)

DP <- read.delim("filter_stats_DP.txt", stringsAsFactors = F, header = TRUE)
DP[2001,3] <-2000
DP %>%
  ggplot() +
  geom_line(aes(x=X.3.bin, y=X.6.number.of.sites), stat="identity") +
  ylim(0,1000) +
  theme_classic()

DP[,3] <-as.numeric(DP[,3])
DP[2001,3] <-2000

DP_sites <- DP %>% 
  group_by(gr=cut(X.3.bin, breaks= seq(0, 2000, by = 50)) ) %>% 
  summarise(total_sites=sum(X.6.number.of.sites))

DP_sites %>%
  ggplot() +
  geom_bar(aes(x=gr, y=total_sites), stat="identity") +
  scale_x_discrete(name ="Depth", 
                   breaks=c("0","500","1000","1500","2000") ) + 
  scale_y_continuous(name="Number of Sites") +
  theme_classic()


QUAL <- read.delim("filter_stats_QUALITY.txt", stringsAsFactors = F, header = TRUE)
QUAL %>%
  ggplot() +
  geom_line(aes(x=X.3.Quality, y=X.4.number.of.SNPs), stat="identity") +
  scale_y_continuous(name="Number of SNPs", trans='log10') + 
  scale_x_continuous(name="Quality") +
  theme_classic()
  











  
  
  scale_x_discrete(name ="Depth", 
                   breaks=c("0","500","1000","1500","2000") ) + 
  theme_classic()
  
  

scale_x_continuous(name="Depth of SNP site", limits=c(0, 2000)) +
  theme_classic()


scale_x_discrete(name ="Depth", 
                 limits=c("0","500","1000","1500","2000")) + 
