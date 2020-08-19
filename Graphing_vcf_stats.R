library(ggplot2)
library(tidyverse)
library(knitr)

o_DP <- read.delim("Work/vcf_stats/original_DP_stats.txt", stringsAsFactors = F, header = TRUE)
p_DP <- read.delim("Work/vcf_stats/pilot_DP_stats.txt", stringsAsFactors = F, header = TRUE)
o_Q <- read.delim("Work/vcf_stats/original_QUAL_stats.txt", stringsAsFactors = F, header = TRUE)
p_Q <- read.delim("Work/vcf_stats/pilot_QUAL_stats.txt", stringsAsFactors = F, header = TRUE)
o_DP[2001,3] <-2000
o_DP[,3] <-as.numeric(o_DP[,3])
o_Q[,3] <-as.numeric(o_Q[,3])
p_DP[,3] <-as.numeric(p_DP[,3])
p_Q[,3] <-as.numeric(p_Q[,3])


Q <- cbind(p_Q, o_Q)
names(Q) <- c("QUAL","p_ID","p_QUAL","p_SNPs","p_trans1","p_trans2","p_indel","QUAL2","o_ID","o_QUAL","o_SNPs","o_trans1","o_trans2","o_indel")
DP <- merge(p_DP, o_DP,by=c("X.3.bin"), all=TRUE)
names(DP)<-c("DEPTH","name1","id1","p_genotypes","p_f_geno","p_sites","p_f_sites","name2","id2","o_genotypes","o_f_geno","o_sites","o_f_sites")



colors <- c("Transcript" = "salmon", "Pilot" = "lightblue")
Q %>%
  ggplot() +
  geom_line(aes(x=p_QUAL, y=o_SNPs, color="Transcript"), stat="identity") +
  geom_line(aes(x=p_QUAL, y=p_SNPs, color="Pilot"), stat="identity") +
  ylim(0,6000) +
  labs(x = "Quality",
       y = "SNP sites",
       color = "Legend") +
  theme_classic()

DP %>%
  ggplot() +
  geom_line(aes(x=DEPTH, y=o_sites, color="Transcript"), stat="identity") +
  geom_line(aes(x=DEPTH, y=p_sites, color="Pilot"), stat="identity") +
  ylim(0,6000) +
  labs(x = "Depth",
       y = "SNP sites",
       color = "Legend") +
  theme_classic()




# Example of getting manual named colors in gg plot with legend
colors <- c("Sepal Width" = "blue", "Petal Length" = "red", "Petal Width" = "orange")
ggplot(iris, aes(x = Sepal.Length)) +
  geom_line(aes(y = Sepal.Width, color = "Sepal Width"), size = 1.5) +
  geom_line(aes(y = Petal.Length, color = "Petal Length"), size = 1.5) +
  geom_line(aes(y = Petal.Width, color = "Petal Width"), size = 1.5) +
  labs(x = "Year",
       y = "(%)",
       color = "Legend") +
  scale_color_manual(values = colors)
