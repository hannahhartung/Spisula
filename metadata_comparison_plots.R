####  Ploting Solidissima Inshore vs Out  ####
#### And Coast vs Shelf from that one PCA ####
library(ggplot2)

#Date: Jan 8th
sol_lengths <- as.data.frame(tmp<-read.table("ShellLengthSolidis.txt", header=TRUE))

#Updated Mar 2022 with some info from which samples were given to Matt
sol_lengths <- as.data.frame(tmp<-read.table("ShellLengthSolidis_mod0314.txt", header=TRUE))
#Updated July 2022 after more aging
sol_lengths <- as.data.frame(tmp<-read.table("ShellLengthSolidis_mod0724.txt", header=TRUE))

#Plot length by age
ggplot(sol_lengths) +
  geom_point(aes(x=AgeUpdate, y=Length_mmUpdate, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering and Site of Collection \n(SsLI data from update)", x="Age (years)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) +
  theme_bw()

#Plot growth function
#from: https://danstich.github.io/we-r-nycafs/fishStats.html
library(FSAdata) # for data
library(dplyr)   # for filter(), mutate()
library(ggplot2)
library(tidyverse)
library(FSA)
library(nlstools)
library(plotrix)

vbmod <- tl ~ Linf * (1 - exp(-K * (age - t0)))
#Where t1 and the ade are the collumns in your dataframe
vbmod <- Length_mmUpdate ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
starts <- vbStarts(formula = Length_mmUpdate ~ AgeUpdate, data = sol_lengths)
mymod <- nls(vbmod, data = sol_lengths, start = starts)
summary(mymod) #Mine did not have significance on t0 while theirs did
  #Calculates estimates for parameters
pred <- predict(mymod)
vbO <- vbFuns("typical")
vb_fit <- nls(Length_mmUpdate ~ vbO(AgeUpdate,Linf,K, t0), data=sol_lengths, start=starts)
boot_fit <- nlsBoot(vb_fit) #No working but I don't need 95% conf ints

#From: https://sizespectrum.org/mizer/reference/plotGrowthCurves.html
#Not for my purpose

#I likely need to do this separately for each group
#Get parameters for each
lengths_IC <- sol_lengths %>%
  filter(ShoreDistance=="Inshore") %>%
  filter(PCACatagory=="Coast")
lengths_IS <- sol_lengths %>%
  filter(ShoreDistance=="Inshore") %>%
  filter(PCACatagory=="Shelf")
lengths_OS <- sol_lengths %>%
  filter(ShoreDistance=="Offshore") %>%
  filter(PCACatagory=="Shelf")
vbmod <- Length_mmUpdate ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
starts_IC <- vbStarts(formula = Length_mmUpdate ~ AgeUpdate, data = lengths_IC)
mymod_IC <- nls(vbmod, data = lengths_IC, start = starts_IC)
summary(mymod_IC)
starts_IS <- vbStarts(formula = Length_mmUpdate ~ AgeUpdate, data = lengths_IS)
mymod_IS <- nls(vbmod, data = lengths_IS, start = starts_IS)
summary(mymod_IS)
starts_OS <- vbStarts(formula = Length_mmUpdate ~ AgeUpdate, data = lengths_OS)
mymod_OS <- nls(vbmod, data = lengths_OS, start = starts_OS)
summary(mymod_OS)

#Can I just plot a function on my r graph using those parameters?
base + geom_function(fun = ~ 0.5*exp(-abs(.x)))

#For inshore coast
geom_function(fun = ~ 0.5*exp(-abs(.x)))
fun = ~145.069*(1-EXP(-0.2375*(.x+1.90632)))
# Using a custom named function
f <- function(x) 0.5*exp(-abs(x))
f_IC <- function(x) 145.069*(1-exp(-0.2375*(x+1.90632)))
base <-
  ggplot() +
  xlim(0, 25)
base + geom_function(fun = f_IC) + geom_function(fun = f_OS)

f_IS <- function(x) 197.7921*(1-exp(-0.1585*(x+0.8733)))
f_OS <- function(x) 134.5084*(1-exp(-0.4801*(x-0.6519)))

ggplot(sol_lengths) +
  geom_point(aes(x=AgeUpdate, y=Length_mmUpdate, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering and Site of Collection \n(SsLI data updateed 3/14)", x="Age (years)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) +
  theme_bw() + geom_function(fun = f_IC, aes(color = "Coast", linetype = "Inshore"))  +
  geom_function(fun = f_IS, aes(color = "Shelf", linetype = "Inshore")) + geom_function(fun = f_OS, aes(color = "Shelf", linetype = "Offshore"))


geom_function(aes(colour = "normal"), fun = dnorm) +
  geom_function(aes(colour = "t, df = 1"), fun = dt, args = list(df = 1))

#Repeat but just separate by Genotype, not connection location
#Date: Mar 22nd 2022
lengths_A <- sol_lengths %>%
  filter(PCACatagory=="Coast")
lengths_B <- sol_lengths %>%
  filter(PCACatagory=="Shelf")
vbmod <- Length_mmUpdate ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
starts_A <- vbStarts(formula = Length_mmUpdate ~ AgeUpdate, data = lengths_A)
mymod_A <- nls(vbmod, data = lengths_A, start = starts_A)
summary(mymod_A)
#145.80861, 0.25648, -0.29142
    #old #145.06897, 0.23750, -1.90632
starts_B <- vbStarts(formula = Length_mmUpdate ~ AgeUpdate, data = lengths_B)
mymod_B <- nls(vbmod, data = lengths_B, start = starts_B)
summary(mymod_B)
#135.14065, 0.48754, 0.73684
    #old #135.1262, 0.4920, 0.7346
#f_A <- function(x) 145.06897*(1-exp(-0.23750*(x+1.90632)))
#f_B <- function(x) 135.1262*(1-exp(-0.4920*(x-0.7346)))

f_A <- function(x) 145.80861*(1-exp(-0.25648*(x+0.29142)))
f_B <- function(x) 135.14065*(1-exp(-0.48754*(x-0.73684)))


ggplot(sol_lengths) +
  geom_point(aes(x=AgeUpdate, y=Length_mmUpdate, color=PCACatagory)) +
  labs(title="Solidissima Growth by OTU Genotype", x="Age (years)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) +
  theme_bw() + geom_function(fun = f_A, aes(color = "Coast"))  +
  geom_function(fun = f_B, aes(color = "Shelf")) +
# Modify legend titles
  labs(color = "OTU") +
  scale_color_manual(labels = c("A", "B"), values = c("#F8766D", "#00BFC4"))
#Have to manually set colors
library(scales)
show_col(hue_pal()(4))

#Mar 2022
#Redo Depth Collection Graphs with all data
####Depth
ggplot(sol_lengths) +
  geom_density(aes(x=Depth_ft, color=PCACatagory)) +
  labs(title="Solidissima Samples Collection Depth", x="Collection Depth (feet)") +
  theme_bw() + labs(color='Genotype Category', linetype="Collection Distance") +
  labs(color = "Genotype") +
  scale_color_manual(labels = c("A", "B"), values = c("#F8766D", "#00BFC4"))
####Size vs Depth
ggplot(sol_lengths) +
  geom_histogram(aes(x=Length_mmUpdate, fill=PCACatagory),color="white") +
  labs(title="Coast vs Shelf Genetic Clutstering and Site of Collection \n(for 2012 Long Island Samples - Metadata Table)", x="Shell Length (mm)") +
  theme_bw() + labs(fill='Genotype Category')
ggplot(sol_lengths) +
  geom_point(aes(x=Depth_ft, y=Length_mmUpdate, color=PCACatagory)) +
  labs(title="Solidissima Shell Length from Collection Depths", x="Depth of Collection (ft)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) +
  theme_bw()
#Add line
ggplot(sol_lengths) +
  geom_point(aes(x=Depth_ft, y=Length_mmUpdate, color=PCACatagory)) +
  labs(title="Solidissima Shell Length from Collection Depths", x="Depth of Collection (ft)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) + geom_smooth(aes(x=Depth_ft, y=Length_mmUpdate, color=PCACatagory), method='lm', se=F, size=0.65, linetype = 2) +
  theme_bw() +
  labs(color = "Genotype") +
  scale_color_manual(labels = c("A", "B"), values = c("#F8766D", "#00BFC4"))
#Add line inshore/outshore - graph colors from multiple variables
ggplot(sol_lengths) +
  geom_point(aes(x=Depth_ft, y=Length_mmUpdate, color=PCACatagory:ShoreDistance)) +
  labs(title="Solidissima Shell Length from Collection Depths", x="Depth of Collection (ft)", y="Shell Length (mm)") +
  geom_smooth(aes(x=Depth_ft, y=Length_mmUpdate, color=PCACatagory:ShoreDistance), method='lm', se=F, size=0.65, linetype = 2) +
  theme_bw() +
  labs(color = "Genotype") +
  scale_color_manual(labels = c("A", "B Found Near Shore","B Found Off Shore"), values = c("#F8766D", "#00BFC4","#92d6ac"))


#Date: Mar 23rd
#Compare length to hinge height
sol_heightlengths <- as.data.frame(tmp<-read.csv("Work/length_height_sol_SLI_NAN.csv", header=TRUE))
ggplot(sol_heightlengths) +
  geom_point(aes(x=Height_mm, y=Length_mm, color=Pop)) +
  labs(x="Other Metric (mm)", y="Shell Length (mm)") + #scale_shape_manual(values=c(16, 3)) +
  theme_bw()





#





#Original Graphs
ggplot(sol_lengths) +
  geom_density(aes(x=Length_mm, linetype=ShoreDistance)) +
  labs(title="Solidissima Based on Site Distance from Shore \n(including NAN & 1999 - SsLI data from metatable)", x="Shell Length (mm)") +
  theme_bw()

ggplot(sol_lengths) +
  geom_density(aes(x=Length_mm, color=PCACatagory)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering  \n(including NAN & 1999 - SsLI data from metatable)", x="Shell Length (mm)") +
  theme_bw()
ggplot(sol_lengths) +
  geom_density(aes(x=Length_mmReport, linetype=ShoreDistance)) +
  labs(title="Solidissima Based on Site Distance from Shore  \n(including NAN & 1999 - SsLI data from report)", x="Shell Length (mm)") +
  theme_bw()
ggplot(sol_lengths) +
  geom_density(aes(x=Length_mm, color=PCACatagory, linetype=ShoreDistance)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering and Site of Collection \n(including NAN & 1999 - SsLI data from metatable)", x="Shell Length (mm)") +
  theme_bw()

ggplot(sol_lengths) +
  geom_density(aes(x=Length_mmReport, color=PCACatagory, linetype=ShoreDistance)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering and Site of Collection \n(including NAN & 1999 - SsLI data from report)", x="Shell Length (mm)") +
  theme_bw()

#Plot Age vs Length
ggplot(sol_lengths) +
  geom_point(aes(x=AgeTable, y=Length_mm, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering and Site of Collection \n(SsLI data from metatable)", x="Age (years)", y="Shell Length (mm)") + geom_smooth(aes(x=AgeTable, y=Length_mm, color=PCACatagory, linetype =ShoreDistance), method='lm', se=F, size=0.65) + scale_shape_manual(values=c(16, 6)) +
  theme_bw()

ggplot(sol_lengths) +
  geom_point(aes(x=AgeReport, y=Length_mmReport, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="Solidissima Based on Coast vs Shelf Clutstering and Site of Collection \n(SsLI data from report)", x="Age (years)", y="Shell Length (mm)") + geom_smooth(aes(x=AgeReport, y=Length_mmReport, color=PCACatagory, linetype =ShoreDistance), method='lm', se=F, size=0.65) + scale_shape_manual(values=c(16, 6)) +
  theme_bw()

#Plot Depth
ggplot(sol_lengths) +
  geom_density(aes(x=Depth_ft, color=PCACatagory, linetype=ShoreDistance)) +
  labs(title="Solidissima Collection Depth (Inshore Data Only)", x="Collection Depth (feet)") +
  theme_bw()


#Date: Jan 11th
#Plot just SsLI
SsLI_lengths <- sol_lengths %>% 
  filter(Population == "SsLI")
####Depth
ggplot(SsLI_lengths) +
  geom_density(aes(x=Depth_ft, color=PCACatagory, linetype=ShoreDistance),size=1.215) +
  labs(title="Long Island 2012 Samples Collection Depth", x="Collection Depth (feet)") +
  theme_bw() + labs(color='Genotype Category', linetype="Collection Distance")
####Size
ggplot(SsLI_lengths) +
  geom_histogram(aes(x=Length_mm, fill=PCACatagory),color="white") +
  labs(title="Coast vs Shelf Genetic Clutstering and Site of Collection \n(for 2012 Long Island Samples - Metadata Table)", x="Shell Length (mm)") +
  theme_bw() + labs(fill='Genotype Category')
ggplot(SsLI_lengths) +
  geom_histogram(aes(x=Length_mmReport, fill=PCACatagory),color="white") +
  labs(title="Coast vs Shelf Genetic Clutstering and Site of Collection \n(for 2012 Long Island Samples)", x="Shell Length (mm)") +
  theme_bw() + labs(fill='Genotype Category')
####Size by Age
ggplot(SsLI_lengths) +
  geom_point(aes(x=AgeTable, y=Length_mm, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="2012 Long Island Solidissima Size to Age \nBased on Coast vs Shelf Clutstering and Site of Collection \n(Data from metadata table)", x="Age (years)", y="Shell Length (mm)") + geom_smooth(aes(x=AgeTable, y=Length_mm, color=PCACatagory, linetype =ShoreDistance), method='lm', se=F, size=0.65) + scale_shape_manual(values=c(16, 6)) +
  theme_bw() + labs(color='Genotype Category', linetype="Collection Distance", shape="Collection Distance")
ggplot(SsLI_lengths) +
  geom_point(aes(x=AgeTable, y=Length_mmReport, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="2012 Long Island Solidissima Size to Age \nBased on Coast vs Shelf Clutstering and Site of Collection \n(Data from DEC Report)", x="Age (years)", y="Shell Length (mm)") + geom_smooth(aes(x=AgeTable, y=Length_mmReport, color=PCACatagory, linetype =ShoreDistance), method='lm', se=F, size=0.65) + scale_shape_manual(values=c(16, 6)) +
  theme_bw() + labs(color='Genotype Category', linetype="Collection Distance", shape="Collection Distance")

ggplot(sol_lengths) +
  geom_point(aes(x=AgeTable, y=Length_mm, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="Solidissima Size to Age \nBased on Coast vs Shelf Clutstering and Site of Collection \n(2012 Long Island Data from metadata table)", x="Age (years)", y="Shell Length (mm)") + geom_smooth(aes(x=AgeTable, y=Length_mm, color=PCACatagory, linetype =ShoreDistance), method='lm', se=F, size=0.65) + scale_shape_manual(values=c(16, 6)) +
  theme_bw() + labs(color='Genotype Category', linetype="Collection Distance", shape="Collection Distance")


ggplot(sol_lengths) +
  geom_point(aes(x=AgeTable, y=Length_mmReport, color=PCACatagory, shape=ShoreDistance)) +
  labs(title="Solidissima Size to Age \nBased on Coast vs Shelf Clutstering and Site of Collection", x="Age (years)", y="Shell Length (mm)") + #geom_smooth(aes(x=AgeTable, y=Length_mmReport, color=PCACatagory, linetype =ShoreDistance), method='lm', se=F, size=0.65) + 
  scale_shape_manual(values=c(16, 6)) +
  theme_bw() + labs(color='Genotype Category', linetype="Collection Distance", shape="Collection Distance")
#Try to remove just one of the geom_smooth
#Too hard, plot without lines, add in later and copy key in
