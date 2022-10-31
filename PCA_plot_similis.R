#### Library ####

##PCA for filter_b1215 sim and sol
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(dplyr)

#### Original Data Entry ####

#Because R does not like columns with names that are numbers, it renames
setwd("/Users/hannah/gitHub/Spisula/")
#no LD
#Date: Dec 15th
#o = only, a=all (aka all sequenced in that run)
sim_oLD <-as.data.frame(read.table("/Users/hannah/gitHub/Spisula/April_2022/LD_prune/LD02/SNP_simonly_LD02prune_forR2.txt", header=TRUE))
sim_a <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_sim_noLD_forR.txt", header=TRUE))
sim_o <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_simonly_noLD_forR.txt", header=TRUE))
sol_a <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_sol_noLD_forR.txt", header=TRUE))
sol_o <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_solonly_noLD_forR.txt", header=TRUE))

sol_LEFT <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_solonlyLEFT_noLD_forR.txt", header=TRUE))
addpopsLEFT<-as.data.frame(read.csv("addpopcode_LEFT.csv", header=TRUE),stringsAsFactors = FALSE)
sol_RIGHT <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_solonlyRIGHT_noLD_forR.txt", header=TRUE))
addpopsRIGHT<-as.data.frame(read.csv("addpopcode_RIGHT.csv", header=TRUE),stringsAsFactors = FALSE)

##### Summer Data #####

#solsimALL
#Actually all 500 samples (removing known problem samples = 484)
#Pruned for LD 
#Date: July 15th 2022
solsimALL <- as.data.frame(read.table("Solidissima/forR/SNP_solrefall_0528_LDprune_prep.txt", header=TRUE))
addpopsALL484 <- as.data.frame(read.csv("Solidissima/forR/addpopcodeALL484.csv", header=TRUE),stringsAsFactors = FALSE)
snpset <- solsimALL
addpopsproper <- addpopsALL484
addpopsALL484_hybrids <- as.data.frame(read.csv("Solidissima/forR/addpopcodeALL484_hybrids.csv", header=TRUE),stringsAsFactors = FALSE)
addpopsproper <- addpopsALL484_hybrids

#######
past(pause)

#### Confirming tree mismatch with PCA ####
#Date: Aug 1st
#Confirming tree mismatch with PCA
testsimsubsample <- as.data.frame(read.csv("Aug2022/forR/subsample_sim0801.csv", header=TRUE))
addpops_sub0801 <- as.data.frame(read.csv("Aug2022/forR/addpopcodesub0801.csv", header=TRUE),stringsAsFactors = FALSE)
snpset <- testsimsubsample
addpopsproper <- addpops_sub0801
testsimall <- as.data.frame(read.csv("Aug2022/forR/allsimref0801.csv", header=TRUE))
addpopsALLsim444 <- as.data.frame(read.csv("Aug2022/forR/addpopcodeALLsim444.csv", header=TRUE),stringsAsFactors = FALSE)
#Redo just similis
simo_prune0801 <- as.data.frame(read.csv("Aug2022/forR/SNP_simonly_LDprune0801_forR.csv", header=TRUE))
addpopssimo0801 <- as.data.frame(read.csv("Aug2022/forR/addpopcodesimo125_0801.csv", header=TRUE),stringsAsFactors = FALSE)

# Input
snpset <- simo_prune0801
addpopsproper <- addpopssimo0801

######
#Date: Dec 18th
#LD Pruned
sol_a_prune <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_solall_1fprad_1218_forR.txt", header=TRUE))
sol_o_prune <- as.data.frame(read.table("Solidissima/filter_b1215/forR/SNP_HWE_b1215_solonly_1fprad_1218_forR.txt", header=TRUE))
addpops388_REDUCECOLOR <- as.data.frame(read.csv("addpopcode388_REDUCECOLOR.csv", header=TRUE),stringsAsFactors = FALSE)
addpops388 <-  as.data.frame(read.csv("addpopcode388.csv", header=TRUE),stringsAsFactors = FALSE)
snpset <- sol_o_prune
addpopsproper <- addpops388
  #402 for a and 388 or 388_REDUCECOLOR for o
  #Removed PVT008andH08007 and NAN004 and 005

#FST Filtered
sol_o_FST085 <- as.data.frame(read.table("Solidissima/filter_b1215/forR/filteringFST/SNP_HWE_b1215_solonly_1fprad_1218_FST085_forR.txt", header=TRUE))
sol_o_FST05 <- as.data.frame(read.table("Solidissima/filter_b1215/forR/filteringFST/SNP_HWE_b1215_solonly_1fprad_1218_FST05_forR.txt", header=TRUE))
sol_o_FST01 <- as.data.frame(read.table("Solidissima/filter_b1215/forR/filteringFST/SNP_HWE_b1215_solonly_1fprad_1218_FST01_forR.txt", header=TRUE))
snpset <- sol_o_FST05
addpopsproper <- addpops388
  #Created from vcf already missing PVT008andH08007 and NAN004 and 005
######

#set correct addpops
addpopssimLD <- as.data.frame(read.csv("/Users/hannah/gitHub/Spisula/April_2022/LD_prune/LD02/addpopssim0622.csv", header=TRUE),stringsAsFactors = FALSE)
addpops143 <- as.data.frame(read.csv("addpopcode143.csv", header=TRUE),stringsAsFactors = FALSE)
addpops125 <- as.data.frame(read.csv("addpopcode125.csv", header=TRUE),stringsAsFactors = FALSE)
addpops125_MAMFC <- as.data.frame(read.csv("addpopcode125_MAMFC.csv", header=TRUE),stringsAsFactors = FALSE)
addpops391 <- as.data.frame(read.csv("addpopcode391.csv", header=TRUE),stringsAsFactors = FALSE)
addpops402 <- as.data.frame(read.csv("addpopcode402.csv", header=TRUE),stringsAsFactors = FALSE)
addpops402_MAMFC <- as.data.frame(read.csv("addpopcode402_MAMFC.csv", header=TRUE),stringsAsFactors = FALSE)
addpops402_MAMFC_BONUS <- as.data.frame(read.csv("addpopcode402_MAMFC_BONUS.csv", header=TRUE),stringsAsFactors = FALSE)
addpops402_MAMFC_BONUS2 <- as.data.frame(read.csv("addpopcode402_MAMFC_BONUS_2.csv", header=TRUE),stringsAsFactors = FALSE)
addpops389 <- as.data.frame(read.csv("addpopcode389.csv", header=TRUE),stringsAsFactors = FALSE) #actually now 385
addpops389_REDUCECOLOR <- addpops389 <- as.data.frame(read.csv("addpopcode389_REDUCECOLOR.csv", header=TRUE),stringsAsFactors = FALSE)
addpops389_REDUCECOLOR_BONUS <- addpops389 <- as.data.frame(read.csv("addpopcode389_REDUCECOLOR_BONUS.csv", header=TRUE),stringsAsFactors = FALSE)

snpset <- sol_o
addpopsproper <- addpops389_REDUCECOLOR_BONUS
addpopsproper <-  addpops402_MAMFC_BONUS2
snpset <- sim_o
addpopsproper <-  addpops125
snpset <- sol_a
snpset <- sol_a_prune
addpopsproper <- addpops125_MAMFC
#population factor levels
poplevs_sim <- c("GBE","CCB","SCC","NLI","SLI","SJF","GA")
poplevs_sol <- c("413","GBE","CCB","SCC","NAN","SIM","NLI","SLI","536","FNJ","SJF","GA")
c("X413","GBE","CCB","SCC","NAT","SIM","NLI","SLI","X536","FNJ","SJF")
poplevs_REDUCE <- c("GB","MA","NY","NJ")

#June 22nd
#Sim
poplevs_sim2 <- c("GA","MA","NY")
poplev <- poplevs_sim2
snpset <- sim_oLD
addpopsproper <- addpopssimLD

#Sol
#Reduce Color June
#lumping all offshore Federal samples (NJ, GB, NAN) into one light-colored symbol
addpops388_Rj <- as.data.frame(read.csv("addpopcode388_REDUCEjune.csv", header=TRUE),stringsAsFactors = FALSE)
snpset <- sol_o_prune
addpopsproper <- addpops388_Rj

addpops388_UR <- as.data.frame(read.csv("addpopcode388_UR.csv", header=TRUE),stringsAsFactors = FALSE)
addpopsproper <- addpops388_UR

#######START GRAPHING######
tmp_snpset <- snpset %>%
  dplyr::select(CHROM,POS) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":"))
locus_snpset <- snpset %>%
  dplyr::select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(LOCUS,everything())
t_snpset = setNames(data.frame(t(locus_snpset[,-1]),stringsAsFactors = FALSE), locus_snpset[,1])
t_snpset <- t_snpset[-c(1,2,3,4,5,6,7,8,9), ]
tmp2_snpset <- cbind.data.frame(t_snpset,addpopsproper)
i2_snpset <- tmp2_snpset %>%
  dplyr::select(X,POP,SITE,SPECIES,everything())
tmp3_snpset <- i2_snpset[-c(1:4)]
tmp3_snpset[is.na(tmp3_snpset)] <- 0 #fiddle with other numbers to confirm again that it has minor impact
tmp3_snpset<- as.data.frame(tmp3_snpset,stringsAsFactors = FALSE)
tmp5_snpset <- mutate_all(tmp3_snpset, function(x) as.numeric(as.character(x)))
tmp5_snpset[is.na(tmp5_snpset)] <- 0
frameit_snpset <- summarise_all(tmp5_snpset,n_distinct)
frameT_snpset <- which(frameit_snpset > 1)
tmp6_snpset <- tmp5_snpset[frameT_snpset]
#prepare for plot
#i2_snpset$POP <- factor(i2_snpset$POP, levels = c("X413","GBE","CCB","SCC","NAT","SIM","NLI","SLI","X536","FNJ","SJF","GA"))
i2_snpset$POP <- factor(i2_snpset$POP, levels = c("GB","NAT","CCB","SCC","MA","NY","NLI","SLI","NJ","GA","OTU-A"))
#Main one below:
#i2_snpset$POP <- factor(i2_snpset$POP, levels = c("OFF","GB","MA","CCB","SC","SCC", "NY","NJ", "MA similis", "NY similis","NLI", "SLI","GA", "NAT","Coast","Shelf","Similis"))
#i2_snpset$SPECIES <- factor(i2_snpset$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
pca_res_snpset <- prcomp(tmp6_snpset , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)

################# SAVING SNPSETS ########################
#Start saving pca_res_snpsets so I do not have to rerun all the time - then comment so I don't accidentally do it
####
#set_SNPs_mainsolidis_with_a_few_sim <- pca_res_snpset
#set_SNPs_JUSTsolidis <- pca_res_snpset
#set_All <- pca_res_snpset
#set_All_i2 <- i2_snpset
####
#Loading save
#pca_res_snpset <- set_All
#i2_snpset <- set_All_i2


##### MAIN PLOT #####
autoplot(pca_res_snpset, data = i2_snpset, color = 'POP', size=2) + theme_bw() #scale_shape_manual(values=c(16, 18, 3),"Species") +
  #scale_color_brewer(spectrum) #scale_color_viridis_d(option="viridis","Population") #+ scale_colour_manual(values=c("MA" = "#FDAE61", "NY" = "#ABDDA4", "GA" = "#FEE08B"),"Population")
#plot(x, main="FST 0.85 Filter - Solidissima")

#Plot sample IDs
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape=FALSE, label.size = 1) +
  scale_shape_manual(values=c(16, 18, 3),"Species") + theme_bw()

#Plot for MAMFC similis
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape ='SPECIES', size=2) +
  scale_shape_manual(values=c(16, 17),"Subspecies") +
  theme_bw() + scale_colour_brewer(palette = "Spectral","Population")

####Try to add probability ovals
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm', frame.alpha = 0.05) +
  scale_shape_manual(values=c(16, 17),"Collection Distance") +
  theme_bw() + scale_colour_brewer(palette = "Spectral","Population") + scale_fill_brewer(palette = "Spectral","Population")


#Date: July 15th
########### For Publication ###############
#Recreate the two species PCA in the same colors
#Just subspecies
autoplot(pca_res_snpset, data = i2_snpset,colour ='SPECIES', size=2) +
  theme_bw() + 
  scale_colour_manual(values=c("#95DB75", "#9A70F8"),"Subspecies")

#Divide solidissima
autoplot(pca_res_snpset, data = i2_snpset,colour ='SITE', size=2) +
  theme_bw() + 
  scale_colour_manual(values=c( "#FF6C67","#01BFC4","#95DB75"),"OTU")

#With Outlines - does not look good.
autoplot(pca_res_snpset, data = i2_snpset,fill ='SITE', colour ='SPECIES', size=2, shape = 21) +
  theme_bw() + 
  scale_fill_manual(values=c( "#FF6C67","#01BFC4","#95DB75"),"OTU") +
  scale_colour_manual(values=c("#95DB75", "#9A70F8"),"Subspecies")
#With Hybrids
autoplot(pca_res_snpset, data = i2_snpset,colour ='SITE', size=2) +
  theme_bw() + 
  scale_colour_manual(values=c("#FF6C67","#9A70F8","#01BFC4","#95DB75"),"OTU")


#Plot by Region (internal/supplemental use only)
regionbothsp<-c("#F251F3","#613CF0","#3BA4F4","#33EFFF","#40E650","#DDFF17","#FFE017","#FF8A17","#FF2E17")
#regionbothsp<-c("#F251F3","#1ECE00","#00CFAE","#33EFFF","#FF8A17","#DDFF17","#FFE017","#FF2E17","#613CF0")
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm', frame.alpha = 0.05) +
  theme_bw() + 
  scale_colour_manual(values=regionbothsp,"Region") + scale_fill_manual(values=regionbothsp,"Region")
autoplot(pca_res_snpset, data = i2_snpset,colour ='POP', size=2) +
  theme_bw() + 
  scale_colour_manual(values=regionbothsp,"Region")

#Plot Just Solidissima
regionsol<-c("#D0679C","#3D2C7D","#3BA4F4","#5BECAB","#5ABB64","#FFB357","#FFE017","#FF8A17","#FF2E17")
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2) +
  theme_bw() + 
  scale_colour_manual(values=regionsol,"Region")
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm', frame.alpha = 0.05) +
  scale_shape_manual(values=c(16, 17),"Collection Distance") +
  theme_bw() +
  scale_colour_manual(values=regionsol,"Region") +
  scale_fill_manual(values=regionsol,"Region")

#
#Color plan - south is orange
#GA = red
    #FF2E17
#NJ = orange
    #FF8A17
#NY = NLI = SLI = variations on yellow/llight green
    #FFE017 #souther
    #DDFF17 #nothr
#MA = dark green
    #1ECE00
#SCC = light blue
    #33EFFF
#CCB = dark teal blue
    #00CFAE
#NAT = purple blue
    #613CF0
#GB = bright pink
    #F251F3








#I have two colored as A but pairing with B, did I mislabel some?
autoplot(pca_res_snpset, data = i2_snpset,colour ='SITE', shape=FALSE, label.size = 1) +
  theme_bw() + 
  scale_colour_manual(values=c( "#FF6C67","#01BFC4","#A3D099"),"OTU")
#Ss4106_1 and 2

#2 Subspecies
#addpopsproper <- sol_a_prune
#addpopsproper <- sol_a
#addpopsproper <- addpops125_MAMFC
#
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm') +
  theme_bw() + scale_colour_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Population") + scale_fill_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Population")

autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm') +
  theme_bw() + scale_colour_manual(values=c("#63423F", "#01BFC3", "#F3AE60"),"Population") + scale_fill_manual(values=c("#63423F", "#01BFC3", "#F3AE60"),"Population") +
  theme(text = element_text(size = 12))


#No ovals
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2) +
  theme_bw() + scale_colour_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Genotype",labels=c("A","B","Similis")) + scale_fill_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Genotype",labels=c("A","B","Similis"))

#Just similis
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm') +
  theme_bw() + scale_colour_manual(values=pptcolors,"Population") + scale_fill_manual(values=pptcolors,"Population")#+ scale_colour_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Population") + scale_fill_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Population")
originalcolors <- c("#B3539E","#ABDDA4", "#FEE08B")
pptcolors <- c("#EE766D","#943EAE","#2AB7D6")

pptcolors <- c("#DE6A60","#0097C1","#E7C582")
"#452123"
"#637795"
"#6C6BA5"
"#3C2158"
"#E7C582"

#Just similis color shift
#GA = Purple, NY = Pink, SCC = Bright Pink
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm') +
  theme_bw() + scale_colour_manual(values=c("#B3539E","#FFBF65", "#FD625E"),labels = c("GA", "NY","MA"),"Population") + scale_fill_manual(values=c("#B3539E","#FFBF65","#FD625E"),labels = c("GA", "NY","MA"),"Population")#+ scale_colour_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Population") + scale_fill_manual(values=c("#F8776D", "#01BFC3", "#98D494"),"Population")

pptcolors <- c("#B3539E","#FFBF65","#FD625E")
pptcolors <- c("#AE3535","#FFAC00","#E7C582")
pptcolors <- c("#AE3535","#A66999","#E7C582")
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2, frame = TRUE, frame.type = 'norm') +
  theme_bw() + scale_colour_manual(values=pptcolors,labels = c("GA", "NY","MA"),"Population") + scale_fill_manual(values=pptcolors,labels = c("GA", "NY","MA"),"Population")


#Bonus
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape ='SPECIES', size=2) +
  scale_shape_manual(values=c(16, 17),"Subspecies") +
  theme_bw() + scale_color_manual(values=c("#D7191C", "#FDAE61", "#ABDDA4", "#2B83BA","#e69138","#93c47d"), "Population") + 
  guides(colour = guide_legend(override.aes = list(shape = 17)))

autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape ='SPECIES', size=2) +
  scale_shape_manual(values=c(16, 17),"Distance") + scale_color_manual(values=c("#D7191C", "#FDAE61", "#fde661", "#ABDDA4", "#2B83BA"), "Population") +
  guides(colour = guide_legend(override.aes = list(shape = 16))) +
  theme_bw() 

#Simils
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP', size=2) +
  theme_bw() + scale_colour_manual(values=c("MA" = "#FDAE61", "NY" = "#ABDDA4", "GA" = "#FEE08B"),"Population")
#Get same colors
#Display the Hex for colors in a brewer pallete
brewer.pal(n = 4, name = 'Spectral')
display.brewer.pal(n = 3, name = 'Spectral') #show them
#From: https://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html

#Special plot for making stuff super visible for sol right/left
#Date: Dec 18th
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape ='SPECIES', size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw() + scale_colour_manual(values=c("X413" = "#F8766D", "GBE" = "#CD9600", "CCB" = "#7CAE00", "SCC" = "steelblue4", "NAT" = "#00BFC4", "SLI" = "orangered4", "X536" = "#C77CFF", "FNJ" = "#FF61CC"))
#Print the brewer for all the other ones but then set SCC and SLI to black or maybe their border <- this doesn't work because they are set to be the shapes where the line is the whole thing, i'd have to make it shade 1 and then fill=POP
library(RColorBrewer)
#display.brewer.pal(n = 8)
#Not a brewer palette but the basic:
library(scales)
show_col(hue_pal()(8))
"#F8766D"
"#CD9600"
"#7CAE00"
"#00BE67"
"#00BFC4"
"#00A9FF"
"#C77CFF"
"#FF61CC"
#Plot sample IDs
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape=FALSE, label.size = 2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw()

#For FST01
#Remove FNJ028 and 029 just to make it look nice
#These also had to be removed before, I assume they are highly related
i2_snpset <- i2_snpset %>%
  filter(X != "FNJ_029") %>%
  filter(X != "FNJ_028")


#PLot other axes
autoplot(pca_res_snpset, data = i2_snpset, colour = 'POP',shape ='SPECIES', size=2,x=2,y=3) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw() + scale_color_viridis_d(option="viridis","Population")

#Sol_only
#Two samples are big outliers <- remove
#NAN_004 and NAN_005

#After doing that, there are two new outliers: 
#PVT08H07 and PVT008
#Now two more!
#FNJ_028
#FNJ_029

#Save prcomp to table so I can pick out which samples are what along PC1
str(pca_res_snpset)
pc.new<-cbind(pca_res_snpset,pca_res_snpset$x[,1:3])
pc_sol_o<- unclass(pc.new)[,1:3]
#Doesn't like to save proper so just print out the head and the tail and put in excel
head(unclass(pc.new)[,1:3],400)
tail(unclass(pc.new)[,1:3],60)

#NEW OUTLIERs for RIGHT
#MCX038, SSLI1113_4, SSLI1112_1, MCX_021, BLP_237, SSLI3109_1, SSLI2108_3



####pca for undedup ddoc####
library(ggplot2)
library(ggfortify)
library(tidyverse)

#filter1012
sim1012vcf <- read.table("Apop_struct_analysis/SNP_filter1012sf_noLD_simonly_forR1025.txt", header=TRUE)
svcf <- as.data.frame(sim1012vcf)
tmp <- svcf %>%
  dplyr::select(CHROM,POS) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":"))
locusvcf <- svcf %>%
  dplyr::select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(LOCUS,everything())
tvcf = setNames(data.frame(t(locusvcf[,-1]),stringsAsFactors = FALSE), locusvcf[,1])
tvcf <- tvcf[-c(1,2,3,4,5,6,7,8,9), ]
addpops125 <- as.data.frame(read.csv("addpopcode125.csv", header=TRUE),stringsAsFactors = FALSE)
tmp <- cbind.data.frame(tvcf,addpops125)
iris2 <- tmp %>%
  dplyr::select(X,POP,SITE,SPECIES,everything())
tmpdf <- iris2[-c(1:4)]
tmpdf[is.na(tmpdf)] <- 0
tmpdf<- as.data.frame(tmpdf,stringsAsFactors = FALSE)
tmpdf1 <- mutate_all(tmpdf, function(x) as.numeric(as.character(x)))
tmpdf1[is.na(tmpdf1)] <- 0
#remove columns with all constants
frameit <- summarise_all(tmpdf1,n_distinct) #not ideal
frameT <- which(frameit > 1) #list of columns that have more than one thing
tmpdf2 <- tmpdf1[frameT]
iris2$POP <- factor(iris2$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
iris2$SPECIES <- factor(iris2$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
pca_res <- prcomp(tmpdf2 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_res, data = iris2, colour = 'POP',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw()
#Plot sample IDs
autoplot(pca_res, data = iris2, colour = 'POP',shape=FALSE, label.size = 1) +
  scale_shape_manual(values=c(16, 18, 3),"Species")
#Plot other axes
autoplot(pca_res, data = iris2, colour = 'POP',size=2,x=2,y=3) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw()




#original, less filtered data
#simplifiedvcf <- read.table("dDocenthelp/vcfs_929/undedup_nmq30_forR.txt", header=TRUE)
#import vcf from filterDoc_105
simplifiedvcf <- read.table("dDocenthelp/JP_filter/SNP_JPfilter105a_forR.txt", header=TRUE)

svcf <- as.data.frame(simplifiedvcf)
tmp <- svcf %>%
  dplyr::select(CHROM,POS) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":"))
locusvcf <- svcf %>%
  dplyr::select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(LOCUS,everything())
tvcf = setNames(data.frame(t(locusvcf[,-1]),stringsAsFactors = FALSE), locusvcf[,1])
tvcf <- tvcf[-c(1,2,3,4,5,6,7,8,9), ]
#tmp <- tvcf %>%
#  select("107:10")

#addpops with all 150
#addpops <- as.data.frame(read.csv("addpopcode.csv", header=TRUE),stringsAsFactors = FALSE)
#tmp <- cbind.data.frame(tvcf,addpops)

#add pops with only 143
addpops143 <- as.data.frame(read.csv("addpopcode143.csv", header=TRUE),stringsAsFactors = FALSE)
tmp <- cbind.data.frame(tvcf,addpops143)

iris2 <- tmp %>%
  dplyr::select(X,POP,SITE,SPECIES,everything())
tmpdf <- iris2[5:2914]
tmpdf[is.na(tmpdf)] <- 0
tmpdf<- as.data.frame(tmpdf,stringsAsFactors = FALSE)
tmpdf1 <- mutate_all(tmpdf, function(x) as.numeric(as.character(x)))
tmpdf1[is.na(tmpdf1)] <- 0
#remove columns with all constants
frameit <- summarise_all(tmpdf1,n_distinct) #not ideal
frameT <- which(frameit > 1) #list of columns that have more than one thing
tmpdf2 <- tmpdf1[frameT]
iris2$POP <- factor(iris2$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
iris2$SPECIES <- factor(iris2$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
pca_res <- prcomp(tmpdf2 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_res, data = iris2, colour = 'POP', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  scale_color_viridis_d(option="plasma","Population") +
  theme_bw()

#Plot other PCA axes
autoplot(pca_res, data = iris2, colour = 'POP', shape='SPECIES',size=2,x=2,y=3) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  scale_color_viridis_d(option="plasma","Population")
#Plot sample IDs
autoplot(pca_res, data = iris2, colour = 'POP',shape=FALSE, label.size = 1) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  scale_color_viridis_d(option="plasma","Population")
#Plot just similis
justsim <- iris2 %>%
  dplyr::select(X,POP,SITE,SPECIES,everything()) %>%
  filter(SPECIES == "similis")
tmpdfsim <- justsim[5:2914]
tmpdfsim[is.na(tmpdfsim)] <- 0
tmpdfsim<- as.data.frame(tmpdfsim,stringsAsFactors = FALSE)
tmpdfsim1 <- mutate_all(tmpdfsim, function(x) as.numeric(as.character(x)))
tmpdfsim1[is.na(tmpdfsim1)] <- 0
pca_ressim <- prcomp(tmpdfsim1)
autoplot(pca_ressim, data = justsim, colour = 'POP') +
  theme_bw()
scale_color_viridis_d(option="viridis","Population")#+
autoplot(pca_ressim, data = justsim, colour = 'POP', shape=FALSE, label.size = 2) +
  theme_bw() #+
scale_color_viridis_d(option="viridis","Population")#+


#61	MW_016 was WAAY FAR OUT - remove for graph
tmpdfsim2 <- justsim[5:2914]
tmpdfsim2 <- tmpdfsim2[-61,]
tmpdfsim2[is.na(tmpdfsim2)] <- 0
tmpdfsim2<- as.data.frame(tmpdfsim2,stringsAsFactors = FALSE)
tmpdfsim3 <- mutate_all(tmpdfsim2, function(x) as.numeric(as.character(x)))
tmpdfsim3[is.na(tmpdfsim3)] <- 0
pca_ressim3 <- prcomp(tmpdfsim3)
autoplot(pca_ressim3, data = justsim[-61,], colour = 'POP') +
  theme_bw()
scale_color_viridis_d(option="viridis","Population")#+


### Compare Similis to dDoc and to Transcript ###
sim1012vcf <- read.csv("SimilisCompare2Transcriptome/SNP1012_dD_simonly_forR.csv", header=TRUE)
sim1129vcf <- read.csv("SimilisCompare2Transcriptome/SNP1129_T_simonly_forR.csv", header=TRUE)

Dsvcf <- as.data.frame(sim1012vcf)
Tsvcf <- as.data.frame(sim1129vcf)

locusDvcf <- Dsvcf %>%
  dplyr::select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(LOCUS,everything())
locusTvcf <- Tsvcf %>%
  dplyr::select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(LOCUS,everything())

Dtvcf = setNames(data.frame(t(locusDvcf[,-1]),stringsAsFactors = FALSE), locusDvcf[,1])
Dtvcf <- Dtvcf[-c(1,2,3,4,5,6,7,8,9), ]
Ttvcf = setNames(data.frame(t(locusTvcf[,-1]),stringsAsFactors = FALSE), locusTvcf[,1])
Ttvcf <- Ttvcf[-c(1,2,3,4,5,6,7,8,9), ]


addpops125 <- as.data.frame(read.csv("addpopcode125.csv", header=TRUE),stringsAsFactors = FALSE)
addpops121T <- as.data.frame(read.csv("addpopcode121T.csv", header=TRUE),stringsAsFactors = FALSE)

tmpT <- cbind.data.frame(Ttvcf,addpops121T)
tmpD <- cbind.data.frame(Dtvcf,addpops125)

irisT <- tmpT %>%
  dplyr::select(X,POP,SITE,SPECIES,everything())
tmpTdf <- irisT[-c(1:4)]
tmpTdf[is.na(tmpTdf)] <- 0
tmpTdf<- as.data.frame(tmpTdf,stringsAsFactors = FALSE)
tmpTdf1 <- mutate_all(tmpTdf, function(x) as.numeric(as.character(x)))
tmpTdf1[is.na(tmpTdf1)] <- 0
frameitT <- summarise_all(tmpTdf1,n_distinct)
frametT <- which(frameitT > 1)
tmpTdf2 <- tmpTdf1[frametT]

irisT$POP <- factor(irisT$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
irisT$SPECIES <- factor(irisT$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
pca_resT <- prcomp(tmpTdf2 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_resT, data = irisT, colour = 'POP',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw() + ggtitle("Transcriptome SNP Similis")
autoplot(pca_resT, data = irisT, colour = 'POP',shape=FALSE, label.size = 1) +
  scale_shape_manual(values=c(16, 18, 3),"Species")
    #Still has some sneaky sol <- remove
      #Remove from 


irisD <- tmpD %>%
  dplyr::select(X,POP,SITE,SPECIES,everything())
tmpDdf <- irisD[-c(1:4)]
tmpDdf[is.na(tmpDdf)] <- 0
tmpDdf<- as.data.frame(tmpDdf,stringsAsFactors = FALSE)
tmpDdf1 <- mutate_all(tmpDdf, function(x) as.numeric(as.character(x)))
tmpDdf1[is.na(tmpDdf1)] <- 0
frameitD <- summarise_all(tmpDdf1,n_distinct)
frametD <- which(frameitD > 1)
tmpDdf2 <- tmpDdf1[frametD]

irisD$POP <- factor(irisD$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
irisD$SPECIES <- factor(irisD$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
pca_resD <- prcomp(tmpDdf2 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_resD, data = irisD, colour = 'POP',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw() + ggtitle("dDocent SNP Similis")
autoplot(pca_resD, data = irisD, colour = 'POP',shape=FALSE, label.size = 1) +
  scale_shape_manual(values=c(16, 18, 3),"Species")



#remove columns with all constants
frameit <- summarise_all(tmpdf1,n_distinct) #not ideal
frameT <- which(frameit > 1) #list of columns that have more than one thing
tmpdf2 <- tmpdf1[frameT]
iris2$POP <- factor(iris2$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
iris2$SPECIES <- factor(iris2$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
pca_res <- prcomp(tmpdf2 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_res, data = iris2, colour = 'POP',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  theme_bw()
#Plot sample IDs
autoplot(pca_res, data = iris2, colour = 'POP',shape=FALSE, label.size = 1) +
  scale_shape_manual(values=c(16, 18, 3),"Species")

############# Haplotype Data##############
#Input genepop files for haplotype data
library(devtools)
library(hierfstat)
library(graph4lg)
library(vcfR)
library(adegenet)
library(tidyverse)

###Just similis
haps <- read.genepop("ima/haps_simONLY_0126.gen", ncode=3)
obj2 <- genind2genpop(haps)
X <- tab(haps, NA.method="mean")
tmp<-tab(haps)
pop(haps) <- rep(c('GA','NLI','SCC'),c(23,13,88)) #add populations
summary(haps)

pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(haps))
myCol <- transp(c("blue","red","green"),.7)[temp]
plot(pca1$li, col=myCol)
legend("topleft", pch=15, col=transp(c("blue","red", "green"),.7),
       leg=c("GA","NLI","SCC"), pt.cex=2)

###Both
haps <- read.genepop("ima/haps_simsol2sim.gen", ncode=3)
X <- tab(haps, NA.method="mean")
tmp<-tab(haps)
pop(haps) <- rep(c('CO','GA','NLI','SCC',"SH"),c(27,23,13,81,47)) #add populations
summary(haps)

pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(haps))
myCol <- transp(c("orange","blue","red","green","purple"),.7)[temp]
plot(pca1$li, col=myCol)
legend("bottomright", pch=15, col=transp(c("orange","blue","red","green","purple"),.7),
       leg=c("Coast","GA","NLI","SCC","Shelf"), pt.cex=2)

#

















####################original notes below###################


###adegenet###








####original notes####
library(ggplot2)
library(ggfortify)
library(tidyr)
df <- iris[1:4]
pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res, data = iris, colour = 'Species')
#So pull out the non numerical values out for the prcomp but then to plot with color use the other

#Let's get to a place where I set up the pca from the vcf

autoplot(pca_res, data = iris, colour = 'Species', label = TRUE, label.size = 3)

autoplot(pca_res, data = iris, colour = 'Species',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


simplifiedvcf <- read.table("vcfwithnanforR.txt", header=TRUE)

#https://github.com/nt246/NTRES6940-data-science/blob/f43638d0616d8f39c3fda6588a85b346a2b86a61/analysis/lesson10-tidy-data.md

#first mutate a collumn that is c___CHROM_p___POS so each loci has a unique identifier
#then make those the collumn names

svcf <- as.data.frame(simplifiedvcf)

tmp <- svcf %>%
  dplyr::select(CHROM,POS) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":"))

locusvcf <- svcf %>%
  dplyr::select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(LOCUS,everything())
#this is the same as the ID but since I was too lazy to fix the IDs in excel

#set the headers when I transpose
tvcf = setNames(data.frame(t(locusvcf[,-1]),stringsAsFactors = FALSE), locusvcf[,1])
tvcf <- tvcf[-c(1,2,3,4,5,6,7,8,9), ]
#remove old collumns
#add column with population

tmp <- tvcf %>%
  dplyr::select("107:10")
#write.csv(tmp,"addpopcode.csv", row.names = TRUE)
#add in POP, SITE and SPECIES
addpops <- as.data.frame(read.csv("addpopcode.csv", header=TRUE),stringsAsFactors = FALSE)
#If you come up with more things to add, you can do this again, accept it will create another X collumn so reaname or remove that one
#And change the bounds on iris2[#,#]

tmp <- cbind.data.frame(tvcf,addpops)
iris4[150,7443]<-"similis"

#reorder
iris2 <- tmp %>%
  select(X,POP,SITE,SPECIES,everything())
#iris2 <- iris4[,-7442] %>%  select(X,POP,SITE,SPECIES,everything())

pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res, data = iris, colour = 'Species')

tmpdf <- iris2[5:7442]
#Okay, this is gross but for now, replace NA with 0 #huh - it basically doesn't change at all if I set NAs to 1 or 2 either cause I guess theres so much other over riding data
tmpdf[is.na(tmpdf)] <- 0
#generated errors that made it turn back into NA?
#FIXED: Always turn into dataframes with with ,stringsAsFactors = FALSE
#replaced 0s have different spacing, no not numerical or something
tmpdf<- as.data.frame(tmpdf,stringsAsFactors = FALSE)
tmpdf1 <- mutate_all(tmpdf, function(x) as.numeric(as.character(x)))
#loses the sample ID but I still have that in a different column for later

pca_res <- prcomp(tmpdf1 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_res, data = iris2, colour = 'POP', shape='SPECIES') #+
  #scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07","#A4A4A4","#E7B803", "#2E1FDF", "#FC4E08"))

#Messing with the population code coloring
theme_set(
  theme_classic() 
)

#Reorder the factor levels:
  # Default order
  levels(iris2$POP)
## [1] "CCB" "GA"  "GBE" "NLI" "SCC" "SJF" "SLI"
# Reverse the order as follow
# iris$Species <- factor(iris$Species, levels = rev(levels(iris$Species)))
# Or specify the factor levels in the order you want
iris2$POP <- factor(iris2$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
iris2$SPECIES <- factor(iris2$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
autoplot(pca_res, data = iris2, colour = 'POP', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  scale_color_viridis_d(option="plasma","Population")  #setting pallete to unknown pallete "v" makes a nice green when done with scale_color_brewer #direction = -1 to flip the color order

#Next I could try adding in the circles that incase populations but also they'd be slightly messy unless I could tell it to do them by site code or by species.

autoplot(pca_res, data = iris2, colour = 'SITE', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(1, 5, 3),"Species")

autoplot(pca_res, data = iris2, colour = 'SITE', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species")


#just similis
justsim <- iris4[,-7442] %>%
  select(X,POP,SITE,SPECIES,everything()) %>%
  filter(SPECIES == "similis")

tmpdfsim <- justsim[5:7442]
tmpdfsim[is.na(tmpdfsim)] <- 0
tmpdfsim<- as.data.frame(tmpdfsim,stringsAsFactors = FALSE)
tmpdfsim1 <- mutate_all(tmpdfsim, function(x) as.numeric(as.character(x)))

pca_ressim <- prcomp(tmpdfsim1)#, center = TRUE, scale = TRUE, na.action = na.omit)
#cannot scale because some of the collumns are all 0 because the similis all match eachother but it should not need to scale because they are all on the same scale
autoplot(pca_ressim, data = justsim, colour = 'POP') #+
  scale_color_viridis_d(option="viridis","Population")#+

autoplot(pca_ressim, data = justsim, colour = 'POP', shape=FALSE, label.size = 1) #+
  scale_color_viridis_d(option="viridis","Population")#+

#the sample labeled as 120 is the NLI sample in SCC territory
#RP20_009_S30





#Misc

#retry adegetnet with new files


datapop <- read.genepop('dDocenthelp/vcfs_929/dedup_929b.gen', ncode=3, quiet = FALSE)

#from
#https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html

library('vcfR')
vcf <- read.vcfR("dDocenthelp/vcfs_929/undedup_nmq30.recode.vcf")
x <- vcfR2genlight(vcf)
#More than two alleles found in some
#removed those ~200
dvcf <- read.vcfR("dDocenthelp/vcfs_929/dedup_nmq30.recode.vcf")
dx <- vcfR2genlight(dvcf)
#removed ~300
library(adegenet)

#give each sample in list a population
pop(x) <- as.factor(c("SFJ", "SFJ", "CCB", "SLI", "SLI", "SLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GBE", "GBE", "GBE", "GBE", "GBE", "NLI", "NLI", "NLI", "NLI", "NLI", "NLI", "NLI", "NLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "NLI", "NLI", "NLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "CCB", "NLI", "NLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC"))
popNames(x)
pop(dx) <- as.factor(c("SFJ", "SFJ", "CCB", "SLI", "SLI", "SLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GBE", "GBE", "GBE", "GBE", "GBE", "NLI", "NLI", "NLI", "NLI", "NLI", "NLI", "NLI", "NLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "NLI", "NLI", "NLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "CCB", "NLI", "NLI", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC"))

ploidy(x) <- 2
ploidy(dx) <- 2


x.dist <- dist(x)
x.dist <- poppr::bitwise.dist(x)

install.packages("poppr")
#?
glPca(x)

#now 
#https://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
#glPca (adegenet): implements PCA for genome-wide SNP data stored as genlight objects; see dedicated tutorial (genomics).
pca1 <- glPca(x, nf=2)
## plot eigenvalues
barplot(pca1$eig, main="eigenvalues", col=heat.colors(length(pca1$eig)))
## basic plot
scatter(pca1, ratio=.2)
## plot showing groups
s.class(pca1$scores, pop(x))
title("Undedup from ~1300 loci with no missing data")


pca2 <- glPca(dx, nf=2)
scatter(pca2, ratio=.2)
s.class(pca2$scores, pop(x))
title("Dedup from ~1600 loci with no missing data")





pca1 <- glPca(x, nf=2)
s.class(pca1$scores, pop(x))
title("Undedup from ~1300 loci with no missing data")


pca2 <- glPca(dx, nf=2)
s.class(pca2$scores, pop(x))
title("Dedup from ~1600 loci with no missing data")

barplot(pca1$eig)

#Then: https://www.rdocumentation.org/packages/adegenet/versions/2.0.1/topics/glPca

## simulate a toy dataset
x <- glSim(50,4e3, 50, ploidy=2)
x
plot(x)

## perform PCA
pca1 <- glPca(x, nf=2)

## plot eigenvalues
barplot(pca1$eig, main="eigenvalues", col=heat.colors(length(pca1$eig)))

## basic plot
scatter(pca1, ratio=.2)

## plot showing groups
s.class(pca1$scores, pop(x), col=colors()[c(131,134)])
title("Undedup from ~1300 loci with no missing data")
add.scatter.eig(pca1$eig,2,1,2)

s.class(pca1$scores, pop(x), col=colors()[c(131,134)])
title("Undedup from ~1300 loci with no missing data")
add.scatter.eig(pca1$eig,2,1,2)




install.packages("viridis")
library(viridis)
library(tidyverse)
library(janitor)
library(dplyr)



library(adegenet)
sim_str  <- read.structure("populations.STR")

sim_test  <- read.structure("populations.stru")
#If I can figure out the answers to all those questions, this could work

sim_test2  <- import2genind("populations.snps.gen")
sim_test3  <- import2genind("populations.haps.gen")

sim_snps <- read.genepop("populations.snps.gen")
sim_haps <- read.genepop("populations.haps.gen")

#It is convertign to agenind object, why not leave it as a genepop

#Lots of NAs?
#straight from the paper
pop(sim_snps) = sapply(strsplit(indNames(sim_snps), '_'), function(sim_snps){sim_snps[1]})
sim_snps.pca = dudi.pca(sim_snps, scannf=FALSE)
s.class(sim_snps.pca$li, pop(sim_snps), col=rainbow(nPop(sim_snps)))
add.scatter.eig(sim_snps.pca$eig[1:10], xax=1, yax=2)

x <- read.genepop("populations.snps.gen")
pop(x) = sapply(strsplit(indNames(x), '_'), function(x){x[1]})
x.pca = dudi.pca(x, scannf=FALSE)
s.class(x.pca$li, pop(x), col=rainbow(nPop(x)))
add.scatter.eig(x.pca$eig[1:10], xax=1, yax=2)


#####  reading from adegenet manual  #####

#To continue with the toy example, we can perform a simple PCA. Allele presence absence
#data are extracted and NAs replaced using tab:
x <- tab(sim_snps, NA.method="mean")  #maybe there are other NA methods. Could I set r = 1.0

## make PCA
pca1 <- dudi.pca(x,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(sim_snps))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]
## basic plot
plot(pca1$li, col=myCol, cex=3, pch=myPch)
## use wordcloud for non-overlapping labels
library(wordcloud)
## Loading required package: RColorBrewer
textplot(pca1$li[,1], pca1$li[,2], words=rownames(x), cex=1.4, new=FALSE) #I have way too many points for this
## legend the axes by adding loadings
abline(h=0,v=0,col="grey",lty=2)
s.arrow(pca1$c1*.5, add.plot=TRUE)
legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7),
       leg=c("Group A","Group B"), pt.cex=2)

#Why didn't it name the axes? How do I know the %variation
#Looking lower in the doc, this is not producing the look of the pca plot I want so find that in the ade manual


tab(obj)
#might let me look at it
#is this snp data? Yes.
#So then. Nevermind, that was for importing


### manipulating
#data(sim_snps)

toto <- genind2genpop(sim_snps)
popNames(toto)  #OOOHH NO THAT NOT RIGHT!!! It's giving sample IDs as pop names
#Altough I do think its probably just the first individual per pop and there's still seven pops, just by the wrong names
nInd(toto)  #And it won't give anything here

#Plot summary data
#TAKES A LONG TIME TO RUN, RETRY AT SOME POINT
titi <- summary(toto) #assume in this case that nancycats is a pop, but it might be ind
names(titi)
## [1] "n" "n.by.pop" "loc.n.all" "pop.n.all" "NA.perc" "Hobs"
## [7] "Hexp"
par(mfrow=c(2,2))
plot(titi$n.by.pop, titi$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(titi$n.by.pop,titi$pop.n.all,lab=names(titi$n.by.pop))
barplot(titi$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(titi$Hexp-titi$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(titi$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)



#PCA from manual page 52, section 6.2
#from geneind
sum(!is.na(sim_snps$tab))

#There are 22,675,576 not missing data.... What the crap? r = 0.8 my butt
#There are 60,179,152 missing data, which will be replaced by tab:
X <- tab(sim_snps, freq = TRUE, NA.method = "mean")
class(X)
dim(X)

pca2 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
#you can plot some stats from making the eigen vales in pca
s.label(pca2$li)
title("PCA of similis SNP datas\nAxes 1-2")
#add.scatter.eig(pca1$eig[1:20], 3,1,2)   #add eigen values in the bottom corner

#instead plot with 
s.class(pca2$li, pop(sim_snps))
title("PCA of similis SNP datas\nAxes 1-2")


colorplot(pca2$li, pca2$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
abline(v=0,h=0,col="grey", lty=2)
title("PCA of microbov dataset\naxes 1-2")

col <- funky(15)
s.class(pca2$li, pop(sim_snps), xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

s.class(pca2$li, pop(sim_snps), col=transp(col,.6), xlab="PC 1", ylab="PC 2")


#This tool is not great, what can I use instead?


#### Now with non dedup


usimplifiedvcf <- read.table("simvcfsimplereadyforR_nocolon.txt", header=TRUE)
#https://github.com/nt246/NTRES6940-data-science/blob/f43638d0616d8f39c3fda6588a85b346a2b86a61/analysis/lesson10-tidy-data.md
#first mutate a collumn that is c___CHROM_p___POS so each loci has a unique identifier
#then make those the collumn names
usvcf <- as.data.frame(usimplifiedvcf)
tmp <- usvcf %>%
  select(CHROM,POS) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":"))
ulocusvcf <- usvcf %>%
  select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  select(LOCUS,everything())
#this is the same as the ID but since I was too lazy to fix the IDs in excel

#filter for loci with % NA before transposing
count_na <- function(x) sum(is.na(x))    

df1 <- ulocusvcf %>%
  mutate(count_na = apply(., 1, count_na))

less <- df1 %>%
  filter(count_na < 10)
#Has a lot more missing data that when I removed the duplicates
#It really feels like pop -r 0.8 aint doing jack shit

#set the headers when I transpose
utvcf = setNames(data.frame(t(less[,-1]),stringsAsFactors = FALSE), less[,1])
#utvcf = setNames(data.frame(t(ulocusvcf[,-1]),stringsAsFactors = FALSE), ulocusvcf[,1])
#utvcf <- utvcf[-c(1,2,3,4,5,6,7,8,9), ]
#remove old collumns
utvcf <- utvcf[-c(1,2,3,4,5,6,7,8,9,160), ]
#add column with population

#tmp <- utvcf %>%
#  select("107:10")
#write.csv(tmp,"addpopcode.csv", row.names = TRUE)
#add in POP, SITE and SPECIES
addpops <- as.data.frame(read.csv("addpopcode.csv", header=TRUE),stringsAsFactors = FALSE)
#If you come up with more things to add, you can do this again, accept it will create another X collumn so reaname or remove that one
#And change the bounds on iris2[#,#]

tmp <- cbind.data.frame(utvcf,addpops)
#iris4[150,7443]<-"similis"

#reorder
uiris2 <- tmp %>%
  select(X,POP,SITE,SPECIES,everything())
#iris2 <- iris4[,-7442] %>%  select(X,POP,SITE,SPECIES,everything())

tmpdf <- uiris2[5:727]
#Okay, this is gross but for now, replace NA with 0 #huh - it basically doesn't change at all if I set NAs to 1 or 2 either cause I guess theres so much other over riding data
tmpdf[is.na(tmpdf)] <- 1
#generated errors that made it turn back into NA?
#FIXED: Always turn into dataframes with with ,stringsAsFactors = FALSE
#replaced 0s have different spacing, no not numerical or something
#tmpdf<- as.data.frame(tmpdf,stringsAsFactors = FALSE)
tmpdf1 <- mutate_all(tmpdf, function(x) as.numeric(as.character(x)))
#loses the sample ID but I still have that in a different column for later

pca_res <- prcomp(tmpdf1 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_res, data = uiris2, colour = 'POP', shape='SPECIES') #+
#scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07","#A4A4A4","#E7B803", "#2E1FDF", "#FC4E08"))

#Messing with the population code coloring
theme_set(
  theme_classic() 
)

#Reorder the factor levels:
# Default order
levels(iris2$POP)
## [1] "CCB" "GA"  "GBE" "NLI" "SCC" "SJF" "SLI"
# Reverse the order as follow
# iris$Species <- factor(iris$Species, levels = rev(levels(iris$Species)))
# Or specify the factor levels in the order you want
uiris2$POP <- factor(iris2$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
uiris2$SPECIES <- factor(iris2$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
autoplot(pca_res, data = uiris2, colour = 'POP', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  scale_color_viridis_d(option="plasma","Population")  #setting pallete to unknown pallete "v" makes a nice green when done with scale_color_brewer #direction = -1 to flip the color order

#Next I could try adding in the circles that incase populations but also they'd be slightly messy unless I could tell it to do them by site code or by species.

autoplot(pca_res, data = iris2, colour = 'SITE', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(1, 5, 3),"Species")

autoplot(pca_res, data = iris2, colour = 'SITE', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species")


autoplot(pca_res, data = iris2, colour = 'SITE',size=2, shape=FALSE, label.size = 1) 

#just similis
ujustsim <- uiris2 %>%
  select(X,POP,SITE,SPECIES,everything()) %>%
  filter(SPECIES == "similis")

tmpdfsim <- ujustsim[5:727]
tmpdfsim[is.na(tmpdfsim)] <- 0
tmpdfsim<- as.data.frame(tmpdfsim,stringsAsFactors = FALSE)
tmpdfsim1 <- mutate_all(tmpdfsim, function(x) as.numeric(as.character(x)))

pca_ressim <- prcomp(tmpdfsim1)#, center = TRUE, scale = TRUE, na.action = na.omit)
#cannot scale because some of the collumns are all 0 because the similis all match eachother but it should not need to scale because they are all on the same scale
autoplot(pca_ressim, data = justsim, colour = 'POP') #+
scale_color_viridis_d(option="viridis","Population")#+

autoplot(pca_ressim, data = justsim, colour = 'POP', shape=FALSE, label.size = 1) #+
scale_color_viridis_d(option="viridis","Population")#+





# dDocent raws in normal method #
usimplifiedvcf <- read.table("vcfwithnanforR_raw_a_form.txt", header=TRUE)
usvcf <- as.data.frame(usimplifiedvcf)
tmp <- usvcf %>%
  select(CHROM,POS) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":"))
ulocusvcf <- usvcf %>%
  select(CHROM,POS,everything()) %>%
  mutate(LOCUS = paste(CHROM,POS,sep=":")) %>%
  select(LOCUS,everything())
#this is the same as the ID but since I was too lazy to fix the IDs in excel
#filter for loci with % NA before transposing
count_na <- function(x) sum(is.na(x))    
df1 <- ulocusvcf %>%
  mutate(count_na = apply(., 1, count_na))
less <- df1 %>%
  filter(count_na < 10)
#3700 snps for ddoc

#set the headers when I transpose
utvcf = setNames(data.frame(t(less[,-1]),stringsAsFactors = FALSE), less[,1])
#utvcf = setNames(data.frame(t(ulocusvcf[,-1]),stringsAsFactors = FALSE), ulocusvcf[,1])
#utvcf <- utvcf[-c(1,2,3,4,5,6,7,8,9), ]
#remove old collumns
utvcf <- utvcf[-c(1,2,3,4,5,6,7,8,9,160), ]
#add column with population

#tmp <- utvcf %>%
#  select("107:10")
#write.csv(tmp,"addpopcode.csv", row.names = TRUE)
#add in POP, SITE and SPECIES
addpops <- as.data.frame(read.csv("addpopcode.csv", header=TRUE),stringsAsFactors = FALSE)
#If you come up with more things to add, you can do this again, accept it will create another X collumn so reaname or remove that one
#And change the bounds on iris2[#,#]

tmp <- cbind.data.frame(utvcf,addpops)
#iris4[150,7443]<-"similis"

#reorder
uiris2 <- tmp %>%
  select(X,POP,SITE,SPECIES,everything())
#iris2 <- iris4[,-7442] %>%  select(X,POP,SITE,SPECIES,everything())

tmpdf <- uiris2[5:3734]
#Okay, this is gross but for now, replace NA with 0 #huh - it basically doesn't change at all if I set NAs to 1 or 2 either cause I guess theres so much other over riding data
tmpdf[is.na(tmpdf)] <- 1 #call NA 1 for quickness to make pca when full collmns cannot be 0
#generated errors that made it turn back into NA?
#FIXED: Always turn into dataframes with with ,stringsAsFactors = FALSE
#replaced 0s have different spacing, no not numerical or something
#tmpdf<- as.data.frame(tmpdf,stringsAsFactors = FALSE)
tmpdf1 <- mutate_all(tmpdf, function(x) as.numeric(as.character(x)))
#loses the sample ID but I still have that in a different column for later

pca_res <- prcomp(tmpdf1 , scale = TRUE)#, center = TRUE, scale = TRUE, na.action = na.omit)
autoplot(pca_res, data = uiris2, colour = 'POP', shape='SPECIES') #+
#scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07","#A4A4A4","#E7B803", "#2E1FDF", "#FC4E08"))

#Messing with the population code coloring
theme_set(
  theme_classic() 
)

#Reorder the factor levels:
# Default order
levels(iris2$POP)
## [1] "CCB" "GA"  "GBE" "NLI" "SCC" "SJF" "SLI"
# Reverse the order as follow
# iris$Species <- factor(iris$Species, levels = rev(levels(iris$Species)))
# Or specify the factor levels in the order you want
uiris2$POP <- factor(iris2$POP, levels = c("GBE","CCB","SCC","NLI","SLI","SJF","GA"))
uiris2$SPECIES <- factor(iris2$SPECIES, levels = c("similis","solidissima","possible_heterozygote"))
autoplot(pca_res, data = uiris2, colour = 'POP', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species") +
  scale_color_viridis_d(option="plasma","Population")  #setting pallete to unknown pallete "v" makes a nice green when done with scale_color_brewer #direction = -1 to flip the color order

#Next I could try adding in the circles that incase populations but also they'd be slightly messy unless I could tell it to do them by site code or by species.

autoplot(pca_res, data = iris2, colour = 'SITE', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(1, 5, 3),"Species")

autoplot(pca_res, data = iris2, colour = 'SITE', shape='SPECIES',size=2) +
  scale_shape_manual(values=c(16, 18, 3),"Species")


autoplot(pca_res, data = iris2, colour = 'SITE',size=2, shape=FALSE, label.size = 1) 

#just similis
ujustsim <- uiris2 %>%
  select(X,POP,SITE,SPECIES,everything()) %>%
  filter(SPECIES == "similis")

tmpdfsim <- ujustsim[5:3734]
tmpdfsim[is.na(tmpdfsim)] <- 0
tmpdfsim<- as.data.frame(tmpdfsim,stringsAsFactors = FALSE)
tmpdfsim1 <- mutate_all(tmpdfsim, function(x) as.numeric(as.character(x)))

pca_ressim <- prcomp(tmpdfsim1)#, center = TRUE, scale = TRUE, na.action = na.omit)
#cannot scale because some of the collumns are all 0 because the similis all match eachother but it should not need to scale because they are all on the same scale
autoplot(pca_ressim, data = justsim, colour = 'POP') #+
scale_color_viridis_d(option="viridis","Population")#+

autoplot(pca_ressim, data = justsim, colour = 'POP', shape=FALSE, label.size = 1) #+
scale_color_viridis_d(option="viridis","Population")#+
