# Naming Convension for primers from contigs
Trinity_out <- read.csv("~/Documents/Hare_Lab/Spisula/LongestIso_Sol_naming.csv",stringsAsFactors=F, header=FALSE)

Contigs <- data.frame(Trinity_out)
colnames(Contigs) <- c("Gene","Isoform","Length")
library(stringr)
library(dplyr)
library(tidyr)

Contigs$DN <- NA
Contigs$c <- NA
Contigs$g <- NA
Contigs$i <- NA
Contigs$Primer <- NA
Contigs$Qprimer <-NA

#str_split_fixed(Contigs$Isoform[i], "_", 5)

#Contigs %>%
#  separate(Isoform,c("DN","c","g","i"),"_")


for (i in 1:nrow(Contigs)) {
  Contigs$DN[i] <- sub("DN","",str_split_fixed(Contigs$Isoform[i], "_", 5)[,2])
  Contigs$c[i] <- sub("c","",str_split_fixed(Contigs$Isoform[i], "_", 5)[,3])
  Contigs$g[i] <- sub("g","",str_split_fixed(Contigs$Isoform[i], "_", 5)[,4])
  Contigs$i[i] <- sub("i","",str_split_fixed(Contigs$Isoform[i], "_", 5)[,5])
  print(i)
  Contigs$Primer[i]<-paste(Contigs$DN[i], "_", Contigs$c[i], Contigs$g[i], sep = "")
}

Contigs$Qprimer <- paste(Contigs$DN, Contigs$c, Contigs$g, sep = "")
Contigs$Qprimer <- paste(Contigs$DN, Contigs$c, Contigs$g, sep = "")
#Test for shorter names
duplicated(Contigs$Qprimer)
sum(duplicated(Contigs$Primer))
#Nothing happens to overlap

max(nchar(Contigs$Qprimer))
9
max(nchar(Contigs$Primer))
10

cg_set <- Contigs %>%
  select(5,6) %>%
  unique()

cg_set<-c(1,2,21,41,42,51,11,22,31,61,91,71,81,10,32,12,13,)
#nope doesn't work to sorten it

#Try c vs g count
#c0,g1 -> 1
#c0,g2 -> 2

#c30,g1


summary(nchar(Contigs$Qprimer))
#Min.   1st Qu.  Median    Mean     3rd Qu.    Max. 
#3.000  6.000    7.000     6.755    7.000      9.000 

#Spis#########R
#Max primer name length, not too bad


#Then check for no repeat isoforms
#I suppose by looking for duplicates
duplicated(Contigs$Primer)
sum(duplicated(Contigs$Primer))
#0 no duplicates

write.csv(Contigs,"~/Documents/Hare_Lab/Spisula/Trinity_Sol_Primer_Names.csv", row.names = FALSE)


#What are the most common enzymes?
# https://convert.town/column-to-comma-separated-list , hit setting to auto add ""
Enzymes <- c("BsoBI","MspI","BstNI","Smal","Xbal","HpaII","AvaI","DraI","HinfI","BspDI","EcoRI","BcII","BcII","Ncil","BstNI","BstNI","CviQI","MspI","AseI","MspI","Hpall","Dral","HinfI","Ndel","Hinfl","Xbal","BtsCI","BstEII","BstEII","BstYI","CviQI","DraI","HinP1I","Xbal","Hhal","BtsCI","HinfI","BtsCI","BtsCI","AvaII","HinfI","BsaHI","BtsCI")

for (E in unique(Enzymes)){
  print(paste(E," ",sum(Enzymes==E)))
}
#BtsCI  5 no: most expensive
#HinfI  4 got it
#Xbal   3 got it
#MspI   3 skipped over cause whereever it cuts, one other also generally cuts
#BstNI  3 got it
