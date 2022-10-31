###### SpaceMix Runs ######
#Date: May 6th

vignette("spacemix_vignette")
#### Install ####
require(devtools)
#install_github("gbradburd/SpaceMix",build_vignettes=TRUE)
library(spam)
require(SpaceMix)

#View tutorial
vignette("spacemix_vignette")

#You can see all documented functions with
??SpaceMix

#Set workdir
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/")
getwd()

#### Input File Creation ####
#Make yout own allele counts matrix instructions
  #Start with filtered vcf (1fp or LD pruned)
  #Remove Header and Transpose
  #Convert vcf to 0/0 format (using left 3)
  #Sum the counts of 1 for each allele for each population
    #SUMPRODUCT(LEN(Sheet3!B10:B19)-LEN(SUBSTITUTE(Sheet3!B10:B19,"1","")))+SUMPRODUCT(LEN(Sheet3!B110:B115)-LEN(SUBSTITUTE(Sheet3!B110:B115,"1","")))
    #Ex:
      #NLI    3 2 1 0 0 2 
      #CTM    2 4 0 1 0 3
  
  #Create a file of the number of alleles with no missing data per locus
    #=32 - SUMPRODUCT(LEN(Sheet3!B10:B19)-LEN(SUBSTITUTE(Sheet3!B10:B19,".","")))+SUMPRODUCT(LEN(Sheet3!B110:B115)-LEN(SUBSTITUTE(Sheet3!B110:B115,".",""))) #Where 32 is 2 x N individuals at that site
    #Ex: (except remove names of locations)
      #NLI    26 26 26 26 23 26 
      #CTM    46 46 45 46 45 46

  #Make Location File (manually)
    #Ex: (except remove names of locations)
      #NLI    -70 30
      #CTM    -72 40


#### Similis Input ####
#Import
allelecounts1 <- read.csv("./Input/allelecounts1_sim.csv", header=FALSE)
samplecounts1 <- read.csv("./Input/samplecount1_sim_REAL.csv", header=FALSE)
locations1 <- read.csv("./Input/poplocations1_sim.csv", header=FALSE)

#Print dimensions
dim(allelecounts1)
dim(samplecounts1)
dim(locations1)

#Convert from list to the same type as the original
#Input dimensions from above
allelecounts1 <- matrix(unlist(allelecounts1), ncol =1508, nrow =6)
samplecounts1 <- matrix(unlist(samplecounts1), ncol =1508, nrow =6)
locations1 <- matrix(unlist(locations1), ncol =2, nrow =6)

#### Similis Run ####
#Options
# 1. "no_movement" 
# 2. "source"
# 3. "target"
# 4. "source_and_target"

#Run
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0522/Output")
run.spacemix.analysis(n.fast.reps = 10,
                      fast.MCMC.ngen = 1e5,
                      fast.model.option = "target",
                      long.model.option = "target",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = allelecounts1,
                      sample.sizes = samplecounts1,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = locations1[,1],
                      spatial.prior.Y.coordinates = locations1[,2],
                      round.earth = FALSE,
                      long.run.initial.parameters=NULL,
                      k = nrow(allelecounts1),
                      loci = ncol(samplecounts1),
                      ngen = 1e6,
                      printfreq = 1e3,
                      samplefreq = 1e3,
                      mixing.diagn.freq = 50,
                      savefreq = 1e5,
                      directory=NULL,
                      prefix = "test1")

#### Understand Output ####
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")

#Graph
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
population.coordinates <- locations1
sample.names <- unlist(c("NLI","GA","CTM","WFH","ELPPB","MW"))
sample.colors <- rainbow(n=6,start=0.1,end=0.8)[as.numeric(cut(population.coordinates[,1],6))]
#repalce 6 above with the number of populations
sim1_source.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_204492/simtest_run1_sorce_LongRun/simtest_run1_sorce_space_MCMC_output1.Robj",
                                                    geographic.locations = population.coordinates,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)
sim1_st.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_753303/simtest_run1_st_LongRun/simtest_run1_st_space_MCMC_output1.Robj",
                                                        geographic.locations = population.coordinates,
                                                        name.vector = sample.names,
                                                        color.vector = sample.colors,
                                                        quantile=0.95,
                                                        burnin=0)
sim1_target.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_445102/simtest_run1_target2_LongRun/simtest_run1_target2_space_MCMC_output1.Robj",
                                                    geographic.locations = population.coordinates,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)
sim1_nm.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_262311/simtest_run1_nm_LongRun/simtest_run1_nm_space_MCMC_output1.Robj",
                                                        geographic.locations = population.coordinates,
                                                        name.vector = sample.names,
                                                        color.vector = sample.colors,
                                                        quantile=0.95,
                                                        burnin=0)
#Similis Highlight 



#Graph the differences between runs
#Target Only
make.spacemix.map(spacemix.map.list = sol2_target.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
make.spacemix.map(spacemix.map.list = sim1_target.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#Source
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
make.spacemix.map(spacemix.map.list = sim1_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#Source and Target
make.spacemix.map(spacemix.map.list = sol2_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
make.spacemix.map(spacemix.map.list = sim1_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#No Movement or Mixture
make.spacemix.map(spacemix.map.list = sol2_nm.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=FALSE,xlim=c(-77,-65),ylim=c(36,44))
make.spacemix.map(spacemix.map.list = sim1_nm.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=FALSE,xlim=c(-82,-68),ylim=c(30,43))

#### Plot Maps of Equal Size ####
library(sf)
library(raster)
library(dplyr)
library(spData)
#library(spDataLarge)
library(tmap)
#library(leaflet) # for interactive maps
library(ggplot2)
library(tidyverse)
#Show full map
tm_shape(us_states) + tm_polygons() + tm_graticules() 
#create my bounds with st_bbox, piped into st as sfc
#eastcoast=st_bbox(c(xmin = -85, xmax = -65,
#                    ymin = 30, ymax = 45)) %>% st_as_sfc()
solterritory=st_bbox(c(xmin = -77, xmax = -65,
                    ymin = 35, ymax = 45)) %>% st_as_sfc()
simterritory=st_bbox(c(xmin = -83, xmax = -67,
                    ymin = 29, ymax = 43)) %>% st_as_sfc()
#color and formatting
#tm_shape(us_states, bbox = eastcoast) +  
#  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
#  tm_graticules(col="white")
tm_shape(us_states, bbox = solterritory) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey99")
bary_iso <- st_read("/Users/hannah/gitHub/Spisula/Work/shape/bathymetry_l_v2.shp")

tm_shape(us_states, bbox = solterritory) +  
  tm_fill(col = "grey99") + tm_borders(col="grey10") + tm_layout(bg.color = "grey99") + 
  tm_polygons() + tm_shape(bary_us) + tm_lines() #+ tm_graticules()
tm_shape(us_states, bbox = simterritory) +  
  tm_fill(col = "grey99") + tm_borders(col="grey10") + tm_layout(bg.color = "grey99") + 
  tm_polygons() + tm_shape(bary_us) + tm_lines() #+ tm_graticules()


#### Testing Graphing #####
#Multiple Graphing Options
make.spacemix.map(spacemix.map.list = example.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=FALSE)
#Add location estimates
make.spacemix.map(spacemix.map.list = example.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#Tweak axes & Visualize addmixture (how is the different from the above?)
make.spacemix.map(example.spacemix.map.list,
                  source.option=TRUE,
                  text=TRUE,xlim=c(-82,-68),ylim=c(30,46))
#Just Location Estimates
make.spacemix.map(example.spacemix.map.list,
                  source.option=FALSE,
                  text=FALSE,xlim=c(-82,-68),ylim=c(30,46))

make.spacemix.map(spacemix.map.list = example.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)

#Try to Add Arrows
#Perhaps they do by default so long as there is actually admixture (also about the coloration of the circles)
#Or try the highlighting function
make.spacemix.map(example.spacemix.map.list,
                  source.option=TRUE,
                  text=TRUE,xlim=c(-82,-68),ylim=c(30,46))
query.spacemix.map(unlist(c("NLI","CTM","WFH","ELPPB","MW", "GA")), #Puts them in WEIRD SPOTS
                   example.spacemix.map.list, 
                   ellipses =TRUE,
                   source.option = TRUE)
make.spacemix.map(example.spacemix.map.list,
                  ellipses = FALSE,
                  source.option=FALSE,
                  text=TRUE,xlim=c(-82,-68),ylim=c(30,46))
query.spacemix.map(unlist(c("GA")), #Puts them in WEIRD SPOTS
                   example.spacemix.map.list, 
                   ellipses =TRUE,
                   source.option = TRUE)

#Doing these on top of eachother are graphic on top of eachother
#Having many of top of eachother shows that there IS some color in my ovals, just very faint
make.spacemix.map(example.spacemix.map.list,
                  source.option=TRUE,
                  ellipses =FALSE,
                  text=TRUE,xlim=c(-82,-68),ylim=c(30,46))

#Pull out likelihood


#### Solidissima Input ####
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/")

#Import
allelecounts2 <- read.csv("./Input/allelecounts2_sol.csv", header=FALSE)
samplecounts2 <- read.csv("./Input/samplecount2_sol_REAL.csv", header=FALSE)
locations2 <- read.csv("./Input/poplocations2_sol.csv", header=FALSE)

#Print dimensions
dim(allelecounts2)
dim(samplecounts2)
dim(locations2)

#Convert from list to the same type as the original
#Input dimensions from above
allelecounts2 <- matrix(unlist(allelecounts2), ncol =2794, nrow =19)
samplecounts2 <- matrix(unlist(samplecounts2), ncol =2794, nrow =19)
locations2 <- matrix(unlist(locations2), ncol =2, nrow =19)

#### Solidissima Run ####
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
run.spacemix.analysis(n.fast.reps = 10,
                      fast.MCMC.ngen = 1e5,
                      fast.model.option = "no_movement",
                      long.model.option = "no_movement",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = allelecounts2,
                      sample.sizes = samplecounts2,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = locations2[,1],
                      spatial.prior.Y.coordinates = locations2[,2],
                      round.earth = FALSE,
                      long.run.initial.parameters=NULL,
                      k = nrow(allelecounts2),
                      loci = ncol(samplecounts2),
                      ngen = 1e6,
                      printfreq = 1e2,
                      samplefreq = 1e3,
                      mixing.diagn.freq = 50,
                      savefreq = 1e5,
                      directory=NULL,
                      prefix = "soltest_run2_nm")

#### Plot Solidissima ####
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
population.coordinates <- locations2
sample.names <- unlist(c('X413','X536','BAR','BLP','BRW','CGC','FNJ','GBE','MCX','MW','NAT','PLY','PT','PVT','SHN','SsLI1k','SsLI23k','SsLI4k','WV'))
sample.colors <- rainbow(n=19,start=0.1,end=0.8)[as.numeric(cut(population.coordinates[,1],19))]
#repalce 6 above with the number of populations

#Source / Target
sol2_st.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_943071/soltest_run2_st_LongRun/soltest_run2_st_space_MCMC_output1.Robj",
                                                    geographic.locations = population.coordinates,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)
#Just #3 Target
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
sol2_target.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_756859/soltest_run2_target_LongRun/soltest_run2_target_space_MCMC_output1.Robj",
                                                    geographic.locations = population.coordinates,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)
#Source...? What are the differences?
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
sol2_source.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_909149/soltest_run2_source_LongRun/soltest_run2_source_space_MCMC_output1.Robj",
                                                        geographic.locations = population.coordinates,
                                                        name.vector = sample.names,
                                                        color.vector = sample.colors,
                                                        quantile=0.95,
                                                        burnin=0)
#NM
sol2_nm.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_275884/soltest_run2_nm_LongRun/soltest_run2_nm_space_MCMC_output1.Robj",
                                                    geographic.locations = population.coordinates,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)


setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#HEY, the locations are CORRECT!!! Nice, I think that target vs source are switched in their descriptions
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)
#First time I've seen arrows! 
#Try highlighting
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=TRUE)
#Plots just arrows, no ellipses!! BELOW
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)
#Trying plotting that with highlighting
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)
query.spacemix.map("NAT",
                   sol2_source.spacemix.map.list,
                   source.option = TRUE)
#Trying plotting that with highlighting with PLY
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)
query.spacemix.map("PLY",
                   sol2_source.spacemix.map.list,
                   source.option = TRUE,
                   ellipses=TRUE) #Adds faint line (can see more clearly if plot multiple times on top)
#Try more
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)
query.spacemix.map("GBE",
                   sol2_source.spacemix.map.list,
                   source.option = TRUE,
                   ellipses=TRUE)
#Expand the view to see the elipses better
make.spacemix.map(spacemix.map.list = sol2_source.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE,xlim=c(-78,-64),ylim=c(36,47))
query.spacemix.map("GBE",
                   sol2_source.spacemix.map.list,
                   source.option = TRUE)
#Try the expanded view on the Soruce and Target
make.spacemix.map(spacemix.map.list = sol2_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=TRUE,xlim=c(-78,-64),ylim=c(36,47))
make.spacemix.map(spacemix.map.list = sol2_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=TRUE,xlim=c(-74,-69),ylim=c(39.5,42))
query.spacemix.map("GBE",
                   sol2_st.spacemix.map.list,
                   source.option = TRUE)


#Target Only
make.spacemix.map(spacemix.map.list = sol2_target.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#Still plots wrong.... - not wrong, just estimating "geogenetic space"
make.spacemix.map(sol2_target.spacemix.map.list,
                  source.option=TRUE,
                  text=TRUE)
make.spacemix.map(sol2_st.spacemix.map.list,
                  source.option=TRUE,
                  text=TRUE)


#Source and Target
#Multiple Graphing Options
make.spacemix.map(spacemix.map.list = sol2_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)
#Did NOT plot locations correctly - brings them super close together - Instead try targ - are they being brought super close together by the model?
library(ggplot2)
tmp <- as.data.frame(population.coordinates)
ggplot(data = tmp, aes(x = V1, y = V2)) +
  geom_point()


make.spacemix.map(sol2_st.spacemix.map.list,
                  source.option=TRUE,
                  ellipses =FALSE,
                  text=TRUE,xlim=c(-76,-66),ylim=c(37,43))

make.spacemix.map(spacemix.map.list = sol2_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
#Tweak axes & Visualize addmixture (how is the different from the above?)
make.spacemix.map(sol2_st.spacemix.map.list,
                  source.option=TRUE,
                  text=TRUE,xlim=c(-82,-68),ylim=c(30,46))
#Just Location Estimates
make.spacemix.map(sol2_st.spacemix.map.list,
                  source.option=FALSE,
                  text=FALSE,xlim=c(-82,-68),ylim=c(30,46))

make.spacemix.map(sol2_st.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=TRUE)



#### Solidissima Alternates - MA ####
#Solidissima doesn't look quite right... Try subsetting the data
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/")

### Massachusetts Only
allelecounts2 <- read.csv("./Input/allelecounts2_sol_MA.csv", header=FALSE)
samplecounts2 <- read.csv("./Input/samplecount2_sol_MA.csv", header=FALSE)
locations2 <- read.csv("./Input/poplocations2_sol_MA.csv", header=FALSE)

#Print dimensions
dim(allelecounts2)
dim(samplecounts2)
dim(locations2)

#Convert from list to the same type as the original
#Input dimensions from above
allelecounts2 <- matrix(unlist(allelecounts2), ncol =2794, nrow =9)
samplecounts2 <- matrix(unlist(samplecounts2), ncol =2794, nrow =9)
locations2 <- matrix(unlist(locations2), ncol =2, nrow =9)

setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output")
run.spacemix.analysis(n.fast.reps = 10,
                      fast.MCMC.ngen = 1e5,
                      fast.model.option = "source_target",
                      long.model.option = "source_target",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = allelecounts2,
                      sample.sizes = samplecounts2,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = locations2[,1],
                      spatial.prior.Y.coordinates = locations2[,2],
                      round.earth = FALSE,
                      long.run.initial.parameters=NULL,
                      k = nrow(allelecounts2),
                      loci = ncol(samplecounts2),
                      ngen = 1e6,
                      printfreq = 1e2,
                      samplefreq = 1e3,
                      mixing.diagn.freq = 50,
                      savefreq = 1e5,
                      directory=NULL,
                      prefix = "sol2_MA_st")
#Does not work- fails object 'accept_rates' not found
#Filters to 2781 loci out of the original 2794 left in curated dataset
#Don't know how to fix


#### Subset Solidissima by A vs B ####
#Also rerun sim <- did not work? Negative infinity?
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0522/")
getwd()
require(devtools)
library(spam)
require(SpaceMix)
#Run and plot both
#Import
allelecounts2 <- read.csv("./Input/counts_sim.csv", header=FALSE)
samplecounts2 <- read.csv("./Input/size_sim.csv", header=FALSE)
locations2 <- read.csv("./Input/locations_sim.csv", header=FALSE)
#Options
# 1. "no_movement" 
# 2. "source"
# 3. "target"
# 4. "source_and_target"

#Print dimensions
dim(allelecounts2)
dim(samplecounts2)
dim(locations2)

#Convert from list to the same type as the original
#Input dimensions from above
allelecounts2 <- matrix(unlist(allelecounts2), ncol =1508, nrow =6)
samplecounts2 <- matrix(unlist(samplecounts2), ncol =1508, nrow =6)
locations2 <- matrix(unlist(locations2), ncol =2, nrow =6)

setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0522/Output")
run_name <- "sim_no_movement"
run.spacemix.analysis(n.fast.reps = 10,
                      fast.MCMC.ngen = 1e5,
                      fast.model.option = "no_movement",
                      long.model.option = "no_movement",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = allelecounts2,
                      sample.sizes = samplecounts2,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = locations2[,1],
                      spatial.prior.Y.coordinates = locations2[,2],
                      round.earth = FALSE,
                      long.run.initial.parameters=NULL,
                      k = nrow(allelecounts2),
                      loci = ncol(samplecounts2),
                      ngen = 1e6,
                      printfreq = 1e3,
                      samplefreq = 1e3,
                      mixing.diagn.freq = 50,
                      savefreq = 1e5,
                      directory= NULL,
                      prefix = run_name)

##### New Plots May 22 #####
#I put all the output files in one folder so I don't have to do the run_### thing

### Input ###
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0522/Output/best_runs")
species <- "sim"
model <- "source_and_target"
folder <- paste(species,"_",model,"_LongRun",sep="")
mcmcfile <- paste("./",folder,"/",species,"_",model,"_space_MCMC_output1.Robj",sep="")

locations2 <- read.csv(paste("../../Input/locations_",species,".csv",sep=""), header=FALSE)
locations2 <- matrix(unlist(locations2), ncol =2, nrow =dim(locations2)[1])
population.coordinates <- locations2
sample.colors <- rainbow(n=dim(locations2)[1],start=0.1,end=0.8)[as.numeric(cut(population.coordinates[,1],dim(locations2)[1]))]

#List of names
samples_B <- unlist(c('X413', 'X536', 'BAR', 'BRW', 'CGC', 'FNJ', 'GBE', 'MW', 'NAT', 'PLY', 'PT', 'PVT', 'SLI'))
samples_A <- unlist(c('BLP', 'MCX', 'SHN', 'SsLI1k', 'SsLI23k', 'SsLI4k', 'WV'))
samples_sim <- unlist(c('NLI', 'GA', 'CTM', 'WFH', 'ELPPB', 'MW'))

sample.names <- get(paste("samples_", species, sep="")) #Use get to turn a pasted string variable name into a output from an actual variable

#Set up Plot
mapname <- paste(species, "_",model,"_","mixmap",sep="")
assign(paste(mapname,sep=""),make.spacemix.map.list(MCMC.output.file = mcmcfile, 
        geographic.locations = population.coordinates,
        name.vector = sample.names,
        color.vector = sample.colors,
        quantile=0.95,
        burnin=0))
#Print Plots - for Real saving
make.spacemix.map(spacemix.map.list = get(mapname),
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=TRUE,
                  #xlim=c(-74,-69),ylim=c(39,43)) #Sol B
                  #xlim=c(-74,-71),ylim=c(40,41.5)) #Sol A basic
                  #xlim=c(-74.5,-71),ylim=c(39,42)) #Sol A zoomed out
                  xlim=c(-82,-68),ylim=c(31,45)) #For sim
#query.spacemix.map("GA",
#query.spacemix.map(c('NLI', 'CTM', 'WFH', 'ELPPB', 'MW'),
#query.spacemix.map(c('BAR', 'BRW', 'CGC', 'MW', 'PLY', 'PT', 'PVT', 'SLI'),
#query.spacemix.map(c('BLP', 'MCX', 'SHN', 'SsLI1k', 'SsLI23k', 'SsLI4k', 'WV'),
query.spacemix.map(c('BLP', 'MCX', 'SHN', 'SsLI1k', 'SsLI23k', 'SsLI4k', 'WV'),
                   get(mapname),
                   ellipses = FALSE,
                   source.option = TRUE)
plot.admix.arrows(admix.source.coords=get(mapname)$MAPP.admix.source.coords,
                  geogen.coords=get(mapname)$MAPP.geogen.coords,
                  admix.proportions = get(mapname)[["MCMC.output"]][["admix.proportions"]],
                  colors=NULL,
                  length=0.2)
#Figured out how to do arrows myself by setting lwd
arrows(x0 = get(mapname)$MAPP.admix.source.coords[,1],
       y0 = get(mapname)$MAPP.admix.source.coords[,2],
       x1 = get(mapname)$MAPP.geogen.coords[,1],
       y1 = get(mapname)$MAPP.geogen.coords[,2],
       lwd = get(mapname)[["MCMC.output"]][["admix.proportions"]] * 100, #Increase the visibility
       col = get(mapname)$color.vector,
       length=0.2)

#col = rep("black",nrow(admix.source.coords=get(mapname)$MAPP.admix.source.coords))

#Save Plot
#jpeg(file=paste("./FiguresOutput/",species,"_",model,".jpeg",sep=""))
#make.spacemix.map(spacemix.map.list = get(mapname),
#                  text=TRUE,
#                  ellipses=TRUE,
#                  source.option=FALSE)
#dev.off()

### Plot Likelihood/ Other Outputs ###
#Load other output data?
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0522/Output/best_runs")
species <- "sim"
model <- "source_and_target"
folder <- paste(species,"_",model,"_LongRun",sep="")
ex.output <- load_MCMC_output(paste("./",folder,"/",species,"_",model,"_space_MCMC_output1.Robj",sep=""))
MCN.freq <- load_MCMC_output(paste("./",folder,"/",species,"_",model,"_MCN.frequencies.list.Robj",sep=""))
MCN.frequencies.list <- MCN.freq$MCN.frequencies.list #Go a level deeper

#Calculate
sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                         use="pairwise.complete.obs")
k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/ex.output$last.params$inv.mean.sample.sizes / 
                                (sum(1/ex.output$last.params$inv.mean.sample.sizes)),
                              nrow=k,ncol=k,byrow=TRUE)
MC.parametric.covariance <- (MC.matrix) %*%     
  ex.output$last.params$admixed.covariance %*% 
  t(MC.matrix)
index.matrix <- upper.tri(sample.covariance,diag=TRUE)

#Plot
plot(sample.covariance[index.matrix], 
     MC.parametric.covariance[index.matrix],
     col=adjustcolor("black",0.3),pch=20,
     xlab="Sample covariance",
     ylab="Parametric covariance",
     main=paste("Model adequacy:\n Matrix Comparison\n","Genotype",species,"   ", model))
abline(0,1,col="red")
#IBD
plot(ex.output$last.params$D[1:k,1:k][index.matrix], 
     sample.covariance[index.matrix],
     pch=19,col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main=paste("Model adequacy:\n IBD patterns\n","Genotype",species,"   ", model))
points(ex.output$last.params$D[1:k,1:k][index.matrix], 
       MC.parametric.covariance[index.matrix],col="red",pch=20)
legend(x="topright",pch=19,col=c(1,2),
       legend=c("observed","model estimate"))
#Trace Plot
plot(ex.output$Prob,xlab="MCMC iterations",ylab="value", type='l',
     main=paste("Posterior probability trace plot\n","Genotype",species,"   ", model))
#Joint Marginal
plot(ex.output$a0,ex.output$a1,xlab="a0",ylab="a1",
     main=paste("Joint marginal of a0 and a1\n","Genotype",species,"   ", model),pch=20,
     col=adjustcolor(rainbow(1000,start=4/6,end=6/6),0.3))
legend(x="bottomright",pch=19,cex=0.8,
       col=rainbow(1000,start=4/6,end=6/6)[c(1,500,1000)],
       legend=c("Sampled MCMC iteration 1",
                "Sampled MCMC iteration 500",
                "Sampled MCMC iteration 1000"))
# Acceptance Rates
plot(ex.output$accept_rates$a0_accept_rate,
     xlab="MCMC iterations",ylab="Acceptance rate",
     main=paste("Acceptance rate of a0\n","Genotype",species,"   ", model),type='l',
     ylim=c(0.35,0.6))
abline(h=0.44,col="gray",lty=2)
matplot(t(ex.output$accept_rates$nugget_accept_rate),
        xlab="MCMC iterations",ylab="Acceptance rate",
        main=paste("Acceptance rates of nuggets\n","Genotype",species,"   ", model),type='l',
        ylim=c(0.3,0.7))
abline(h=0.44,col="gray",lty=2)


##### Figure Out How to Load Files #####
#Start with example output data
vignette("spacemix_vignette")
data(ex.output)
# Trace plot of posterior probability over the MCMC
plot(ex.output$Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot",type='l')
# Acceptance rate of a0 over the course of the 
#   MCMC analysis
plot(ex.output$accept_rates$a0_accept_rate,
     xlab="MCMC iterations",ylab="Acceptance rate",
     main="Acceptance rate of a0",type='l',
     ylim=c(0.35,0.6))
abline(h=0.44,col="gray",lty=2)

data(MCN.frequencies.list.RData)
# first, load the standardized (mean-centered and normalized)
#   allele frequency data object.  This object, which is the 
#   "MCN.frequencies.list" (Mean Centered and Normalized) is 
#   saved in the Long Run directory, and is generated if the 
#   user has specified either allele count or allele frequeny 
#   data. 
#   Note that it is not generated if the user has specified the 
#   sample covariance.
#Rerun with those adjustments??
#I did run with those parameters, just how do I load it?

# Or the patterns of decay of covariance with 
#   geographic distance can be compared between 
#   the data and the model.
k <- 60 #complete guess
#Need to do the above thing in order to do this thing
plot(ex.output$last.params$D[1:k,1:k][index.matrix], 
     sample.covariance[index.matrix],
     pch=19,col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main="Model adequacy:\n IBD patterns")
points(ex.output$last.params$D[1:k,1:k][index.matrix], 
       MC.parametric.covariance[index.matrix],col="red",pch=20)
legend(x="topright",pch=19,col=c(1,2),
       legend=c("observed","model estimate"))


#Try loading the MCMC
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_624100/simtest_run1_st2_LongRun/")
data(simtest_run1_st2_MCN.frequencies.list)
data(ex.output)

weather <- dget("~/ecs_unit5/weather.robj")
testfreqout <- dget("simtest_run1_st2_MCN.frequencies.list.Robj")
load("/Output/theStats.Robj"); str(the.stats)

testfreqout <- load("simtest_run1_st2_MCN.frequencies.list.Robj")
str(testfreqout)
#Doing Open In... Might work

#Try the tutorial
sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies), use="pairwise.complete.obs")

View(simtest_run1_st2_space_MCMC_output1)
#Maybe not, this did not load in

testother <- load("simtest_run1_st2_seed.Robj", verbose = TRUE)
#Look at documentation for make.spacemix.map.list
#How is it reading it back in?
sol2_nm.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "./run_275884/soltest_run2_nm_LongRun/soltest_run2_nm_space_MCMC_output1.Robj",
                                                    geographic.locations = population.coordinates,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)

#MCMC.output <- load_MCMC_output(MCMC.output.file)
load_MCMC_output <- function(MCMC.output.file){
  tmpenv <- environment()
  tmp <- load(MCMC.output.file,envir=tmpenv)
  mcmc.output <- lapply(tmp,get,envir=tmpenv)
  names(mcmc.output) <- tmp
  return(mcmc.output)
}

getwd()
MCMC.output <- load_MCMC_output("simtest_run1_st2_space_MCMC_output1.Robj")
#Okay, so for one, is it possible that that is all I need and I do not need the other outputs
#Secondly, to look at the MCMC freq, can I use this same loading mechanism?

MCMC.freq <- load_MCMC_output("simtest_run1_st2_MCN.frequencies.list.Robj")
#Heck yeah!!! Awesome!!!
#Screw their documentation for that to have taken so long to figure out!


#### Graphs for Report ####
#Highlight/change text color?
#Instead, try to make colors make more sense
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0522/Output/best_runs")
species <- "B"
model <- "source_and_target"
folder <- paste(species,"_",model,"_LongRun",sep="")
mcmcfile <- paste("./",folder,"/",species,"_",model,"_space_MCMC_output1.Robj",sep="")
locations2 <- read.csv(paste("../../Input/locations_",species,".csv",sep=""), header=FALSE)
locations2 <- matrix(unlist(locations2), ncol =2, nrow =dim(locations2)[1])
population.coordinates <- locations2

#Do colors by choice
sample.colors <- rainbow(n=dim(locations2)[1],start=0.1,end=0.8)[as.numeric(cut(population.coordinates[,1],dim(locations2)[1]))]

samples_B <- unlist(c('X413', 'X536', 'BAR', 'BRW', 'CGC', 'FNJ', 'GBE', 'MW', 'NAT', 'PLY', 'PT', 'PVT', 'SLI'))
samples_A <- unlist(c('BLP', 'MCX', 'SHN', 'SsLI1k', 'SsLI23k', 'SsLI4k', 'WV'))
samples_sim <- unlist(c('NLI', 'GA', 'CTM', 'WFH', 'ELPPB', 'MW'))

colors_sim <- unlist(c('#4ffa1b', '#fa991b', '#1b6dfa', '#1b6dfa', '#1b6dfa', '#1b6dfa'))
colors_A <- unlist(c('#9dc922', '#4ffa1b', '#4ffa1b', '#9dc922', '#9dc922', '#4ffa1b', '#1b6dfa'))

colors_B <- unlist(c('X413', 'X536', 'BAR', 'BRW', 'CGC', 'FNJ', 'GBE', 'MW', 'NAT', 'PLY', 'PT', 'PVT', 'SLI'))
colors_B <- unlist(c('#fa741b', '#DC27FF', '#1bb7fa', '#1bb7fa', '#1bb7fa', '#fa741b', '#DC27FF', '#1b6dfa', '#7f18f5', '#1bb7fa', '#1bb7fa', '#1bb7fa', '#4ffa1b'))
#GBE
#fa19eb
#NJ
#'#fa741b'
#Southern coast of MA
#1bb7fa
#Cape Cod Bay
#1b6dfa


sample.colors <- get(paste("colors_", species, sep=""))

sample.names <- get(paste("samples_", species, sep="")) #Use get to turn a pasted string variable name into a output from an actual variable

#Set up Plot
mapname <- paste(species, "_",model,"_","mixmap",sep="")
assign(paste(mapname,sep=""),make.spacemix.map.list(MCMC.output.file = mcmcfile, 
     geographic.locations = population.coordinates,
     name.vector = sample.names,
     color.vector = sample.colors,
     quantile=0.95,
     burnin=0))
#Print Plots - for Real saving
make.spacemix.map(spacemix.map.list = get(mapname),
                  #text=TRUE,
                  ellipses=TRUE,
                  source.option=TRUE,
                  #xlim=c(-74,-69),ylim=c(39,43)) #Sol B
                  xlim=c(-74,-68.5),ylim=c(39,44)) #Sol B zoom out
                  #xlim=c(-74.5,-71),ylim=c(39.5,41.5)) #Sol A basic
                  #xlim=c(-74.5,-71),ylim=c(39,42)) #Sol A zoomed out
                  #xlim=c(-82,-68),ylim=c(31,45)) #For sim
query.spacemix.map("GA",
#query.spacemix.map(c('NLI', 'CTM', 'WFH', 'ELPPB', 'MW'),
#query.spacemix.map(c('BAR', 'BRW', 'CGC', 'MW', 'PLY', 'PT', 'PVT', 'SLI'),
#query.spacemix.map(c('BLP', 'MCX', 'SHN', 'SsLI1k', 'SsLI23k', 'SsLI4k', 'WV'),
#query.spacemix.map(c('BLP', 'MCX', 'SHN', 'SsLI1k', 'SsLI23k', 'SsLI4k', 'WV'),
                   get(mapname),
                   ellipses = FALSE,
                   source.option = TRUE)
plot.admix.arrows(admix.source.coords=get(mapname)$MAPP.admix.source.coords,
                  geogen.coords=get(mapname)$MAPP.geogen.coords,
                  admix.proportions = get(mapname)[["MCMC.output"]][["admix.proportions"]],
                  colors=NULL,
                  length=0.2)
#Figured out how to do arrows myself by setting lwd
arrows(x0 = get(mapname)$MAPP.admix.source.coords[,1],
       y0 = get(mapname)$MAPP.admix.source.coords[,2],
       x1 = get(mapname)$MAPP.geogen.coords[,1],
       y1 = get(mapname)$MAPP.geogen.coords[,2],
       lwd = get(mapname)[["MCMC.output"]][["admix.proportions"]] * 100, #Increase the visibility
       col = get(mapname)$color.vector,
       length=0.2)





#### Output Likelihood / Acceptances ####
#Is there a likelihood/probability output?
#I saw some in the logs for the fast runs, but where is that put later and is that helpful?

###IBD Estimate
#For the st model for similis
getwd()
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_624100/simtest_run1_st2_LongRun")
#Set names to match tutorial for ease
ex.output <- load_MCMC_output("simtest_run1_st2_space_MCMC_output1.Robj")
MCN.freq <- load_MCMC_output("simtest_run1_st2_MCN.frequencies.list.Robj")
MCN.frequencies.list <- MCN.freq$MCN.frequencies.list #Go a level deeper

#and normalized sample allele frequencies.
sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                         use="pairwise.complete.obs")
k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/ex.output$last.params$inv.mean.sample.sizes / 
                                (sum(1/ex.output$last.params$inv.mean.sample.sizes)),
                              nrow=k,ncol=k,byrow=TRUE)

MC.parametric.covariance <- (MC.matrix) %*%     
  ex.output$last.params$admixed.covariance %*% 
  t(MC.matrix)
index.matrix <- upper.tri(sample.covariance,diag=TRUE)
plot(sample.covariance[index.matrix], 
     MC.parametric.covariance[index.matrix],
     col=adjustcolor("black",0.3),pch=20,
     xlab="Sample covariance",
     ylab="Parametric covariance",
     main="Model adequacy:\n matrix comparison:\n Similis S/T Model")
abline(0,1,col="red")
#Heck yeah

#Now the IBD thing
plot(ex.output$last.params$D[1:k,1:k][index.matrix], 
     sample.covariance[index.matrix],
     pch=19,col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main="Model adequacy:\n IBD patterns:\n Similis S/T Model")
points(ex.output$last.params$D[1:k,1:k][index.matrix], 
       MC.parametric.covariance[index.matrix],col="red",pch=20)
legend(x="topright",pch=19,col=c(1,2),
       legend=c("observed","model estimate"))

##Repeat with others
target <-"/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_445102/simtest_run1_target2_LongRun"
source <- "/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_204492/simtest_run1_sorce_LongRun"
nm <- "/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_262311/simtest_run1_nm_LongRun"

setwd(nm)
#Set names to match tutorial for ease
ex.output <- load_MCMC_output("simtest_run1_nm_space_MCMC_output1.Robj")
MCN.freq <- load_MCMC_output("simtest_run1_nm_MCN.frequencies.list.Robj")
MCN.frequencies.list <- MCN.freq$MCN.frequencies.list #Go a level deeper

sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                         use="pairwise.complete.obs")
k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/ex.output$last.params$inv.mean.sample.sizes / 
                                (sum(1/ex.output$last.params$inv.mean.sample.sizes)),
                              nrow=k,ncol=k,byrow=TRUE)

MC.parametric.covariance <- (MC.matrix) %*%     
  ex.output$last.params$admixed.covariance %*% 
  t(MC.matrix)
index.matrix <- upper.tri(sample.covariance,diag=TRUE)
plot(sample.covariance[index.matrix], 
     MC.parametric.covariance[index.matrix],
     col=adjustcolor("black",0.3),pch=20,
     xlab="Sample covariance",
     ylab="Parametric covariance",
     main="Model adequacy:\n matrix comparison:\n Similis No Movement Model")
abline(0,1,col="red")

plot(ex.output$last.params$D[1:k,1:k][index.matrix], 
     sample.covariance[index.matrix],
     pch=19,col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main="Model adequacy:\n IBD patterns:\n Similis No Movement Model")
points(ex.output$last.params$D[1:k,1:k][index.matrix], 
       MC.parametric.covariance[index.matrix],col="red",pch=20)
legend(x="topright",pch=19,col=c(1,2),
       legend=c("observed","model estimate"))

###Run the rest of the stuff from the tutorial
getwd()
setwd("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_624100/simtest_run1_st2_LongRun")
#Set names to match tutorial for ease
ex.output <- load_MCMC_output("simtest_run1_st2_space_MCMC_output1.Robj")
# Trace plot of posterior probability over the MCMC
plot(ex.output$Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot",type='l')
# Trace plots of alpha parameters of the spatial covariance function
matplot(t(ex.output$nugget),type='l',
        xlab="MCMC iterations",ylab="Parameter value",
        main="Trace plot of nugget parameters")
# Joint marginal plot of a0 and a1
#   colored by where in the MCMC these 
#   parameters took their values
plot(ex.output$a0,ex.output$a1,xlab="a0",ylab="a1",
     main="Joint marginal of a0 and a1",pch=20,
     col=adjustcolor(rainbow(1000,start=4/6,end=6/6),0.3))
legend(x="bottomright",pch=19,cex=0.8,
       col=rainbow(1000,start=4/6,end=6/6)[c(1,500,1000)],
       legend=c("Sampled MCMC iteration 1",
                "Sampled MCMC iteration 500",
                "Sampled MCMC iteration 1000"))
# Acceptance rate of a0 over the course of the 
#   MCMC analysis
plot(ex.output$accept_rates$a0_accept_rate,
     xlab="MCMC iterations",ylab="Acceptance rate",
     main="Acceptance rate of a0",type='l',
     ylim=c(0.35,0.6))
abline(h=0.44,col="gray",lty=2)
# Acceptance rates of nugget parameters over the 
#   course of the MCMC analysis
matplot(t(ex.output$accept_rates$nugget_accept_rate),
        xlab="MCMC iterations",ylab="Acceptance rate",
        main="Acceptance rates of nuggets",type='l',
        ylim=c(0.3,0.7))
abline(h=0.44,col="gray",lty=2)

sim_st2.output <- load_MCMC_output("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_624100/simtest_run1_st2_LongRun/simtest_run1_st2_space_MCMC_output1.Robj")
sim_target2.output <- load_MCMC_output("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_445102/simtest_run1_target2_LongRun/simtest_run1_target2_space_MCMC_output1.Robj")
sim_nm.output <- load_MCMC_output("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_262311/simtest_run1_nm_LongRun/simtest_run1_nm_space_MCMC_output1.Robj")
sim_source.output <- load_MCMC_output("/Users/hannah/gitHub/Spisula/April_2022/SpaceMixAttempt/New_Runs_0506/Output/run_204492/simtest_run1_sorce_LongRun/simtest_run1_sorce_space_MCMC_output1.Robj")

plot(sim_st2.output$Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot\n Similis S/T",type='l')
plot(sim_target2.output$Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot\n Similis Target",type='l')
plot(sim_source.output$Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot\n Similis Source",type='l')
plot(sim_nm.output$Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot\n Similis No Movement",type='l')

### Now, what about Likelihoods?
#I should sample more regularly so my plots look more like the examples rather than 
#Average of the Prob from output
mean(sim_nm.output$Prob)
7042.481



