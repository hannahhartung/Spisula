#Neighbor Joining Tree

#Note: You can improve the species coloring with:https://yulab-smu.top/treedata-book/chapter6.html
#Grouping OTUs

### Libraries
library('vcfR')
library(adegenet)
library(ape)
library("treeio")
library("ggtree")
library(ggplot2)
library(dplyr)

#### Create and Plot Population Tree ####
vcf <- read.vcfR("/Users/hannah/gitHub/Spisula/April_2022/divmigrate/SNP_solrefall_0528_LDprune.vcf")
#Use above for normal LD and below for LD 0.2
vcf <- read.vcfR("/Users/hannah/gitHub/Spisula/April_2022/LD_prune/LD02/SNP_solrefall_LD02prune.vcf")
x <- vcfR2genlight(vcf)
#Get a list of popnames in order
df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")
df2 <- df %>%
  filter(label %in% x$ind.names)
df3 <- data.frame(label=x$ind.names)
df4 <- full_join(df3, df2, by='label')
#Check if they are in the same order
df5 <- data.frame(label3=df3$label,label2=df2$label,pop2=df2$POP)
#df5[,1] == df5[,2] #all TRUE
x <- vcfR2genind(vcf, return.alleles = TRUE)
#my_genclone <- poppr::as.genclone(x)
#pop(my_genclone) <- df5[,3] #Works but can I do it in the genind
pop(x) <- df5[,3]
genpopx <- genind2genpop(x)
genpopx.dist <- dist.genpop(genpopx, method = 1, diag = FALSE, upper = FALSE)
#Nei's distance (not Euclidean) = method 1
tr <- nj(genpopx.dist)
ggtree(tr) + 
  #geom_tippoint(size=0.5,shape=15) +
  #Plot sample names to confirm the colors are correct
  geom_tiplab(size=2.5, color="black")
#UPGMA for poptree
library(phangorn)
tr4 <- upgma(genpopx.dist)
ggtree(tr4, layout = 'rectangular', ladderize=TRUE) + 
  #geom_tippoint(size=0.5,shape=15) +
  #Plot sample names to confirm the colors are correct
  geom_tiplab(size=1, color="black") + 
  theme(aspect.ratio = 0.5) + geom_treescale()

##########
#Date: August 1st
#Try using phangorn to make a bootstrapped or ML tree from the TreeAnnotator output
#temp <- read.nexus.data("Aug_2022/SNAPP/c_snapper_annotree.nex")
#Did not work

#Date: Aug 5th
#Try ape
tree_c <- ape::read.nexus("Aug_2022/SNAPP/c_snapper_annotree.nex")
#Loaded in!
plot(tree_c)
#Just shows branch lengths like densitree did. What could I do in ape with the output from SNAPP.
#Xml output
#Convert from xml to ape with:
# install.packages("devtools")
#devtools::install_github("USCBiostats/rphyloxml")
#https://uscbiostats.github.io/rphyloxml/
#Too glitchy

#########

#one of 'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' or 'ape'

### Create Ind Tree with ape ###
x.dist <- poppr::bitwise.dist(x)
tr2 <- nj(x.dist)
df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")
df2 <- df %>%
  filter(label %in% x$tip.label)
#d$var1 <- df2[,2]
tr3 <- full_join(tr2, df2, by='label')
cols <- c("#3399ff","#cc66ff","#6600ff","#ffe680","#ccff33","#009900")
#Plot
options(ignore.negative.edge=TRUE)
ggtree(tr3) + 
  geom_tippoint(aes(colour=POP),size=0.5,shape=15) +
  #Plot sample names to confirm the colors are correct
  #geom_tiplab(size=0.5, color="black") +
  scale_color_manual(values=cols)


#### Plot Individual Tree from vcf2PopTree Newick Output ####
#NJ plot
setwd("/Users/hannah/gitHub/Spisula/April_2022/neighboor")
x <- read.tree("nwk_tree1.txt")
#Since it is unrooted, the first value is NaN so set to 0?
  #If you don't do this, it is a cladogram with no lengths
x$edge.length[1] <- 0
#Import population information
df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")
df2 <- df %>%
  filter(label %in% x$tip.label)
#d$var1 <- df2[,2]
tree <- full_join(x, df2, by='label')
#Set colors that make more sense
cols <- c("#3399ff","#cc66ff","#6600ff","#ffe680","#ccff33","#009900")
#Plot
options(ignore.negative.edge=TRUE)
ggtree(tree) + 
  geom_tippoint(aes(colour=POP),size=0.5,shape=15) +
#Plot sample names to confirm the colors are correct
  #geom_tiplab(size=0.5, color="black") +
  scale_color_manual(values=cols)



##############################################################


#### Tutorials/Notes ####
### Create tree at population level
#From: https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html
library('vcfR')
library(adegenet)
library(ape)
vcf <- read.vcfR("/Users/hannah/gitHub/Spisula/April_2022/divmigrate/SNP_solrefall_0528_LDprune.vcf")
vcf
x <- vcfR2genlight(vcf)
gt <- extract.gt(vcf, element = "GT")
gt[c(2,6,18), 1:3]
t(as.matrix(x))[c(1,5,17), 1:3]

#Get a list of popnames in order
df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")
df2 <- df %>%
  filter(label %in% x$ind.names)
df3 <- data.frame(label=x$ind.names)
df4 <- full_join(df3, df2, by='label')
#Check if they are in the same order
  #df5 <- data.frame(label3=df3$label,label2=df2$label,pop2=df2$POP)
  #df5[,1] == df5[,2] #all TRUE
pop(x) <- as.factor(df5[,3])
popNames(x)
ploidy(x) <- 2
x.dist <- dist(x)
x.dist <- poppr::bitwise.dist(x)

#Go back to poppr
#https://rdrr.io/cran/adegenet/man/dist.genpop.html
#input is genpop
vcf <- read.vcfR("/Users/hannah/gitHub/Spisula/April_2022/divmigrate/SNP_solrefall_0528_LDprune.vcf")
x <- vcfR2genind(vcf, return.alleles = TRUE)
my_genclone <- poppr::as.genclone(x)
pop(my_genclone) <- df5[,3] #Works but can I do it in the genind
pop(x) <- df5[,3]
genpopx <- genind2genpop(x)
genpopx.dist <- dist.genpop(genpopx, method = 1, diag = FALSE, upper = FALSE)
#Nei's distance (not Euclidean) = method 1
tr <- nj(genpopx.dist)
ggtree(tr) + 
  #geom_tippoint(size=0.5,shape=15) +
  geom_tiplab(size=6) +
  scale_color_manual(values=cols) + 
  geom_hilight(node=4, fill="steelblue", alpha=.6) +
  geom_hilight(node=9, fill="darkgreen", alpha=.6) +
  geom_hilight(node=1, fill="grey", alpha=.6) +
  geom_hilight(node=2, fill="lightpink", alpha=.6)
d <- data.frame(node=c(5, 9), type=c("A", "B"))
ggtree(tr) + geom_hilight(data=d, aes(node=node, fill=type),
                          type = "gradient", 
                          gradient.direction = 'rt')
ggtree(tr) + geom_hilight(data=d, mapping=aes(node=node, fill=type),
               type = "gradient", gradient.direction = 'tr',
               alpha = .8)


#Calculate for populations: https://rdrr.io/cran/dartR/man/gl.dist.pop.html
library(dartR)
#gl.install.vanilla.dartR()
glx <- gl.compatability.check(glx)
x.dist <- gl.dist.pop(x)
tr <- nj(x.dist)
#### GO TO BOTTOM - RUN THE SCRIPT LINE BY LINE
#Set my own version of it
library(report)
m <- gl2.dist.pop(x)
#This is TOO MUCH DEBUGGING
#https://dyerlab.github.io/applied_population_genetics/genetic-distances.html
library(devtools)
#install_github("dyerlab/gstudio")
#install_github("dyerlab/popgraph")
library(gstudio)
data(arapat)
D1 <- genetic_distance(x,stratum = "POP",mode = "Nei")
#NO MANUAL - ALSO NO GO - COULD EMAIL THEM MAYBE


df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")
df2 <- df %>%
  filter(label %in% x$ind.names)
df3 <- data.frame(label=x$ind.names)
df4 <- full_join(df3, df2, by='label')
#Check if they are in the same order
#df5 <- data.frame(label3=df3$label,label2=df2$label,pop2=df2$POP)
#df5[,1] == df5[,2] #all TRUE
pop(x) <- as.factor(df5[,3])
popx <- genind2genpop(x, as.factor(df5[,3]))


df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")
df2 <- df %>%
  filter(label %in% tr$tip.label)
tree <- full_join(tr, df2, by='label')
#Set colors that make more sense
cols <- c("#3399ff","#cc66ff","#6600ff","#ffe680","#ccff33","#009900")
#Plot
options(ignore.negative.edge=TRUE)
ggtree(tree) + 
  geom_tippoint(aes(colour=POP),size=0.5,shape=15) +
  #Plot sample names to confirm the colors are correct
  #geom_tiplab(size=0.5, color="black") +
  scale_color_manual(values=cols)

#
tree <- full_join(x, df2, by='label')



plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- nj(dist.dna(woodmouse))
plot(trw)

### Individual Level
#Do it from vcf2poptree output and ggtree
#From: https://yulab-smu.top/treedata-book/chapter4.html
library("treeio")
library("ggtree")
library(ggplot2)
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
ggtree(tree, color="firebrick", size=2, linetype="dotted")
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
#Woot!
#Now input mine as possible
setwd("/Users/hannah/gitHub/Spisula/April_2022/neighboor")
tree <- read.tree("nwk_tree1.txt")
ggtree(tree, color="firebrick", size=0.5, linetype="dotted")
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
#Aw hell yeah! We got 'em now 
ggtree(tree, color="firebrick", size=0.5, linetype="dotted",branch.length="none", ladderize=FALSE)
#The vcf2poptree does not have edge.legth information


#Practice merging trait/population data with it

#instead color the tips
#install.packages("ggnewscale")
library(ggnewscale)

set.seed(2020)
x <- rtree(30)
d <- data.frame(label=x$tip.label, var1=abs(rnorm(30)), var2=abs(rnorm(30)))
tree <- full_join(x, d, by='label')
trs <- list(TREE1 = tree, TREE2 = tree)
class(trs) <- 'treedataList'
ggtree(trs) + facet_wrap(~.id) + 
  geom_tippoint(aes(subset=.id == 'TREE1', colour=var1)) + 
  scale_colour_gradient(low='blue', high='red') +  
  ggnewscale::new_scale_colour()  + 
  geom_tippoint(aes(colour=var2), data=td_filter(.id == "TREE2")) + 
  scale_colour_viridis_c()

#Now me
setwd("/Users/hannah/gitHub/Spisula/April_2022/neighboor")
x <- read.tree("nwk_tree1.txt")
df <- data.frame(read.table("samplist_tree0621.txt"))
rownames(df) <- df$V1
colnames(df) <- c("label","POP")

#Try using their method for trait values
d <- data.frame(label=x$tip.label, var1=abs(rnorm(484)), var2=abs(rnorm(484)))
  #Worked, now set var1 to pop
  #The issue could be having too many values

#Filter a dataframe by values only present in the other one
library(dplyr)
df2 <- df %>%
  filter(label %in% d$label)

d$var1 <- df2[,2]
tree <- full_join(x, d, by='label')

trs <- list(TREE1 = tree, TREE2 = tree)
class(trs) <- 'treedataList'
ggtree(trs) + facet_wrap(~.id) + 
  geom_tippoint(aes(subset=.id == 'TREE1', colour=var1)) + 
  #scale_colour_gradient(low='blue', high='red') +  
  ggnewscale::new_scale_colour()  + 
  geom_tippoint(aes(colour=var2), data=td_filter(.id == "TREE2")) + 
  scale_colour_viridis_c()

#Plot just the one
ggtree(tree) + geom_tippoint(aes(colour=var1),size=0.2)



#Runs for a VERY long time, maybe try turning the pop names into numbers
df1 <- df
df1[df1[,2] == "Boff",2] <- 1.0
df1[df1[,2] == "Bin",2] <- 2.0
df1[df1[,2] == "MA",2] <- 4.0
df1[df1[,2] == "GA",2] <- 5.0
df1[df1[,2] == "NLI",2] <- 6.0
df1[df1[,2] == "A",2] <- 3.0 #Gotta do A after GA and MA
tree2 <- full_join(x, df1, by='label')
ggtree(tree2)
#Changing it to numbers still takes forever
tree <- read.tree("nwk_tree1.txt")
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
ggplot(tree2, aes(x, y)) + geom_tree() + theme_tree()

trs <- list(TREE1 = tree2, TREE2 = tree2)
ggtree(trs) + facet_wrap(~.id) + 
  geom_tippoint(aes(subset=.id == 'TREE1', colour=POP)) + 
  scale_colour_gradient(low='blue', high='red') +  
  ggnewscale::new_scale_colour()  + 
  geom_tippoint(aes(colour=POP), data=td_filter(.id == "TREE2"))  
  #scale_colour_viridis_c()


#Don't do this cause it colors the whole branch
#install.packages("TDbook")
#install.packages("phytools")

library(TDbook)
svl <- as.matrix(df_svl)[,1]
fit <- phytools::fastAnc(tree_anole, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree_anole, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree_anole, d, by = 'node')
p1 <- ggtree(tree, aes(color=trait), layout = 'circular', 
             ladderize = FALSE, continuous = 'colour', size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 
p1
#Color colors the whole branch from ancestral rather than just doing the last branch - that might be unchangable
#However, I could instead do the tips
p <- ggtree(tree) + 
  geom_nodepoint(color="#b5e521", alpha=1/4, size=10) 
p + geom_tippoint(color="#FDAC4F", shape=8, size=3)

#But I still need to know how to add other information to the tree info for aes
#svl is a named number (created by converting with as.matrix)
#Otherwise first it can just be imported as a df
df_svl
df <- read.table("samplist_tree0621.txt")
#Set names
rownames(df) <- df$V1
svl1 <- as.matrix(df)[,2]
#Put together
td <- data.frame(node = nodeid(tree, names(svl1)),
                 trait = svl1)
nd <- data.frame(node = names(fit$ace), trait = fit$ace) #Uhh, doesn't do anything? But still works?
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree2 <- full_join(tree, d, by = 'node')
p1 <- ggtree(tree2, aes(color=trait), 
            continuous = 'colour', size=0.5) #+
  #scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  #geom_tiplab(hjust = -.1) + 
  #xlim(0, 1.2) + 
  #theme(legend.position = c(.05, .85)) 
p1


#Figuring out the installation and updating
#Warning in install.packages :
#'lib = "ggtree"' is not writable
#Would you like to use a personal library instead? (yes/No/cancel) yes
#Warning in install.packages :
#  package ‘treeio’ is not available (for R version 3.6.2)
#Oh! When it says things aren't avaliable, maybe I don't have thier repository selected:
setRepositories()
#Ah, now it's saying ggtree was built for R 3.6.3
#Still doesn't work after updating R
#On it's rdrr.io page, it has special instructions for installing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
#Requires compilation/reinstallation on a bunch of other stuff too cause I updated
  #Try installing cmake to make compliation installs more likely to work?
library("treeio")
#Cite: LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package for phylogenetic tree input and output with richly annotated and associated data. Molecular Biology and Evolution 2019, accepted. doi: 10.1093/molbev/msz240
library("ggtree")
#Error: package or namespace load failed for ‘ggtree’:
    #object ‘get_aes_var’ is not exported by 'namespace:rvcheck'
#The reason you're running into this error is because the latest version of rvcheck (0.2.0) has removed the get_aes_var function. The current ggtree version (3.0.4) is aware of this change, and doesn't look for it. However, because you're using an old version of R, you're using also using an old version of both Bioconductor and subsequently the ggtree package. This older version still thinks get_aes_var should be present, and then fails when it isn't.

#Okay, so maybe do the big update but before I try that, just see if it could magically work...?
# Error in base::nchar(wide_chars$test, type = "width") : 
    #lazy-load database '/Users/hannah/Library/R/3.6/library/cli/R/sysdata.rdb' is corrupt
#So nah, okay

setwd("/Users/hannah/gitHub/Spisula/April_2022/")
nwk <- system.file("neighboor","nwk_tree1.txt", package="treeio")
#Oh system file is to read their example from the package itself
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
ggtree(tree, color="firebrick", size=2, linetype="dotted")
library(ggplot2)
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

#Do it from plink data
setwd("/Users/hannah/gitHub/Spisula/April_2022/neighboor")
M <- read.table("plink.dist")
ids <- read.table("plink.real.ids.txt")
dimnames(M) <- list(ids[,1],ids[,1])
M <- as.matrix.data.frame(M)
tr <- nj(M)
tre <- ladderize(tr)
plot(tre, cex=.6)
myBoots <- boot.phylo(tre, dna, function(e) root(nj(dist.dna(e, model = "TN93")),1))

#https://search.r-project.org/CRAN/refmans/ape/html/nj.html
#https://www.molecularecologist.com/2016/02/26/quick-and-dirty-tree-building-in-r/

library(ape)
### From Saitou and Nei (1987, Table 1):
x <- c(7, 8, 11, 13, 16, 13, 17, 5, 8, 10, 13,
       10, 14, 5, 7, 10, 7, 11, 8, 11, 8, 12,
       5, 6, 10, 9, 13, 8)
M <- matrix(0, 8, 8)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(1:8, 1:8)
tr <- nj(M)
plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- nj(dist.dna(woodmouse))
plot(trw)

## Line by line pop
#gl.dist.pop https://rdrr.io/cran/dartR/src/R/gl.dist.pop.r

gl2.dist.pop <- function(x,
                        method = "euclidean",
                        plot.out = TRUE,
                        scale = FALSE,
                        output="dist",
                        plot_theme = theme_dartR(),
                        plot_colors = two_colors,
                        save2tmp = FALSE,
                        verbose = NULL) {
  
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "reshape2"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it."
    ))
  }
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  
  # CHECK DATATYPE
  datatype <-
    utils.check.datatype(x, accept = c("SNP","SilicoDArT"), verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # DO THE JOB
  
  available_methods <-
    c(
      "euclidean",
      "nei",
      "reynolds",
      "chord",
      "fixed-diff"
    )
  
  if (!(method %in% available_methods)) {
    cat(error("Fatal Error: Specified distance method is not among those 
                available.\n"))
    stop("Specify one of ",paste(available_methods, 
                                 collapse = ", ")," or fixed-diff.\n")
  }
  # hard.min.p <- 0.25
  
  nI <- nInd(x)
  nL <- nLoc(x)
  nP <- nPop(x)
  dd <- array(NA, c(nPop(x), nPop(x)))
  
  # Calculate distances 
  # if (method %in% distmethod) {
  if (verbose >= 2) {
    cat(report(paste(
      "  Calculating distances: ", method, "\n"
    )))
    cat(report(
      "  Refer to the dartR Distance Analysis tutorial for algorithms\n"
    ))
  }
  
  # Calculate allele frequencies for each population and locus
  f <- gl.percent.freq(x, verbose = 0)
  # Select only pop, locus, frequency columns
  f <- f[, c(1, 2, 6)]
  # Convert to a pop x locus matrix
  f <- reshape2::dcast(f, popn ~ locus, value.var = "frequency")
  # Reassign names to the populations, and convert from percentages to proportions
  row.names(f) = f[, 1]
  f <- f[,-c(1)]
  p <- f / 100
  
  # For both DArTseq and SilicoDArT
  if (method == "euclidean") {
    for (i in (1:(nP - 1))) {
      for (j in ((i + 1):nP)) {
        p_ind1 <- p[i,]
        p_ind2 <- p[j,]
        sq <- (p_ind1-p_ind2)**2
        sq <- sq[!is.na(sq)]
        L <- length(sq)
        if(scale==TRUE){
          if(datatype=="SNP"){
            dd[j,i] <- 0.5*sqrt(sum(sq)/L)
          } else {
            dd[j,i] <- sqrt(sum(sq)/L)
          }
        } else {
          dd[j,i] <- sqrt(sum(sq))
        }
      }
    }
  }
  # # Test code
  # x <- dartR::gl2gi(testset.gl)
  # x <- adegenet::genind2genpop(x)
  # D_check <- adegenet::dist.genpop(x,4) # Rogers D
  # hist(D_check,breaks=50)
  # D <- dartR::gl.dist.pop(testset.gl, method='euclidean',output="matrix",scale=TRUE)
  # D[upper.tri(D)] <- t(D)[upper.tri(D)]
  # hist(D/2,breaks=50)
  # #VALIDATED [with minor differences, missing handling?]
  
  # For DArTseq only
  if (method == "reynolds") {
    if(datatype=="SilicoDArT"){
      stop(error("Fatal Error: Reynolds Distance is not available 
                       for presence-absence data\n"))
    }
    for (i in (1:(nP - 1))) {
      for (j in ((i + 1):nP)) {
        # Pull the loci for individuals i and j
        pind1 <- p[i,]
        pind2 <- p[j,]
        # Delete the pairwise missing
        tmp <- pind1+pind2
        pind1 <- pind1[!is.na(tmp)]
        pind2 <- pind2[!is.na(tmp)]
        # Squares
        psq <- (pind1-pind2)**2
        # Repeat for q
        qind1 <- 1-pind1
        qind2 <- 1-pind2
        qsq <- (qind1-qind2)**2
        # Cross products
        p12 <- pind1*pind2
        q12 <- qind1*qind2
        # Non-missing loci
        #L <- length(psq)
        
        #dd[j,i] <- sqrt(sum(psq+qsq)/(2*sum(1-p12-q12)))
        dd[j,i] <- -log(1-sqrt(sum(psq+qsq)/(2*sum(1-p12-q12))))
        #dd[j,1] <- sqrt(sum(psq)/(sum(1-p12-q12)))
      }
    }
  }
  # # Test code
  # x <- dartR::gl2gi(testset.gl)
  # x <- adegenet::genind2genpop(x)
  # D_check <- adegenet::dist.genpop(x,3) # Reynolds in common use
  # D_check <- -log(1-D_check) # Proportional to divergence time
  # hist(D_check,breaks=50)
  # D <- dartR::gl.dist.pop(testset.gl, method='reynolds',output='matrix',scale=TRUE)
  # D[upper.tri(D)] <- t(D)[upper.tri(D)]
  # hist(D,breaks=50)
  # #VALIDATED [with minor difference, missing handling?]
  
  if (method == "nei") {
    if(datatype=="SilicoDArT"){
      stop(error("Fatal Error: Nei Standard Distance is not available
                       for presence-absence data\n"))
    }
    for (i in (1:(nP - 1))) {
      for (j in ((i + 1):nP)) {
        # Pull the loci for individuals i and j
        prow1 <- p[i,]
        prow2 <- p[j,]
        # Delete the pairwise missing
        tmp <- prow1+prow2
        prow1 <- prow1[!is.na(tmp)]
        prow2 <- prow2[!is.na(tmp)]
        # Squares
        p1sq <- prow1*prow1
        p2sq <- prow2*prow2
        # Repeat for q=1-p
        qrow1 <- 1-prow1
        qrow2 <- 1-prow2
        q1sq <- qrow1*qrow1
        q2sq <- qrow2*qrow2
        # Cross products
        p12 <- prow1*prow2
        q12 <- qrow1*qrow2
        # Number of non-missing loci
        L <- length(p12)
        
        dd[j,i] <- -log(sum(p12+q12)/(sqrt(sum(p1sq+q1sq))*sqrt(sum(p2sq+q2sq))))
      }
    }
  }
  # # Test code
  # x <- dartR::gl2gi(testset.gl)
  # x <- adegenet::genind2genpop(x)
  # D_check <- adegenet::dist.genpop(x,1) 
  # hist(D_check,breaks=50)
  # D <- dartR::gl.dist.pop(testset.gl, method='nei',output='matrix',scale=TRUE)
  # hist(D,breaks=50)
  # #VALIDATED [with minor difference, missing handling?]
  
  if (method == "chord") {
    if(datatype=="SilicoDArT"){
      stop(error("Fatal Error: Czfordi-Edwards Chord Distance is not available
                       for presence-absence data\n"))
    }
    for (i in (1:(nP - 1))) {
      for (j in ((i + 1):nP)) {
        # Pull the loci for individuals i and j
        prow1 <- p[i,]
        prow2 <- p[j,]
        # Delete the pairwise missing
        tmp <- prow1+prow2
        prow1 <- prow1[!is.na(tmp)]
        prow2 <- prow2[!is.na(tmp)]
        # create proportions for allele 2
        qrow1 <- 1-prow1
        qrow2 <- 1-prow2
        # Cross products
        p12 <- prow1*prow2
        q12 <- qrow1*qrow2
        # Non-missing Loci
        L <- length(p12)
        
        dd[j,i] <- (2/pi)*sqrt(2*(1 - (sum(sqrt(p12))/L + sum(sqrt(q12)/L))))
      }
    }
  }
  # # Test code
  # x <- dartR::gl2gi(testset.gl)
  # x <- adegenet::genind2genpop(x)
  # D_check <- adegenet::dist.genpop(x,2) # Angular or Edwards?
  # #D_check <- -log(1-D_check) # Proportional to divergence time
  # hist(D_check,breaks=50)
  # D <- dartR::gl.dist.pop(testset.gl, method='chord',output='matrix',scale=TRUE)
  # D[upper.tri(D)] <- t(D)[upper.tri(D)]
  # hist(D,breaks=50)
  # #VALIDATED [with minor difference, missing handling?]
  
  if (method == "fixed-diff") {
    dd <- gl.fixed.diff(x, verbose = 0)[[3]]/100
    if (verbose >= 2) {
      cat(report("  Calculating proportion of fixed differences\n"))
      cat(
        warn(
          "Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n"
        )
      )
    }
  }
  
  # # Revert to original order ord <- rank(popNames(x)) mat <- as.matrix(dd)[ord, ord] dd <- as.dist(mat)
  
  if(method != "fixed-diff") {
    dimnames(dd) <- list(popNames(x), popNames(x))
  }
  
  # PLOT Plot Box-Whisker plot
  
  if (plot.out) {
    if (datatype == "SNP") {
      title_plot <- paste0("SNP data\nUsing ", method, " distance")
    } else {
      title_plot <-
        paste0("Tag P/A data (SilicoDArT)\nUsing ",
               method,
               " distance")
    }
    values <- NULL
    val <- as.vector(dd)
    val <- val[!is.na(val)]
    df_plot <- data.frame(values = val)
    
    # Boxplot
    p1 <- ggplot(df_plot, aes(y = values)) +
      geom_boxplot(color = plot_colors[1], fill = plot_colors[2]) +
      coord_flip()  +
      plot_theme  +
      xlim(range = c(-1,1)) + 
      ylim(min(df_plot$values, na.rm = TRUE),max(df_plot$values, na.rm = TRUE)) + 
      ylab(" ") + 
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) + 
      ggtitle(title_plot)
    
    # Histogram
    p2 <- ggplot(df_plot, aes(x = values)) +
      geom_histogram(bins = 20,color = plot_colors[1], fill = plot_colors[2]) +
      xlim(min(df_plot$values, na.rm = TRUE), max(df_plot$values, na.rm = TRUE)) +
      xlab("Distance Metric") +
      ylab("Count") +
      plot_theme
  }
  
  # SUMMARY Print out some statistics
  if (verbose >= 3) {
    cat("  Reporting inter-population distances\n")
    cat("  Distance measure:", method, "\n")
    cat("    No. of populations =", nPop(x), "\n")
    cat("    Average no. of individuals per population =",
        round(nInd(x) / nPop(x),1),
        "\n")
    cat("    No. of loci =", nLoc(x), "\n")
    cat("    Minimum Distance: ", round(min(dd,na.rm=TRUE), 2), "\n")
    cat("    Maximum Distance: ", round(max(dd,na.rm=TRUE), 2), "\n")
    cat("    Average Distance: ", round(mean(dd,na.rm=TRUE), 3), "\n")
  }
  
  # SAVE INTERMEDIATES TO TEMPDIR
  
  # creating temp file names
  if (save2tmp) {
    if (plot.out) {
      temp_plot <- tempfile(pattern = "Plot_")
      match_call <-
        paste0(names(match.call()),
               "_",
               as.character(match.call()),
               collapse = "_")
      # saving to tempdir
      saveRDS(list(match_call, p3), file = temp_plot)
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call, dd), file = temp_table)
    if (verbose >= 2) {
      cat(report("  Saving tabulation to session tempfile\n"))
      cat(
        report(
          "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
        )
      )
    }
  }
  
  # PRINTING OUTPUTS
  if (plot.out) {
    # using package patchwork
    p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
    suppressWarnings(print(p3))
  }
  
  if(output=="dist"){
    dd <- as.dist(dd)
    if(verbose >= 2){cat(report("  Returning a stats::dist object\n"))}
  } else {
    dd <- as.matrix(dd)
    if(verbose >= 2){cat(report("  Returning a square matrix object\n"))}
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(dd)
}
