source("/Users/hannah/gitHub/Spisula/src/plotting_funcs.R")
#New plots
#Date: June 13th
setwd("/Users/hannah/gitHub/Spisula/April_2022/treemix_out")
#Import
set <- "solmain"
run <- "_bm2" #or "_bm0globe" "_server"

stempath <- paste("/Users/hannah/gitHub/Spisula/April_2022/treemix_out/",set,run,sep="")

#Plot
plot_tree(stempath)


### Old plots ###
plot_tree("treemix/out_stem")
plot_tree("treemix/out_stem_allSNPs")

plot_resid("treemix/out_stem_allSNPs", "treemix/poporder.txt")

plot_resid("treemix/outstem", "treemix/testpop.txt")

#Test with no K and no rooting
plot_tree("treemix/out_stem_nr0nk") + title("nr0nk")
plot_resid("treemix/out_stem_nr0nk", "treemix/poporder.txt") + title("nr0nk")
plot_tree("treemix/out_stem_nr0nk_all") + title("nr0nk all")
plot_resid("treemix/out_stem_nr0nk_all", "treemix/poporder.txt") + title("nr0nk all")

#No root and k = 10
plot_tree("treemix/out_stem_nr0k10_all") + title("nr0k10 all")
plot_resid("treemix/out_stem_nr0k10_all", "treemix/poporder.txt") + title("nr0k10 all")
plot_tree("treemix/out_stem_nr0k10") + title("nr0k10")
plot_resid("treemix/out_stem_nr0k10", "treemix/poporder.txt") + title("nr0k10")
#Rooted on NLI

#Try what I actually intend
#Root on GA
plot_tree("treemix/out_stem_rGA_m0") + title("rGA_m0")
plot_resid("treemix/out_stem_rGA_m0", "treemix/poporder.txt") + title("rGA_m0")
plot_tree("treemix/out_stem_rGA_m1") + title("rGA_m1")
plot_resid("treemix/out_stem_rGA_m1", "treemix/poporder.txt") + title("rGA_m1")
plot_tree("treemix/out_stem_rGA_m2") + title("rGA_m2")
plot_resid("treemix/out_stem_rGA_m2", "treemix/poporder.txt") + title("rGA_m2")
plot_tree("treemix/out_stem_rGA_m3") + title("rGA_m3")
plot_resid("treemix/out_stem_rGA_m3", "treemix/poporder.txt") + title("rGA_m3")
plot_tree("treemix/out_stem_rGA_m4") + title("rGA_m4")
plot_resid("treemix/out_stem_rGA_m4", "treemix/poporder.txt") + title("rGA_m4")
plot_tree("treemix/out_stem_rGA_m5") + title("rGA_m5")
plot_resid("treemix/out_stem_rGA_m5", "treemix/poporder.txt") + title("rGA_m5")

#With one bootstrap
plot_tree("treemix/out_stem_rGA_m0b") + title("rGA_m0b")
plot_resid("treemix/out_stem_rGA_m0b", "treemix/poporder.txt") + title("rGA_m0b")
plot_tree("treemix/out_stem_rGA_m1b") + title("rGA_m1b")
plot_resid("treemix/out_stem_rGA_m1b", "treemix/poporder.txt") + title("rGA_m1b")
plot_tree("treemix/out_stem_rGA_m2b") + title("rGA_m2b")
plot_resid("treemix/out_stem_rGA_m2b", "treemix/poporder.txt") + title("rGA_m2b")
plot_tree("treemix/out_stem_rGA_m3b") + title("rGA_m3b")
plot_resid("treemix/out_stem_rGA_m3b", "treemix/poporder.txt") + title("rGA_m3b")
plot_tree("treemix/out_stem_rGA_m4b") + title("rGA_m4b")
plot_resid("treemix/out_stem_rGA_m4b", "treemix/poporder.txt") + title("rGA_m4b")


#With solidissima
plot_tree("treemix/out_stem_simsol") + title("simsol_m0")
plot_resid("treemix/out_stem_simsol", "treemix/popordersimsol.txt") + title("simsol_m0")
plot_tree("treemix/out_stem_simsolm1") + title("simsol_m1")
plot_resid("treemix/out_stem_simsolm1", "treemix/popordersimsol.txt") + title("simsol_m1")
plot_tree("treemix/out_stem_simsolm2") + title("simsol_m2")
plot_resid("treemix/out_stem_simsolm2", "treemix/popordersimsol.txt") + title("simsol_m2")
plot_tree("treemix/out_stem_simsolm3") + title("simsol_m3")
plot_resid("treemix/out_stem_simsolm3", "treemix/popordersimsol.txt") + title("simsol_m3")
library(RColorBrewer)
library(R.utils)

#SNE
plot_tree("treemix/out_solSNEm0") + title("solSNEm0")
plot_resid("treemix/out_solSNEm0", "treemix/popordersolSNE.txt") + title("simsol_m0")

plot_tree(as.factor("treemix/out_solSNEm1r")) + title("solSNEm1r")

plot_tree(as.factor("treemix/out_solSNEm1rb")) + title("solSNEm1rb")
plot_tree(as.factor("treemix/out_solSNE_7km0rb")) + title("solSNE_7km0rb")
plot_tree(as.factor("treemix/out_solSNE_7km1rb")) + title("solSNE_7km1rb")
plot_tree(as.factor("treemix/out_solSNE_7k3m1rb")) + title("solSNE_7k3m1rb")
plot_tree(as.factor("treemix/out_solSNE_7k3m1b")) + title("solSNE_7k3m1b")


#Test
plot_tree("treemix/intest_m2") + title("test2") #gave same errors

