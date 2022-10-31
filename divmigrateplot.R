###Import genepop
###Divmigrate

### Plot LD Prune Sets ###
#Date: June 13th
library(diveRsity)
setwd("/Users/hannah/gitHub/Spisula/April_2022/divmigrate")
setwd("/Users/hannah/gitHub/Spisula/April_2022/LD_prune/LD02")
#Import
set <- "solonly"
stats <- "gst" #d, gst, nm
bootstraps <- 1000
run <- paste(set,"_",stats,"b",bootstraps,sep = "")
divpath <- paste(set,"_LD02.txt",sep = "")
divout <- paste(run)

divMigrate(infile = divpath, outfile = divout, boots = bootstraps, stat = stats,  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)





### Prior Notes ###
#Genepop
test_haps <- read.genepop("ima/haps_simONLY_0126.gen", ncode=3)
#install.packages('diveRsity')
library(diveRsity)
divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest", boots = 0, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest", boots = 100, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#Solidissima
divMigrate(infile = "ima/solonly_SvC_genfromSNPs1.txt", outfile = "ima/snp_SvC_gen_divtest", boots = 100, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)


#Similis
divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_2", boots = 2000, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_tg", boots = 10, stat = "gst",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_5g", boots = 5000, stat = "g",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_3", boots = 4000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_6", boots = 10000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)


divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_6b", boots = 10000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "ima/haps_simONLY_0126.gen", outfile = "ima/haps_simONLY_0126_divtest_6c", boots = 10000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)


test_haps_SNE <- read.genepop("ima/haps_simONLY_222_SNE.gen", ncode=3)
#install.packages('diveRsity')
divMigrate(infile = "ima/haps_simONLY_222_SNE.gen", outfile = "ima/haps_222_SNE_divtest_1", boots = 10, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

divMigrate(infile = "ima/haps_simONLY_222_SNE.gen", outfile = "ima/haps_222_SNE_divtest_4", boots = 1000, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

divMigrate(infile = "ima/haps_simONLY_222_SNE.gen", outfile = "ima/haps_222_SNE_divtest_20", boots = 2000, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)




#Date: Mar 22
#Solidissima from all SNPs (specific rather than relative to similis)
#install.packages('diveRsity')
library(diveRsity)

#Just A vs B
divMigrate(infile = "Solidissima/haps_sol_SvC_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_SvC_divtest_b0", boots = 0, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_SvC_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_SvC_divtest_b400d", boots = 400, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_SvC_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_SvC_divtest_b1000all", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#Across Regions
divMigrate(infile = "Solidissima/haps_sol_region_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_region_divtest_b0", boots = 0, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_region_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_region_divtest_b400d", boots = 400, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_region_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_region_divtest_b1000all", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_region_0314.gen", outfile = "Solidissima/divMigrate22/haps_sol_region_divtest_b1000g", boots = 1000, stat = "g",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)


#Date: Apr 1

#Within Genotypes Across Regions
#A
divMigrate(infile = "Solidissima/haps_sol_A_region_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_A_region_divtest_b0", boots = 0, stat = "d",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

divMigrate(infile = "Solidissima/haps_sol_A_region_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_A_region_divtest_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_A_region_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_A_region_divtest_gb1000", boots = 1000, stat = "g",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#B (with SCC + NAN = SMA)
divMigrate(infile = "Solidissima/haps_sol_B_region_SMA_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_B_region_SMA_divtest_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
divMigrate(infile = "Solidissima/haps_sol_B_region_SMA_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_B_region_SMA_divtest_gb1000", boots = 1000, stat = "g",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#B with (SCC + SLI = Sound)
divMigrate(infile = "Solidissima/haps_sol_B_region_sound_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_B_region_sound_divtest_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)
#Same results

#B with all separate
divMigrate(infile = "Solidissima/haps_sol_B_region_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_B_region_separate_divtest_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#B with no 1999 (SMA)
divMigrate(infile = "Solidissima/haps_sol_B_region_no1999_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_B_region_no1999_divtest_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#B with ONLY 1999
divMigrate(infile = "Solidissima/haps_sol_1999only_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_B_1999onlyR_divtest_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

#A with no CCB
divMigrate(infile = "Solidissima/haps_sol_A_region_noCCB_0401.gen", outfile = "Solidissima/divMigrate0104/haps_sol_A_region_divtest_noCCB_b1000", boots = 1000, stat = "all",  filter_threshold = 0, plot_network = TRUE,  plot_col = "darkblue", para = TRUE)

