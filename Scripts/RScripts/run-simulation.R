##drift simulation script 
#part to run from bash based on the sim.R script
#09.02.2024


conda activate r_env
R

```R
library(poolSeq)
#make one file for each pop

subSync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/subset.sync", polarization="reference", 
  gen=c(1,3,1,2,1,4,1,3),repl=c(1,1,2,2,3,3,4,4)) #import only the start/end points of interest

t1Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2009.sync", polarization="reference", 
  gen=c(2,4),repl=c(1,2)) #import only the start/end points of interest

t2Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2011.sync", polarization="reference", 
  gen=c(1,2),repl=c(1,2)) #import only the start/end points of interest

t3Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2015.sync", polarization="reference", 
  gen=c(1,4),repl=c(1,2)) #import only the start/end points of interest

t4Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2022.sync", polarization="reference", 
  gen=c(1,3),repl=c(1,2)) #import only the start/end points of interest

#estimate Ne
af <- as.data.frame(subSync@alleles)
estNe_2009 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="JR.planII")
estNe_2011 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=3, poolSize=c(50, 50), method="JR.planII")
estNe_2015 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="JR.planII")
estNe_2022 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="JR.planII")

listSync <- list(subSync, t1Sync, t2Sync, t3Sync, t4Sync)
names(listSync) <- c("subSync", "t1Sync", "t2Sync", "t3Sync", "t4Sync")


# run simulation based on the starting allele frequency of each year (check this correlation)
## Ne is based on above

af_2009 <- af(t1Sync,,, 1, 2)
af_2011 <- af(t2Sync,,, 1, 1)
af_2015 <- af(t3Sync,,, 1, 1)
af_2022 <- af(t4Sync,,, 1, 1)

startaf_2009 <- data.frame(snp = names(af_2009), af = as.numeric(af_2009))
startaf_2011 <- data.frame(snp = names(af_2011), af = as.numeric(af_2011))
startaf_2015 <- data.frame(snp = names(af_2015), af = as.numeric(af_2015))
startaf_2022 <- data.frame(snp = names(af_2022), af = as.numeric(af_2022))

#same for coverage
cov_2009 <- t1Sync@alleles[,c(8,10)]
cov_2011 <- t2Sync@alleles[,c(8,10)]
cov_2015 <- t3Sync@alleles[,c(8,10)]
cov_2022 <- t4Sync@alleles[,c(8,10)]

#########################################
###now simulate the data for each year###
#########################################

#####################
####Simulation 1#####
#####################
########2009#########

file1 <- t1Sync@alleles$chr
file2 <- t2Sync@alleles$pos
start_af <- startaf_2009$af
ne_in <- round(estNe_2009)
f2_cov <- cov_2009[,1]
f4_cov <- cov_2009[,2]
samp_size <- 50

n=5
for (i in 1:n){
  outt <- simulate_allele_freq(startingAF= start_af, effective_pop=round(ne_in), T0=0, T_N=2, samp_size=50,
                          F0_depth = round(f2_cov), FN_depth = round(f4_cov))
  write.table(outt, file=paste('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/simulated.2009.',i, '.sync',sep=""), sep = '\t',
                            col.names = FALSE, row.names = FALSE, quote=FALSE)

if(i%%4 ==0){print(i)}
}

#####################
########2011#########

file1 <- t2Sync@alleles$chr
file2 <- t2Sync@alleles$pos
start_af <- startaf_2011$af
ne_in <- round(estNe_2011)
f2_cov <- cov_2011[,1]
f4_cov <- cov_2011[,2]
samp_size <- 50

n=10
for (i in 1:n){
  outt <- simulate_allele_freq(startingAF= start_af, effective_pop=round(ne_in), T0=0, T_N=1, samp_size=50,
                          F0_depth = round(f2_cov), FN_depth = round(f4_cov))
  write.table(outt, file=paste('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/simulated.2011.',i, '.sync',sep=""), sep = '\t',
                            col.names = FALSE, row.names = FALSE, quote=FALSE)

if(i%%4 ==0){print(i)}
}

#####################
########2015#########

file1 <- t23Sync@alleles$chr
file2 <- t3Sync@alleles$pos
start_af <- startaf_2015$af
ne_in <- round(estNe_2015)
f2_cov <- cov_2015[,1]
f4_cov <- cov_2015[,2]
samp_size <- 50

n=10
for (i in 1:n){
  outt <- simulate_allele_freq(startingAF= start_af, effective_pop=round(ne_in), T0=0, T_N=4, samp_size=50,
                          F0_depth = round(f2_cov), FN_depth = round(f4_cov))
  write.table(outt, file=paste('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/simulated.2015.',i, '.sync',sep=""), sep = '\t',
                            col.names = FALSE, row.names = FALSE, quote=FALSE)

if(i%%4 ==0){print(i)}
}

#####################
########2022#########

file1 <- t4Sync@alleles$chr
file2 <- t4Sync@alleles$pos
start_af <- startaf_2022$af
ne_in <- round(estNe_2022)
f2_cov <- cov_2022[,1]
f4_cov <- cov_2022[,2]
samp_size <- 50

n=10
for (i in 1:n){
  outt <- simulate_allele_freq(startingAF= start_af, effective_pop=round(ne_in), T0=0, T_N=3, samp_size=50,
                          F0_depth = round(f2_cov), FN_depth = round(f4_cov))
  write.table(outt, file=paste('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/simulated.2022.',i, '.sync',sep=""), sep = '\t',
                            col.names = FALSE, row.names = FALSE, quote=FALSE)

if(i%%4 ==0){print(i)}
}