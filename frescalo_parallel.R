
# A script to calculate the alpha values for Frescalo 
# The estimation of trends in species occurence has not been implemented
#
# This has been checked against Frescalo output and it agrees
#
# Written by Jon Yearsley (jon.yearsley@ucd.ie) 2nd Sept 2016

# Modified:
# 6/9/2016 (JY): Data sent to cluster nodes in chunks (variable=chunkSize)
# 8/9/2016 (JY): Save species frequency file and include expected species richness in output
# 6/10/2016 (JY): Include validated trend analysis

#setwd('/home/jon/MEGA/Jack/Frescalo')
#setwd('/home/jon/MEGA/MainFolder/Jack/Frescalo')
#rm(list=ls())  # Remove all variables from the memory
dir <- getwd()
setwd(dir)
require(foreach, quietly=T)
require(doParallel, quietly=T)
source('scripts/frescalo_functions.R')

# Set the size of the cluster (number of nodes). 
# cores=NULL will automatically pick the number of cores per node 
# (by defauult this is half the number of available cores)
#cl <- makeCluster(12)
cl <- makeCluster(1)
registerDoParallel(cl, cores = NULL)

R_star = 0.27
trend_analysis = TRUE

# weight.file = './bsbi_hectads_2000_2006_2010_weights.txt'
# species.file = './bsbi_hectads_2000_2006_2010_sample.txt'
# outputFilename = './bsbi_hectads_frescalo_out.txt'

#weight.file = './TestData/weights.txt'
#species.file = './TestData/Test.txt'

outputPrefix = 'outputs/testApr24_'

Phi = 0.74    # The standardisation for recorder effort
chunkSize = 250 # Number of hectads to pass to each compute node


# Import data
#d = read.table(weight.file, header=F, col.names=c('location1','location2','w','w1','w2','nsim','ndist'), stringsAsFactors=F)
#s = read.table(species.file, header=F, col.names=c('location','species','time'), stringsAsFactors=F)
#write.csv(dat, file = "outputs/cluster2kmData.csv")
#write.csv(final_weights, file = "outputs/cluster2kmWeights.csv")
#d <- read.csv(file = "data/cluster2kmWeights.csv", stringsAsFactors = F, header = T)
#d <- read.csv(file = "data/vcDatTetNeighs_allSites_499.csv", stringsAsFactors = F, header = T)
#s <- read.csv(file = "data/cluster2kmData.csv", stringsAsFactors = F, header = T)

# unicorns test data to compare with project 227 (also see https://github.com/sacrevert/fRescalo)
s <- read.csv(file = "data/clusterTestDat.csv", stringsAsFactors = F, header = T)
s <- s[,c(4,2,3)]
d <- read.delim(file = "data/GB_LC_Wts.txt", header = F, sep = "")
d <- d[,c("V1","V2","V3")]
names(d) <- c("location1", "location2", "w")

##############################################################
spLocations = unique(s$location)
d <- d[d$location1 %in% spLocations,] # just keep relevant neighbourhoods
speciesNames = as.character(unique(s$species))   # Create list of unique species
# The variable species could be made numerical to make it more efficient

# For each region record presence/absence of each species (a_ij in Hill 2011)
locationGroups = as.factor(rep(c(1:ceiling(length(spLocations)/chunkSize)),each=chunkSize))
sSplit = split(s, locationGroups[match(s$location, spLocations)])  # Split species data up into hectads

#idx = iter(sSplit)
#speciesList <- foreach(spList = 1:length(sSplit), .inorder=T, .combine='c') %dopar% {
speciesList <- foreach(i = 1:length(sSplit), .inorder=T, .combine='c') %dopar% {
  speciesListFun(spList = sSplit[[i]], species = speciesNames) # which species are in each location?
}
#head(speciesList[[1]])
# Add an additional species list where everything is absent -- assume this was about testing missing.data option? -- OLP, March 2023
#speciesList[[length(speciesList)+1]] = rep(0, times=length(speciesNames))
#spLocations = c(spLocations,'null_location')

#################################################################
# For each focal regional calculate the sampling effort multipler
dSub = d[,1:3] # not really needed if d already has 3 columsn
location1List = as.character(unique(dSub$location1))    # Create list of unique focal regions
location1Groups = as.factor(rep(c(1:ceiling(length(location1List)/chunkSize)),each=chunkSize))

dSplit = split(dSub, location1Groups[match(dSub$location1, location1List)])  # Split neighbourhood data up into focal regions
output <- foreach(i=1:length(dSplit), .inorder=T, .combine='cfun', .multicombine=TRUE) %dopar% {
  frescalo(data_in = dSplit[[i]], speciesList, spLocations, speciesNames, Phi, R_star=0.27, missing.data = 2)
}
head(output$frescalo.out)
head(output$freq.out)

# Write the output to a text file
write.table(format(output$frescalo.out[order(output$frescalo.out$location),], digits=4,zero.print=T, width=10), 
            file=paste(outputPrefix,'_frescalo_out.txt',sep=''), col.names=T, row.names=F, quote=F, sep=' ')

write.table(format(output$freq.out[order(output$freq.out$location, output$freq.out$rank),], digits=4, zero.print=T, width=10, scientific=F, justify='left'), 
            file=paste(outputPrefix,'_frescalo_freq.txt',sep=''), col.names=T, row.names=F, quote=F, sep=' ')

#################################################################
# Trend analysis
if (trend_analysis) {
  # Do the Frescalo trend analysis if there are more than 1 year bins (use same location groups as sSplit)
  sSplit2 = split(s, as.factor(s$time))  # Split species data up into year bins
  #idx3 = iter(sSplit2)
  #trend.out <- foreach(s_data=idx3, .inorder=T, .combine='rbind') %dopar% {
  trend.out <- foreach(i=1:length(sSplit2), .inorder=T, .combine='rbind') %dopar% {
    #trend(s_data = sSplit2[[i]], output$freq.out, calcSD = T)
    trend(s_data = sSplit2[[i]], output$freq.out)
    #trendResult <- trend(s, outTestAll$freq.out)
    #trend(s_data, output$freq.out)
  }
}
head(trend.out, n = 25)

################################################################

if (trend_analysis) {
  write.table(format(trend.out[order(trend.out$species,trend.out$time),], digits=4, zero.print=T, width=10, scientific=F, justify='left'), 
              file=paste(outputPrefix,'_frescalo_trend.txt',sep=''), col.names=T, row.names=F, quote=F, sep=' ')
}

stopCluster(cl)
gc()

## Check against fortran outputs
load(file = "data/unicorn_TF.rda")
fortranTrend <- unicorn_TF$trend
names(fortranTrend)[1:2] <- c("species","time")
trend.out$time <- ifelse(trend.out$time == 1, 1984.5, 1994.5)
## Note that the parallel frescalo currently drops species/time period combinations where the species was completely absent
tfCompare1 <- merge(fortranTrend, trend.out, 
                   by = c("time", "species"), all = T)
head(tfCompare1)
tfCompare2 <- merge(fortranTrend, trend.out, 
                   by = c("time", "species"), all.y = T) # drop NAs in x for correlation
cor(tfCompare2$tFactor, tfCompare2$TFactor) # 0.996
plot(tfCompare2$tFactor, tfCompare2$TFactor, main = "Time factors") # 
abline(a = 0, b = 1)



# source('../R_Scripts/os2eastnorth.R')
# 
# en=os2eastnorth(GR=output$frescalo.out$location, hectad=T)$en
# output$frescalo.out$x = en[,1]
# output$frescalo.out$y = en[,2]
# 
# library(ggplot2)
# 
# ggplot(data=output$frescalo.out, aes(x=x, y=y, colour=alpha)) + geom_point()
