# A script to calculate the alpha values for Frescalo 
# The estimation of trends in species occurence has not been implemented
#
# This has been checked against Frescalo output and it agrees
#
# Written by Jon Yearsley (jon.yearsley@ucd.ie) 2nd Sept 2016, modified by Oli Pescott (olipes@ceh.ac.uk), April 2024

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
cl <- makeCluster(3)
registerDoParallel(cl, cores = NULL)

R_star = 0.2703
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
  frescalo(data_in = dSplit[[i]], speciesList, spLocations, speciesNames, Phi, R_star=R_star, missing.data = 2)
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
  #trend.out <- foreach(i=1:length(sSplit2), .inorder=T, .combine='rbind') %dopar% {
  # Amend to return more than one list
  trend.out <- foreach(i=1:length(sSplit2), .inorder=T, .combine='cfunTrend', .multicombine=TRUE) %dopar% {
    trend(s_data = sSplit2[[i]], output$freq.out)
  }
}
head(trend.out, n = 10)
head(trend.out$trend.out, n = 10)
site.time.out <- trend.out$site.time.out
head(site.time.out, n = 10)
#View(trend.out)

################################################################

# If we were saving output:
# if (trend_analysis) {
#   write.table(format(trend.out$trend.out[order(trend.out$trend.out$species,trend.out$time),], digits=4, zero.print=T, width=10, scientific=F,   justify='left'), 
#               file=paste(outputPrefix,'_frescalo_trend.txt',sep=''), col.names=T, row.names=F, quote=F, sep=' ')
# }
stopCluster(cl)
gc()

## Check against fortran outputs
library(reshape2)
load(file = "data/unicorn_TF.rda")

## Site-based stuff ##
# Alpha
#par(mfrow=c(1,1))
head(unicorn_TF$stat)
alphaDF <- data.frame(alpha = output$frescalo.out$alpha, Location = output$frescalo.out$location)
alphaCompare <- merge(unicorn_TF$stat, alphaDF, by = c("Location"))
head(alphaCompare)
cor(alphaCompare$alpha, alphaCompare$Alpha) #0.9998
plot(alphaCompare$alpha, alphaCompare$Alpha, main = "Alpha") #
#plot(alphaCompare$alpha, alphaCompare$Alpha, main = "Alpha", xlim = c(15,30), ylim = c(15,30)) # 0.993
abline(a = 0, b = 1)

# Estimated species richness
#par(mfrow=c(1,1))
head(unicorn_TF$stat)
spnumDF <- data.frame(spnum = output$frescalo.out$spnum_out, Location = output$frescalo.out$location)
spnumCompare <- merge(unicorn_TF$stat, spnumDF, by = c("Location"))
head(spnumCompare)
cor(spnumCompare$spnum, spnumCompare$Spnum_out) # 0.9999
plot(spnumCompare$spnum, spnumCompare$Spnum_out, main = "Predicted \n site richness")
abline(a = 0, b = 1)

## Species x site based stuff ##
# Local frequency
#par(mfrow=c(1,1))
head(unicorn_TF$freq)
fo <- output$freq.out[,c(1,2,4)]
names(fo) <- c("Location", "Species", "freq")
lwfCompare <- merge(unicorn_TF$freq, fo, by = c("Species", "Location"))
head(lwfCompare)
cor(lwfCompare$freq, lwfCompare$Freq) # 0.9997
plot(lwfCompare$freq, lwfCompare$Freq, main = "Raw species' weighted \n freqs")
abline(a = 0, b = 1)

# Ranks
#par(mfrow=c(1,1))
ro <- output$freq.out[,c(1,2,7)]
names(ro) <- c("Location", "Species", "rank1")
lwfRank1Compare <- merge(unicorn_TF$freq, ro, by = c("Species", "Location"))
head(lwfRank1Compare)
cor(lwfRank1Compare$rank1, lwfRank1Compare$Rank1) # 0.9987
plot(lwfRank1Compare$rank1, lwfRank1Compare$Rank1, main = "Species ranks") # inflated ranks presumablt reflected in rescaled rank below
abline(a = 0, b = 1)
# Discrepancies
lwfRank1Compare$diff <- abs(round(lwfRank1Compare$Rank1, digits = 1) - round(lwfRank1Compare$rank1, digits = 1))
lwfRank1Compare[order(lwfRank1Compare$diff, decreasing = T),]
sum(lwfRank1Compare$diff==0)/nrow(lwfRank1Compare)*100 # 93.7% agreement to 1 decimal place
# Seem to be some cases of switching, e.g. Species 47 and 57 in NH53?
chk1 <- lwfRank1Compare[lwfRank1Compare$Location=="NH53",]
chk1 <- chk1[order(chk1$Rank),] ## looks like it is to do with tied frequencies
# Check another with a relatively large discrepancy in rank1 (0.7)
chk <- lwfRank1Compare[lwfRank1Compare$Location=="NT51",]
chk <- chk[order(chk$Rank),] # disagreements are mainly around where species are absent and therefore presumably assigned ranks based on ordering of names only (no effect on benchmarks)

# Rescaled freq
#par(mfrow=c(1,1))
rf <- output$freq.out[,c(1,2,5)]
names(rf) <- c("Location", "Species", "freq_1")
jDat2DFCompare <- merge(unicorn_TF$freq, rf, by = c("Species", "Location"))
head(jDat2DFCompare)
cor(jDat2DFCompare$freq_1, jDat2DFCompare$Freq1) # 0.9994
plot(jDat2DFCompare$freq_1, jDat2DFCompare$Freq1, main = "Rescaled species' \n weighted frequencies (f_ij)")
abline(a = 0, b = 1)

## Time factor and SD
fortranTrend <- unicorn_TF$trend
names(fortranTrend)[1:2] <- c("species","time")
para.trend <- trend.out$trend.out
para.trend$time <- ifelse(para.trend$time == 1, 1984.5, 1994.5)
tfCompare <- merge(fortranTrend, para.trend, by = c("time", "species"))
# Time factors
cor(tfCompare$tFactor, tfCompare$TFactor) # 0.998
# Time factors SDs
cor(tfCompare$StDev.x, tfCompare$StDev.y) # 0.999
par(mfrow=c(1,2))
plot(tfCompare$tFactor, tfCompare$TFactor, main = "Time factors") # 
abline(a = 0, b = 1)
plot(tfCompare$StDev.x, tfCompare$StDev.y, main = "Time factor SDs") # 
abline(a = 0, b = 1)

## Site x time period x species predicted occupancies
# Following Biljsma (2013), this part takes the rescaled species frequencies (freq1) from the frescalo() function, and the 
# mean species time factor for a time period from the trend() function, and uses them to predict
# species occupancy under standard recording effort (i.e. when s_it = 1)

head(output$freq.out)
head(output$freq.out[output$freq.out$species=="Species 1",])
head(para.trend)
list_data2 <- by(data = output$freq.out, as.factor(output$freq.out$species), function(x) createPijt(x, trend = para.trend))
options(scipen = 999)
occ.out <- do.call(rbind, list_data2)
row.names(occ.out) <- 1:nrow(occ.out)
occ.out$P_ijt <- round(occ.out$P_ijt, digits = 3)
hist(occ.out$P_ijt)

