##########################
# Function definitions for frescalo_parallel R Script
#
# This script is sourced in the frescalo_parallel.R script
#
# Jon Yearsley (jon.yearsley@ucd.ie) Oct 2016
# All these scripts have been checked against the original FRESCALO FORTRAN program
#
# Minor amendments by Oli Pescott (olipes@ceh.ac.uk) Apr 2024, to produce time factor SDs, and to match Hill outputs for sparse datasets
# Checks were made against both fortran outputs and independent R coding of fortran approach here: https://github.com/sacrevert/fRescalo
# Unicorns dataset from sparta used to check all outputs (see https://github.com/biologicalRecordsCentre/sparta) and https://github.com/sacrevert/fRescalo 
#
##########################

require(foreach, quietly=T)

min.fun <- function(alpha,fij,Phi) {
  # The function to minimise to fit Frescalo model to the data 
  # Page 5 column 2 in Hill (2011)
  F=1-(1-fij)^alpha
  #F[abs(fij-1)<1e-10] = 1 # Yearsley
  F[abs(fij-1)<1.0E-10] = 0.99999 # to match Hill (2012), although not needed for Yearlsey calc strictly
  F[abs(fij-1)>0.99999] = 1.0E-10
  return(sum(F^2)/sum(F) - Phi)
}

speciesListFun = function(spList, species) {
  # A function that returns the presences and absences of species per region
  # spList should be a data frame with columns location, species and time  

  # Create list of all locations
  locationList = unique(spList$location)
  
  # Initialize output list
  out = lapply(c(1:length(locationList)), function(i) {i})
  
  for (i in 1:length(locationList)) {
    spListSub = subset(spList, location==locationList[i])
    out[[i]] = as.integer(species %in% spListSub$species)
  }
  return(out)
}

frescalo = function(data_in, speciesList, spLocations, speciesNames, Phi=0.74, R_star=0.2703, missing.data = 2) {
  ## testing
  #data_in = dSplit[[1]]; Phi=0.74; R_star=0.27; missing.data = 1
  #missing.data = 2
  
  # This function calculates the sampling effort multiplier, alpha, that would equate sampling effort across all regions
  # This is the main method of Frescalo.
  
  # The script returns a list containing two data frames: frescalo.out, freq.out
  # frescalo.out contains the recorder-effort multiplier for each focal location, and other related information
  # freq.out contains the corrected local frequency of each species in each nehgbourhood of a focal region, 
  # along with other related information
  # Species identified as a benchmark in freq.out (benchmark=1) can be used to assess trends
  
  # Argument missing.data defines what to do if data in the neighbourhood are missing.
  # missing.data = 1 do not proceed with the calculation 
  # missing.data = 2 (default) set missing species data to all absences and proceed
  
  # R_star defines the corrected frequency threshold for benchmark species 
  # (default = 0.27, means that after correction for recorder effort the top 27% of species are used as benchmarks)
  
  location1List = unique(data_in$location1)
  nSpecies = sapply(speciesList[match(location1List,spLocations,nomatch=length(spLocations))], FUN=sum, simplify=T)
  out = data.frame(location=location1List, 
                   nSpecies=nSpecies, 
                   phi_in=NA, alpha=NA, phi_out=Phi, spnum_in=NA, spnum_out=NA, iter=NA)
  
  freq.out = data.frame(location=c(), species=c(), pres=c(), freq=c(), rank=c(), rank_1=c(), benchmark=c())
  
  for (f in 1:length(location1List)) {
    focal = location1List[f]
    focal_d = subset(data_in, location1==location1List[f])
    
    # Identify species in neighbourhood of focal region
    # If data are mssing assign it to the empty species list (last element of speciesList)
    neighbourhood = match(focal_d$location2, spLocations, nomatch=length(speciesList)+1) # empty neighbourhood sites are NULL
    speciesRegional = speciesList[neighbourhood]
    speciesRegional <- lapply(speciesRegional, function(x) if(is.null(x)) rep(0, times = length(speciesNames)) else x) # correct issue with last neighbourhood
    
    missingData <- neighbourhood==length(speciesList)+1 # length(speciesList) is number of locations; missing data is neighbourhoods with
    if (any(missingData) & missing.data==1) {
      # Species data missing from a location in the neighbourhood. Ignore this focal location
      warning(paste('Removing location ',focal,'. Missing data for locations in the neighbourhood.', sep=''))
      
      out$alpha[f] = NA
      out$iter[f] = NA
      out$phi_in[f]= NA 
    } else {
      # Calculate weights of locations in the neighbourhood
      weights = focal_d$w/(sum(focal_d$w)+1.0E-10)
      
      # Create weighted neighbourhood frequencies (checked against Frescalo)
      frequency = Reduce('+',Map('*',as.list(weights), speciesRegional))  
      # Added by Oli Pescott, Apr 2024, to avoid zeros and match Hill (2012) outputs
      frequency <- ifelse(frequency==0, frequency+1.0E-10, frequency)
      phi_in = sum(frequency^2) / sum(frequency)
      
      # Calculate the multiplier (alpha) that equalises recording effort 
      alpha.min = 1  # Minimum alpha ( =1 means no correction required)
      alpha.max = 5
      # Increase alpha.max until min.fun() becomes positive (i.e. ensure there is a zero)
      while (min.fun(alpha.max, frequency, Phi)<0) { alpha.max = alpha.max + 5}  
      while (min.fun(alpha.min, frequency, Phi)>0) { alpha.min = alpha.min/2}  
      
      # Find sampling-effort multiplier
      sol=uniroot(min.fun,interval=c(alpha.min,alpha.max), tol=0.0003, frequency, Phi)
      out$alpha[f] = sol$root
      out$phi_in[f] = phi_in
      out$iter[f] = sol$iter
      # Expected species richness before recorder effort correct
      out$spnum_in[f] = sum(frequency)   
      # Expected species richness after recorder effort correct
      out$spnum_out[f] = sum(1-exp(sol$root*log(1-frequency)))
      
      # Create the data frame with the local species frequencies after correction
      freq.order = order(frequency, decreasing=T) # no explicit method to deal with ties, so any ties just arranged using original name order
      focal.ind = match(focal,spLocations,nomatch=length(spLocations))
      
      # Pick out benchmark species (assumes that species are ordered by rank)
      benchmarkSpecies = rep(0, times=length(frequency))
      R_prime = c(1:length(speciesList[[focal.ind]]))/out$spnum_out[f] # rescaled rank
      # Benchmark species are either those where rescaled rank < R_star, or are ranked number 1 even tho all R_prime > R_star
      benchmarkSpecies[R_prime<R_star | c(1:length(speciesList[[focal.ind]]))==1 ] = 1
      #benchmarkSpecies[R_prime<R_star] = 1 ## Yearsley original only had one condition. Adding second increases correlation with Hill fortran marginally
      
      # freq_1 is the corrected neighbourhood frequencies
      freq.out = rbind(freq.out, data.frame(location=focal, 
                                            species=speciesNames[freq.order], 
                                            pres=speciesList[[focal.ind]][freq.order],
                                            freq=frequency[freq.order],
                                            freq_1=1-exp(sol$root*log(1-frequency[freq.order])),
                                            rank=c(1:length(speciesList[[focal.ind]])),
                                            rank_1=R_prime,
                                            benchmark=benchmarkSpecies))
    }
    #progress::setTxtProgressBar(progress_bar, value = f)
  }
  return(list(frescalo.out=out, freq.out=freq.out))
}


trend = function(s_data, freq.out) {
  # Frescalo trend analysis
  # Function to calculate the time scaling factor for each species in a single year bin
  #
  # s_data are the records from a single time period. Each row is a reocrd of a species. 
  # The data frame needs at least 3 columns:
  #      location  an identifier for the location of the record
  #      species   an identifier for the species
  #      time      an identifier for the time bin
  #
  # freq.out is a data frame with the corrected neighbourhood frequencies from a FRESCALO analysis (see the frescalo function). 
  # This data frame should have columns:
  #      location  an identifier for the location
  #      species   an identifier for the species
  #      freq_1    the corrected neighbourhood frequency
  #      rank_1    the rank of the corrected neighbourhood frequency
  #      benchmark =1 if this species should be used as a benchmark species
  
  # Check that the input data (s_data) is from a single time bin
  timeBin = unique(s_data$time)
  if (length(timeBin)>1) {
    warning('More than one time bin supplied to trend() function')
  }
  #freq.out <- output$freq.out
  #sSplit2 = split(s, as.factor(s$time))
  #s_data = sSplit2[[1]]
  #Species 5 in time period has zero TF, check how I cna keep this
  locationList = as.character(unique(freq.out$location))
  #spList = unique(s_data$species)
  # OLP: use full list of species so that time-period specific zeros are kept
  spList = unique(s$species)
  
  # Calculate the proportion of benchmark species in each hectad (for this time bin)
  focal_s = split(s_data, factor(s_data$location, levels=locationList))
  focal_bench = split(freq.out, factor(freq.out$location, levels=locationList))
  s_it = mapply(FUN=function(x,y) {sum(x$benchmark[x$species %in% y$species])/sum(x$benchmark)}, 
                x=focal_bench, y=focal_s, SIMPLIFY=T, USE.NAMES=F)
  s_it[is.na(s_it)] <- 1.0E-7 # added by OL Pescott Apr 2024, as apparently sometimes x$benchmark can be zero (and so s_it is NaN)
  # Calculate weights to downweight infrequenctly sampled locations 
  # i.e. seeing fewer than 9.95% of the benchmark species (s_it<0.0995)
  w = rep(1, each=length(locationList))
  w[s_it<0.0995] = 10*s_it[s_it<0.0995]+0.005 
  # Add site x time period recording effort info to output
  site.time.out <- data.frame(location = locationList, 
                                          time = rep(unique(s_data$time), times = length(locationList)), 
                                          s_it = s_it, 
                                          w = w)
  
  focal_s2 = split(s_data, factor(s_data$species, levels=spList))
  sumP_ijtw =  unlist(lapply(FUN=function(X){sum(w[locationList%in%X$location])}, 
                     X=focal_s2), use.names=F)

  focal_f = split(freq.out, factor(freq.out$species, levels=spList))
  sf = lapply(FUN=function(x) {x$freq_1[match(locationList,x$location)]*s_it}, focal_f)

  # Calculate Q_ijt. If sf>0.98 Set Q_ijt=-log(1-0.98)=3.912023
  # P_ijt = 1-exp(-Q_ijt x_jt) where x_jt is the time factor for species j at time t
  # and P_ijt is prob of observing species j in hectad i at time t
  Q_ijt = lapply(FUN=function(x){y=-log1p(-x); y[y>3.912023]=3.912023; return(y)}, X=sf)

  StDev = x = xSD = rep(NA, times=length(spList))   # Vector to contain the time factors
  sptot1 = estvar = rep(0, times=length(spList))   # Vector to contain the time factor standard deviations
  for (i in 1:length(spList)) {
    # Rescale frequencies by effort for all hectads
    # tmp = subset(freq.out, species==spList[s])
    # sf_tmp = tmp$freq_1[match(locationList,tmp$location)]*s_it

    # # Calculate Q_ijt   
    # # P_ijt = 1-exp(-Q_ijt x_jt) where x_jt is the time factor for species j at time t
    # # and P_ijt is prob of observing species j in hectad i at time t
    # Q_ijt = rep(-log(1-0.98), each=length(locationList))
    # Q_ijt[sf[[s]]<0.98] = -log(1-sf[[s]][sf[[s]]<0.98])
    
    # Select a x max that ensire a sign change in min_trend_fun
    if (any(Q_ijt[[i]]>0)) {
      xMax = 5 # Yearsley original
      #xMax = 1
      while(min_trend_fun(xMax, Q_ijt[[i]], w, sumP_ijtw[i])<0) {xMax = xMax+5}
      sol = uniroot(min_trend_fun,interval=c(0,xMax), tol=0.0005, Q_ijt[[i]], w, sumP_ijtw[i])
      x[i] = sol$root
      #x = sol$root
      #i = 52
      #####################################################
      ## Following is for calculation of time factor SDs ##
      ############ Oli Pescott, Apr 2024 ##################
      #####################################################
      estvar[i] = sum((1-exp(-Q_ijt[[i]]*x[i])) * (1-(1-exp(-Q_ijt[[i]]*x[i]))) * w * w)
      sptot1[i] = sumP_ijtw[i] + sqrt(estvar[i]) 
      xSD[i] <- sptot1[i]/(sumP_ijtw[i]+0.0000001)
      while(min_trend_fun(xSD[i], Q_ijt[[i]], w, sptot1[i])<0) {xSD[i]  = xSD[i] + 1}
      # Small value added to xSD[i] to accommodate species with no data in a time bin (stops uniroot failing)
      sol2 = uniroot(min_trend_fun,interval=c(0,xSD[i]+1.0E-10), tol=0.0005, Q_ijt[[i]], w, sptot1[i])
      xSD[i] = sol2$root
      StDev[i] = abs(xSD[i] - x[i])
    }
  }
  df <- data.frame(species=spList, time=timeBin, tFactor=x, StDev = StDev, estvar=estvar, sptot1=sptot1)
  df <- df[order(df$species, df$time),]
  return(list(trend.out = df, site.time.out = site.time.out ))
}

min_trend_fun = function(x, Q, w, sumPw) {
  # A function to solve sum_i P_ijt * w_i = sum_i (1-exp(-Q_ijt*x))*w_i in the Frescalo trend analysis Hill (2012)
  # P_ijt = 1-exp(-Q_ijt*x)
  # w_i is a weighting factor to downweight area with low fraction of benchmark species (s_it<0.095)
  # P_ijt = probability of observing species j in region i at time t
  return(sum((1-exp(-Q*x))*w) - sumPw)
}

cfun = function(...) {
  # Bespoke function to combine the output from the frescalo() function
  input_list <- list(...)
  return(list(frescalo.out=Reduce('rbind',Map(function(x){x[[1]]},input_list)), 
              freq.out=Reduce('rbind',Map(function(x){x[[2]]},input_list))))
}

cfunTrend = function(...) {
  # Another bespoke function to combine the output from the trend() function
  input_list <- list(...)
  return(list(trend.out=Reduce('rbind',Map(function(x){x[[1]]},input_list)), 
              site.time.out=Reduce('rbind',Map(function(x){x[[2]]},input_list))))
}

# Addtional function for creating P_ijt values
createPijt <- function(x, trend){
  list_data <- by(data = x, as.factor(x$location), function(x) (1-(1-x$freq_1)^trend[trend$species==x$species,]$tFactor))
  # Initialize variables
  species <- vector("list", length(list_data))
  IDs <- vector("list", length(list_data))
  Values <- vector("list", length(list_data))
  Labels <- vector("list", length(list_data))
  labels <- unique(trend$time) # Label vector
  # Populate the variables with data from the list
  for (i in seq_along(list_data)) {
    species[[i]] <- rep(unique(x$species), times = length(list_data[[i]]))
    IDs[[i]] <- rep(names(list_data)[i], times = length(list_data[[i]]))
    Values[[i]] <- list_data[[i]]
    Labels[[i]] <- labels
  }
  # Combine all elements into vectors
  species <- unlist(species)
  IDs <- unlist(IDs)
  Values <- unlist(Values)
  Labels <- unlist(Labels)
  df <- data.frame(species = species, location = IDs, time = Labels, P_ijt = Values)
  return(df)
}  



#######################################################
# This function was written to validate the approach of FRESCALO
# trend2 = function(s_data, freq.out) {
#   # Function to calculate the scaling factor for each species in a single year bin
#   # This function uses the same algorithm as the Frescalo FORTRAN code.
#   # This agrees with the output from FRESCALO
#   
#   # s_data is species data from a single time bin
#   # freq.out are the corrected species frequencies averaged acrosss all times
#   
#   timeBin = unique(s_data$time)
#   print(paste('Time bin =', timeBin))
#   
#   if (length(timeBin)>1) {
#     warning('More than one time bin supplied to trend() function')
#   }
#   # Split up data into different species
#   s_data$species = as.factor(s_data$species)
#   spList = levels(s_data$species)
# 
#   # List all the focal locations
#   locationList = unique(freq.out$location)
#   
#   
#   # Calculate the proportion of benchmark species in each hectad (for this time bin)
#   s_it = rep(0, times=length(locationList))
#   for (f in 1:length(locationList)) {
#     focal_s = s_data$species[s_data$location==locationList[f]]
#     focal_bench = freq.out$species[freq.out$location==locationList[f] & freq.out$benchmark==1]
#     
#     s_it[f] = sum(focal_bench %in% focal_s) / length(focal_bench)
#   }
#   
#   # Calculate weights to downweight infrequenctly sampled locations (s_it<0.1)
#   wgt = rep(1, each=length(locationList))
#   wgt[s_it<0.0995] = 10*s_it[s_it<0.0995]+0.005
#   
#   x = rep(NA, times=length(spList))   # Vector to contain the time factors
#   
#   for (s in 1:length(spList)) { # Loop through every species
#     kmax=100
#     tf=1
#     loop=TRUE
#     k=1
#     iocc = locationList %in% s_data$location[s_data$species==spList[s]] # =1 if species s is found at location i
#     
#     
#     while (k<=kmax & loop) {
#       esttot=0
#       sptot=0
#       jtot=0
#       ic1=0
#       ic2=0
#       
#       # Rescale frequencies by effort for all hectads
#       # pfac = f_ij * s_it
#       focal_f = subset(freq.out, species==spList[s])
#       pfac = focal_f$freq_1[match(locationList,focal_f$location)]*s_it
#       for (i in 1:length(locationList)) {
#         
#         if (pfac[i]>0) { ic1=ic1+1}
#         if (pfac[i]>0.98) {
#           pfac[i]=0.98
#           ic2=ic2+1
#         }
#         plog=-log(1-pfac[i])
#         estval=1-exp(-plog*tf)
#         esttot=esttot+wgt[i]*estval
#         sptot=sptot+wgt[i]*iocc[i]
#         jtot=jtot+iocc[i]
#       }
#       
#       
#       k=k+1
#       loop = abs(sptot-esttot)>=0.0005
#       if (loop) {
#         tf=tf*sptot/(esttot+0.0000001)
#       }
#     }
#     x[s] = tf
#   }
#   return(data.frame(species=spList, time=timeBin, tFactor=x))
# }

