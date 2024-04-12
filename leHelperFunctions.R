#### Helper functions ####
# iterate function across rows of dataframe
f_lapply_row <- function(df) {
  lapply(seq_len(nrow(df)), function(i) as.list(df[i,,drop=FALSE]))
}
# extract 100 fitted values (= rows) per species (49 years = columns)
extractDf <- function(x) {do.call(rbind, lapply(x, function(x) t(as.data.frame(x)[[1]]))) 
}
# fit 100 linear models to draws from mean and sd-specified distributions from Frescalo output
fitLinModSp <- function(x, timescale = timescale) { 
  t1 <- f_lapply_row(x) # convert data frame to list of lists (one list for each time period of Frescalo results)
  t2 <- lapply(t1, function(x) rnorm(100, mean = x$tFactor, sd = x$StDev))
  t3 <- data.frame(do.call(cbind, t2))
  names(t3) <- c("1", "2")
  t4 <- f_lapply_row(t3)
  t5 <- list()
  ## Extract predicted values from model fit for ALL years
  t5 <- lapply(t4, function(x, timescale = timescale) {
    lTimes <- data.frame(year = c(1,2))
    tmp <- data.frame(tfactor = unlist(x), year = lTimes)
    lMod <- lm(tmp$tfactor ~ tmp$year)
    lMod$coefficients
  })
}
## end of helper functions
##########################