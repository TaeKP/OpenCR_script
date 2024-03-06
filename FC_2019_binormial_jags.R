# Spatial capture-recapture for fishing cat population density
# Khao Sam Roi Yot dataset 2019
# Bayesian analysis with JAGS
# Binomial detection model with habitat patch
# Author: TaeKP
# 5 March 2024
# ----------------------------------------------------------------------------

library(jagsUI)
library(mcmcOutput)

# Get the detection data and check
# --------------------------------
FC_raw <- read.csv("FC_2019_jags.csv", comment="#")
head(FC_raw)
fc2019 <- as.matrix(FC_raw[, 1:50])
dim(fc2019)  # 34 x 50
table(fc2019)  
#    0    1    2    3    4    5    6    7    8    9   12   16 
# 1651   20   10    2    4    2    1    4    2    2    1    1 

# Get the trap locations
# ----------------------
traps <- read.csv("FCTraps_2019.csv", comment="#")
head(traps)  # rows for 50; coordinates already in km, we'll center as well

Detlocs <- as.matrix(scale(traps, scale=FALSE))
plot(Detlocs)
( nDetlocs <- nrow(Detlocs) )

# Decide on state space
# ---------------------
# allow a buffer around the detector array

# make habitat mask
#library(raster)
#library(makeJAGSmask)
#?makeJAGSmask
#require(raster)
#KSRY.raster <- raster('KSRY_hab_raster')
#KSRY.hab <- convertRaster(KSRY_hab_raster, Detlocs, plot = TRUE)

buffer <- 8  # kilometres

# outer edges of the state space
( xmin <- min(Detlocs[,1]) - buffer )
( xmax <- max(Detlocs[,1]) + buffer )
( ymin <- min(Detlocs[,2]) - buffer )
( ymax <- max(Detlocs[,2]) + buffer )

xlim <- c(xmin, xmax)                                       # x limits of state-space
ylim <- c(ymin, ymax)                                       # y limits of state-space
(A <- (max(xlim) - min(xlim)) * (max(ylim) - min(ylim)))      # State-space area
#1744

# Plot it:
plot(Detlocs, xlim=xlim, ylim=ylim,
     bty='n', pch=16, col='black')
rect(xmin, ymin, xmax, ymax, border='red', lwd=3)

# All the detectors were active for 75 days
nOcc <- 75

# Augment and bundle the data
# ---------------------------
( nCaps <- nrow(fc2019) )
nAug <- 200  # number of rows to add  #looks good at 200 inds
yAug <- rbind(fc2019, matrix(0, nAug, nDetlocs))
( M <- nrow(yAug) )

# Bundle data for JAGS
jagsData <- list(y = yAug, M = M, nOcc = 75, A=A,
                 w = c(rep(1, nCaps), rep(NA, M-nCaps)),
                 Detlocs = Detlocs, nDetlocs = nDetlocs,
                 xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
# CRN: "Glue together" data from all 3 years into arrays of the correct dimension
# (see in jags code below)

str(jagsData)

# List of 10
# $ y       : num [1:x, 1:50, 1:3] 0 0 0 0 0 0 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:50] "trap01" "trap02" "trap03" "trap04" ...
# $ M       : num [1:3]
# $ nOcc    : num [1:3]
# $ w       : num [1:x, ] 1 1 1 1 1 1 1 1 1 1 ...
# $ Detlocs : num [1:50, 1:2] -6.25 -7.79 -5.63 -3.36 -1.35 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:2] "x" "y"
# ..- attr(*, "scaled:center")= Named num [1:2] 600 1347
# .. ..- attr(*, "names")= chr [1:2] "x" "y"
# $ nDetlocs: int 50
# $ xmin    : num -17.8
# $ xmax    : num 17.1
# $ ymin    : num -31.1
# $ ymax    : num 29.5

# Write JAGS model file
cat(file="SCRfc_yrLoop_binom.jags", "
model{
  
  for(t in 1:Tmax){
  
    # Priors
    p0[t] ~ dbeta(1, 1)          # Baseline detection probability
    sigma[t] ~ dunif(0, 3)       # Half-normal scale
    omega[t] ~ dbeta(1, 1)       # Data augmentation parameter
  
    # Likelihood
    for(i in 1:M[t]){             # Loop over all M individuals
    
      w[i,t] ~ dbern(omega[t])      # w = 1 if animal is real/present
      
      # State model: point process model
      AC[i, 1, t] ~ dunif(xmin[t], xmax[t]) # x-coord of activity centre
      AC[i, 2, t] ~ dunif(ymin[t], ymax[t]) # y coord of activity centre
      
      # Observation model: p ~ distance between trap and estimated AC
      for(j in 1:nDetlocs[t]){           # Loop over all detectors
        d2[i, j, t] <- (AC[i, 1, t] - Detlocs[j, 1, t])^2 +
        (AC[i, 2, t] - Detlocs[j, 2, t])^2           # distance^2
      p[i, j, t] <- p0[t] * exp(- d2[i, j, t]/(2*sigma[t]^2)) # Detection prob
      y[i, j, t] ~ dbin(p[i, j, t] * w[i, t], nOcc)     # The observed data
      }
    }

    # Derived quantities
    N[t] <- sum(w[t])                       # Population size in state-space (=area)
    D[t] <- N[t] / A[t]                        # Density over state-space
  }
}
")

# Run the model
# -------------
# Parameters to save
parameters <- c("p0", "sigma", "omega", "N", "D", "AC", "w")

# MCMC settings
ni <- 1000; nb <- 50; nc <- 3; nt <- 5; na <- 1000

#load("FC2019_jagsUI.RData")
out <- jags(jagsData, NULL, parameters, "SCRfc_yrLoop_binom.jags", n.iter=ni, n.burnin=nb, n.chains=nc,
              n.thin=nt, n.adapt=na, DIC=FALSE, parallel=TRUE)

#save(out, file="FC2019_jagsUI.RData")
# takes 7 mins.

# check convergence and summarize posteriors 
traceplot(out) # Not shown
print(out, 3)

# Check the output
# ----------------
library(mcmcOutput)
( mco <- mcmcOutput(out) )

diagPlot(mco, params=1:4)
plot(mco, 1:4)
View(summary(mco))

# check that N << M
max(mco$N)
M

#plot population with HDI
N <- mco$N
postPlot(N, xlab="Population", showCurve=TRUE)
