# Spatial capture-recapture for fishing cat population density
# Khao Sam Roi Yot dataset 2021
# Bayesian analysis with JAGS
# Binomial detection model 
# Author: TaeKP
# 5 March 2024
# ----------------------------------------------------------------------------

library(jagsUI)
library(mcmcOutput)

# Get the detection data and check
# --------------------------------
FC_raw <- read.csv("FC_2021_jags.csv", comment="#")
head(FC_raw)
fc2021 <- as.matrix(FC_raw[, 1:50])
dim(fc2021)  # 55 x 50
table(fc2021)  
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   23   24   26 
# 2672   19   17   11    5    6    4    3    1    2    1    1    1    2    1    1    1    1    1 

# Get the trap locations
# ----------------------
traps <- read.csv("FCTraps_2021.csv", comment="#")
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

#create habitat buffer
buffer <- 8  # kilometres

# outer edges of the state space
( xmin <- min(Detlocs[,1]) - buffer )
( xmax <- max(Detlocs[,1]) + buffer )
( ymin <- min(Detlocs[,2]) - buffer )
( ymax <- max(Detlocs[,2]) + buffer )

xlim <- c(xmin, xmax)                                       # x limits of state-space
ylim <- c(ymin, ymax)                                       # y limits of state-space
(A <- (max(xlim) - min(xlim)) * (max(ylim) - min(ylim)))      # State-space area
#1743

# Plot it:
plot(Detlocs, xlim=xlim, ylim=ylim,
     bty='n', pch=16, col='black')
rect(xmin, ymin, xmax, ymax, border='red', lwd=3)

# All the detectors were active for 148 days
nOcc <- 148

# Augment and bundle the data
# ---------------------------
( nCaps <- nrow(fc2021) )
nAug <- 350  # number of rows to add  #looks good at 200 inds
yAug <- rbind(fc2021, matrix(0, nAug, nDetlocs))
( M <- nrow(yAug) )       #405

# Bundle data for JAGS
jagsData <- list(y = yAug, M = M, nOcc = nOcc, A=A,
                 w = c(rep(1, nCaps), rep(NA, M-nCaps)),
                 Detlocs = Detlocs, nDetlocs = nDetlocs,
                 xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
str(jagsData)

# List of 11
# $ y       : num [1:405, 1:50] 0 0 0 0 0 0 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:50] "trap01" "trap02" "trap03" "trap04" ...
# $ M       : int 255
# $ nOcc    : num 148
# $ A       : num 1743
# $ w       : num [1:255] 1 1 1 1 1 1 1 1 1 1 ...
# $ Detlocs : num [1:50, 1:2] -5.65 -8.05 -6.21 -3.38 -1.37 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:2] "x" "y"
# ..- attr(*, "scaled:center")= Named num [1:2] 600 1348
# .. ..- attr(*, "names")= chr [1:2] "x" "y"
# $ nDetlocs: int 50
# $ xmin    : num -16.1
# $ xmax    : num 15
# $ ymin    : num -28.7
# $ ymax    : num 27.4

# Write JAGS model file
cat(file="SCRfc_binorm.jags", "
model{
   # Priors
  p0 ~ dbeta(1, 1)          # Baseline detection probability
  sigma ~ dunif(0, 3)       # Half-normal scale
  omega ~ dbeta(1, 1)       # Data augmentation parameter
  
  # Likelihood
  for(i in 1:M){             # Loop over all M individuals
    w[i] ~ dbern(omega)      # w = 1 if animal is real/present
    # State model: point process model
    AC[i, 1] ~ dunif(xmin, xmax) # x-coord of activity centre
    AC[i, 2] ~ dunif(ymin, ymax) # y coord of activity centre
    # Observation model: p ~ distance between trap and estimated AC
    for(j in 1:nDetlocs){           # Loop over all detectors
      d2[i,j] <- (AC[i,1] - Detlocs[j,1])^2 +
        (AC[i,2] - Detlocs[j,2])^2           # distance^2
     p[i,j] <- p0 * exp(- d2[i,j]/(2*sigma^2)) # Detection prob
     y[i,j] ~ dbin(p[i,j] * w[i], nOcc)     # The observed data
    }
  }

  # Derived quantities
  N <- sum(w)                       # Population size in state-space (=area)
  D <- N / A                        # Density over state-space
  }
")

# Run the model
# -------------
# Parameters to save
parameters <- c("p0", "sigma", "omega", "N", "D", "AC", "w")

# MCMC settings #try with small number first!!!
ni <- 1000; nb <- 50; nc <- 3; nt <- 5; na <- 1000

#load("FC2021_jagsUI.RData")
out <- jags(jagsData, NULL, parameters, "SCRfc_binorm.jags", n.iter=ni, n.burnin=nb, n.chains=nc,
              n.thin=nt, n.adapt=na, DIC=FALSE, parallel=TRUE)

#save(out, file="FC2021_jagsUI.RData")
# takes 7 mins.

# check convergence and summarize posteriors 
traceplot(out) # Not shown
print(out, 3)

# Check the output
# ----------------
library(mcmcOutput)
( mco <- mcmcOutput(out) )

diagPlot(mco, params=1:5)
plot(mco, 1:5)
View(summary(mco))

# check that N << M
max(mco$N)
M

#plot population with HDI
N <- mco$N
postPlot(N, xlab="Population", showCurve=TRUE)
