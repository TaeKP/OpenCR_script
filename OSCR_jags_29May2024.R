# FC-KSRY data
# 4 dimensional array, i:j:k:t
# prepare data for open-SCR by JAGS
# 3 years period, 2019, 2021, 2023
# totally 127 fishing cats
# 104 trap locations
# aggregated capture history as 5-monthly occasions

# Developed from black bear code Chandler,R.B.,& Clark,J.D.(2014). Spatially explicit integrated population models. Methods in Ecology and Evolution, 5(12), 1351-1360.

# getwd()
#=================================================
library(jagsUI)
library(makeJAGSmask)
library(mcmcOutput)
library(wiqid)
library(spatstat) # Functions to manipulate the polygon
library(MASS)
library(secr)

#remotes::install_github("mikemeredith/makeJAGSmask")
library(raster)

# load input data
load("FC_array_3yrs.RData")

print(FC_arr_3years) 
str(FC_arr_3years) 
 
# check the data by indexing the specific occasion and year

FC_arr_3years[,,"Occ1", "2019"]

FC_arr_3years[,,"Occ1", "2021"]

FC_arr_3years[,,"Occ1", "2023"]

# 127 individuals/104 traps/over 5 occasions/3 years
table(FC_arr_3years)  #total detection
#     0     1 
# 82115   435

# ------------------------------------
# Trap locations:
KSRY.trap <- read.traps("trap3yrs.txt", detector = "proximity")     # meter scale
head(KSRY.trap)

# Convert to a matrix
Detlocs <- as.matrix(KSRY.trap)
# Detlocs <- as.matrix(KSRY.trap)/1000        # convert to km scale 

# just check the trap station is correct!
( nDetlocs <- nrow(Detlocs) )
# 104
str(Detlocs)

# -------------------------------------
# Working with state-space area
# we can generate our KSRY patch here!! follows the closed model

library(scales)
plot(Detlocs, xlab='x coordinate', ylab='y coordinate', frame=FALSE, las=1, pch=16,
     col='#002A64FF', asp=1)

# Define limits of the state-space and add it to the plot
xlims <- c(580000, 610000)
ylims <- c(1380000, 1310000)

# km scale
#xlims <- c(585, 610)
#ylims <- c(1370, 1320)

rect(xlims[1], ylims[1], xlims[2], ylims[2], col=alpha('grey', 0.3), border=NA)

# -------------------------------------
# Create the 'secr' mask: obtain the polygon shape from closed model
require(maptools)

SRY.shape <- readShapePoly('KSRY_mask_hab_22Jul19')
plot(SRY.shape)
points(KSRY.trap, pch=3, col='red')

# Making the mask followed the MLE state space # buffer 6km

### load the result here!!
load("SRY_mask.RData")

#SRY.mask<- make.mask(KSRY.trap, buffer = 6000, spacing = 500, 
#                     type = 'trapbuffer',  poly= SRY.shape)

#save(SRY.mask, file="SRY_mask.RData")  # buffer 6 km

head(SRY.mask)
plot(SRY.mask)
points(KSRY.trap, pch=3, col='red')

# ------------------------------------
# Generate 'habMat':
# Mike packages 'makeJAGSmask'

mymask <- convertMask(SRY.mask, KSRY.trap, plot = TRUE) 
class(mymask)
head(mymask$habMat)
str(mymask)

### load the result here!!
load("SRY_mymask.RData")
# -------------------------------------
# Replace the missing values by zero  # check data

fc3y <- FC_arr_3years

fc.sum <- fc3y
fc.sum[is.na(fc.sum)] <- 0

# Produce table
fcXtrap <- apply(fc.sum, c(1, 2), sum)

## total detection per FC and traps
table(fcXtrap)
#     0     1     2     3     4     5 
# 13001    85    53    36    29     4

## total number of times each FC was detected
table(apply(fcXtrap, 1, sum))
#  1  2  3  4  5  6  7  8  9 10 12 13 
# 32 22 23 17 12  8  5  2  2  2  1  1 
## 32FCs were detected on only single occasion, only 1FC was detected on a total of 13 occasions 

#========================================================================================
# Loop over 127 individuals and map their locations (individually)
# combined over all 3 years (which represent a 5-year duration)
# op <- par(ask=TRUE)                               # Press Esc to interrupt

op <- par(ask=dev.interactive(orNone=TRUE))         # Only ask if plotting on-screen
for (i in 1:nrow(fcXtrap)){
  cat(paste("\n\n*** Plot fishing cat number", i, "***\n\n"))
  lucky.traps <- which(fcXtrap[i,] > 0)
  plot(Detlocs, xlab='x coordinate', ylab='y coordinate', frame=FALSE, pch=1,
       col=rgb(0,0,0,1), cex=1.5, asp=1, main=paste('Fishing cat number', i))
  points((Detlocs)[lucky.traps,1], (Detlocs)[lucky.traps,2], pch=16, col='red', cex=1.5)
}

par(op)

#========================================================================================

# The integrated population model
# ====================================
# Do data augmentation of the SCR data

M <- 600
y <- array(0, dim=c(M, dim(fc3y)[2:4]))
y[1:dim(fc3y)[1],,,] <- fc3y

# Bundle data and produce data overview
jags.data <- with (mymask, list(y=y, T=3, K=dim(y)[3], J=dim(y)[2], M=M, x=Detlocs,
                  habMat=habMat, trapMat=trapMat, upperLimit=upperLimit,
                  ones=rep(1, M)))

str(jags.data)

# ====================================
# Added JAGS model without the occupancy data here!
# add the habmat here!

# # adjust p0 and sigma with informative prior ver.5

cat(file="model_fc_habmat_adjust5", "
model {
  # Priors
  psi ~ dunif(0, 1)                               # M*psi = E[N(1)], for data augmentation
  phi ~ dunif(0, 1)                               # Survival
  gamma ~ dunif(0, 3)                             # Per-capita recruitment
  
  # Prior for baseline cap prob # binomial 0 to 1
  p0 ~ dbeta(1.3, 3.7)
  
  # Prior for scale parameter # poisson dist
  sigma ~ dgamma(3, 1)
  
  # Derived parameters
  for (t in 1:T){
    N[t] <- sum(z[,t])                            # Population size at time t
    EB[t] <- N[t] * gamma                         # Expected number of recruits
    A[t] <- max(M - sum(a[,t]), 0.001)            # Bears available to be recruited
    b[t] <- min(EB[t] / A[t], 0.999)              # Probability of being recruited
  }
  
  # Population model (individual-based model)
  # Initial state
  for (i in 1:M){
    z[i,1] ~ dbern(psi)                           # Is a member of M alive?
    a[i,1] <- z[i,1]                              # Recruited yet?
    
    # State model: point process model
    AC[i, 1] ~ dunif(1, upperLimit[1])    # x-coord of activity centers
    AC[i, 2] ~ dunif(1, upperLimit[2])    # y coord of activity centers
    hab[i] <- habMat[trunc(AC[i, 1]), trunc(AC[i, 2])] # habitat look-up
    ones[i] ~ dbern(hab[i])               # ones[i] = 1, the ones trick
    
    #s[i,1] ~ dunif(xlims[1], xlims[2])            # Activity center, homogeneous
    #s[i,2] ~ dunif(ylims[1], ylims[2])
    
    # Dynamics over time
    for (t in 2:T){
      z[i,t] ~ dbern(z[i,t-1] * phi + (1 - a[i,t-1]) * b[t-1])
      a[i,t] <- max(z[i,1:t])
    } #t
  } #i
  
  # Observation models
  # Spatial capture-recapture data
  for (i in 1:M){
    for (j in 1:J){
      #d2[i,j] <- (s[i,1] - x[j,1])^2 + (s[i,2] - x[j,2])^2  # Distance squared
      d2[i,j] <- (AC[i,1] - trapMat[j,1])^2 + (AC[i,2] - trapMat[j,2])^2           # distance^2
      p[i,j] <- p0 * exp(-d2[i,j] / (2 * sigma^2))
    } #j
    for (j in 1:J){                               # Trap
      for (k in 1:K){                             # Secondary occasions
        # The years with SCR data
        y[i,j,k,1] ~ dbern(p[i,j] * z[i,1])
        y[i,j,k,2] ~ dbern(p[i,j] * z[i,2])
        y[i,j,k,3] ~ dbern(p[i,j] * z[i,3])
      } #k
    } #j
    zi[i] <- (sum(z[i,]) > 0)                     # Was this bear ever alive?
  } #i
  
  # Occupancy data
  #for (j in 1:J){
  #  for (k in 1:K){
  #    o[j,k,1] ~ dbern(1-prod(1-p[,j] * z[,2]))
  #    o[j,k,2] ~ dbern(1-prod(1-p[,j] * z[,4]))
  #    o[j,k,3] ~ dbern(1-prod(1-p[,j] * z[,6]))
  #  } #k
  #} #j
  N.ever <- sum(zi[])                             # Bears ever alive – superpopulation size
}
")
# ====================================
# Initial values
# not to provide z=1 for all individuals, because there would be no individual that is ready to be recruited
# and MCMC will not work properly
# initial values adjusted from black bear code

zin <- array(0, dim=c(dim(y)[1],3))
zin[1:200,] <- 1

inits2 <- function() {list(z=zin, phi=runif(1, 0.1, 0.9), p0=runif(1, 0.01, 0.1),
                           psi=runif(1, 0.2, 0.4), sigma=runif(1, 1, 4), gamma=runif(1, 0.1, 1))}

# Parameters monitored
# parameters <- c("psi", "phi", "gamma", "p0", "sigma", "N", "EB", "b", "A", "s", "zi", "z", "N.ever")
parameters <- c("phi", "gamma", "p0", "sigma", "N", "EB", "b", "A", "N.ever")

# MCMC settings
ni <- 15000; nb <- 5000; nc <- 6; nt <- 1; na <- 2000  # for testing

# ----------------
# Fit with JAGS model
# Adjust p0 and sigma with informative prior (see model_fc_habmat_adjust5)

Sys.time()
out5 <- jags(jags.data, inits2, parameters, "model_fc_habmat_adjust5", 
             n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
Sys.time()

#save(out5, file="28May2024_FC_JAGS_ni15e3_model_fc_habmat_adjust5.RData")

print(out5, 3)
# traceplot(out1)

# ~14 hrs
### load the result here!!
load("28May2024_FC_JAGS_ni15e3_model_fc_habmat_adjust5.RData")

# ----------------
# Check the MCMC output

library(mcmcOutput)
mco <- mcmcOutput(out5) 
summary(mco)                  #see the mcmc error

#              mean     sd   median      l95      u95  Rhat MCEpc
# phi         0.311  0.052    0.310    0.209    0.411 1.004 0.813
# gamma       0.989  0.112    0.984    0.775    1.212 1.001 0.808
# p0          0.379  0.044    0.371    0.306    0.470 1.506 5.177
# sigma       3.046  0.210    3.082    2.639    3.389 1.817 5.697

diagPlot(mco, params=1:7)
plot(mco, 1:7)

# check that N << M
max(mco$N.ever)
M

# plot super population with HDI
superN <- mco$N.ever
postPlot(superN, xlab="Super Population", showCurve=TRUE)

# ----------------
# plot the result of year-specific population

plot(y=out5$mean$N, x=1:3, ylim=c(0,150), type="b", pch=16, axes=FALSE,
     ylab="Population size (N)", xlab="Year")

segments(1:3, out5$q2.5$N, 1:3, out5$q97.5$N)
axis(1, at=1:3, labels = c(2019, 2021, 2023))
axis(2, las=1)
