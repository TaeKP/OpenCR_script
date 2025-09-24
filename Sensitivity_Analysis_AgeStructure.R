
library(popbio)

# SET VITAL RATE & POPULATION PARAMETERS FROM SCR #
#-------------------------------------------------#

## Inital population proportion
init_adultProp_mean <- 0.89 # Mean of Proportion of initial population that is adult; 3 years
init_adultProp_SD <- 0.05 # SD of Proportion of initial population that is adult

## Growth rate (lambda)
lambda <- 1.16
SE_growth_rate <- 0.12
SD_growth_rate <- SE_growth_rate

## Recruitment (f)
f <- 0.76
SE_recruitment <- 0.12
SD_recruitment <- SE_recruitment

## Survival rate (S)
S <- 0.49
SE_survival_rate <- 0.091
SD_survival_rate <- SE_survival_rate


# BUILDING PROJECTION MATRICES #
#-------------------------------#

## Source function for building projection matrix
source("build_popMat.R")

## Set simulation parameters
stochastic <- TRUE
nsim <- 1000

## Run function (for nsim stochastic samples)
popMat <- build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                       survival_rate = S, survival_rate_sd = SD_survival_rate,
                       recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                       init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                       stochastic = stochastic)

A <- popMat$A


# ASYMPTOTIC ANALYSIS #
#---------------------#

## Asymptotic population growth rate 
#  = dominant right eigenvalue of matrix
lambda <- popbio::lambda(A)

## Matrix element sensitivities
# = derivative of lambda with respect to each matrix element
# = potential effect a small absolute change in matrix element could have on population growth rate
sens.ME <- popbio::sensitivity(A)

## Matrix element elasticities (= relative sensitivities)
# = sensitivity multiplied by (matrix element/lambda)
# = potential effect a small relative change in matrix element could have on population growth rate
elas.ME <- popbio::elasticity(A)


## Vital rate sensitivities
# = derivative of lambda with respect to each vital rate
# = potential effect a small absolute change in vital rate could have on population growth rate

# Projection matrix
#A[1, 1] <- 0 
#A[2, 1] <- S1yr - E1yr 
#A[2, 2] <- S1yr - E1yr 
#A[1, 2] <- f1yr_ad


# Derivative matrix for S1yr
derivA_S1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_S1yr[1, 1] <- 0 
derivA_S1yr[2, 1] <- 1 
derivA_S1yr[2, 2] <- 1 
derivA_S1yr[1, 2] <- 0

sens_S1yr <- sum(derivA_S1yr*sens.ME)

# Derivative matrix for E1yr
derivA_E1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_E1yr[1, 1] <- 0 
derivA_E1yr[2, 1] <- -1 
derivA_E1yr[2, 2] <- -1
derivA_E1yr[1, 2] <- 0

sens_E1yr <- sum(derivA_E1yr*sens.ME)

# Derivative matrix for f1yr_ad
derivA_f1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_f1yr[1, 1] <- 0 
derivA_f1yr[2, 1] <- 0 
derivA_f1yr[2, 2] <- 0 
derivA_f1yr[1, 2] <- 1

sens_f1yr <- sum(derivA_f1yr*sens.ME)

## Vital rate elasticities
# = sensitivity multiplied by (vital rate/lambda)
# = potential effect a small relative change in vital could have on population growth rate

elas_S1yr <- sens_S1yr * (popMat$S1yr/lambda)
elas_E1yr <- sens_E1yr * (popMat$E1yr/lambda)
elas_f1yr <- sens_f1yr * (popMat$f1yr_ad/lambda)


