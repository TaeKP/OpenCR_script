
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

#A <- popMat$A

## Run function with 1000 simulation stochastic samples
popMat2 <- build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                                              survival_rate = S, survival_rate_sd = SD_survival_rate,
                                              recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                                              init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                                              stochastic = nsim)

A1 <- popMat2$A

# ASYMPTOTIC ANALYSIS #
#---------------------#

## Asymptotic population growth rate 
#  = dominant right eigenvalue of matrix
#lambda <- popbio::lambda(A)

## Matrix element sensitivities
# = derivative of lambda with respect to each matrix element
# = potential effect a small absolute change in matrix element could have on population growth rate
#sens.ME <- popbio::sensitivity(A)

## Matrix element elasticities (= relative sensitivities)
# = sensitivity multiplied by (matrix element/lambda)
# = potential effect a small relative change in matrix element could have on population growth rate
#elas.ME <- popbio::elasticity(A)

#  = dominant right eigenvalue of matrix for matrix A1
lambda1 <- popbio::lambda(A1)

## Matrix element sensitivities for matrix A1
sens.ME1 <- popbio::sensitivity(A1)

## Matrix element elasticities for matrix A1
elas.ME1 <- popbio::elasticity(A1)


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

sens_S1yr <- sum(derivA_S1yr*sens.ME1)

# Derivative matrix for E1yr
derivA_E1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_E1yr[1, 1] <- 0 
derivA_E1yr[2, 1] <- -1 
derivA_E1yr[2, 2] <- -1
derivA_E1yr[1, 2] <- 0

sens_E1yr <- sum(derivA_E1yr*sens.ME1)

# Derivative matrix for f1yr_ad
derivA_f1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_f1yr[1, 1] <- 0 
derivA_f1yr[2, 1] <- 0 
derivA_f1yr[2, 2] <- 0 
derivA_f1yr[1, 2] <- 1

sens_f1yr <- sum(derivA_f1yr*sens.ME1)

## Vital rate elasticities
# = sensitivity multiplied by (vital rate/lambda)
# = potential effect a small relative change in vital could have on population growth rate

elas_S1yr <- sens_S1yr * (popMat2$S1yr/lambda1)
elas_E1yr <- sens_E1yr * (popMat2$E1yr/lambda1)
elas_f1yr <- sens_f1yr * (popMat2$f1yr_ad/lambda1)

## ----------------------------------------------------------------------------------------
## Set up the scenarios with the matrix perturbation 

# Scenarios 1. Reduce 5% juvenile survival
# Derivative matrix for S1yr_juv
derivA_S1yr_juv <- matrix(NA, nrow = 2, ncol = 2)
derivA_S1yr_juv[1, 1] <- 0 
derivA_S1yr_juv[2, 1] <- 1 * 0.95
derivA_S1yr_juv[2, 2] <- 1 
derivA_S1yr_juv[1, 2] <- 0

sens_S1yr_juv <- sum(derivA_S1yr_juv*sens.ME1)

# Scenarios 2. Reduce 5% adult survival
# Derivative matrix for S1yr_ad
derivA_S1yr_ad <- matrix(NA, nrow = 2, ncol = 2)
derivA_S1yr_ad[1, 1] <- 0 
derivA_S1yr_ad[2, 1] <- 1 
derivA_S1yr_ad[2, 2] <- 1 * 0.95
derivA_S1yr_ad[1, 2] <- 0

sens_S1yr_ad <- sum(derivA_S1yr_ad*sens.ME1)

# Scenarios 3. Reduce 5% adult recruitment
# Derivative matrix for f1yr_ad
derivA_f1yr_ad <- matrix(NA, nrow = 2, ncol = 2)
derivA_f1yr_ad[1, 1] <- 0 
derivA_f1yr_ad[2, 1] <- 0 
derivA_f1yr_ad[2, 2] <- 0 
derivA_f1yr_ad[1, 2] <- 1 * 0.95

sens_f1yr_ad <- sum(derivA_f1yr_ad*sens.ME1)

## Vital rate elasticities for 3 scenarios

elas_S1yr_juv <- sens_S1yr_juv * (popMat2$S1yr/lambda1)
elas_S1yr_ad <- sens_S1yr_ad * (popMat2$S1yr/lambda1)
elas_f1yr_ad <- sens_f1yr_ad * (popMat2$f1yr_ad/lambda1)

## Compare Sensitivity and Elasticity Results

Sensitivity_set <- data.frame(
  Scenario = c("Baseline_S", "Baseline_E", "Baseline_f", "Juvenile Survival", "Adult Survival", "Recruitment"),
  Sensitivity = c(sens_S1yr, sens_E1yr, sens_f1yr, sens_S1yr_juv, sens_S1yr_ad, sens_f1yr_ad)
                              )

Elasticity_set <- data.frame(
  Scenario = c("Baseline_S", "Baseline_E", "Baseline_f", "Juvenile Survival", "Adult Survival", "Recruitment"),
  Elasticity = c(elas_S1yr, elas_E1yr, elas_f1yr, elas_S1yr_juv, elas_S1yr_ad, elas_f1yr_ad)
                              )
