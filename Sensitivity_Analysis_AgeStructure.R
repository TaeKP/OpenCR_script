
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

# Checking the structure
str(popMat)
str(A)

## Run function with 1000 simulation stochastic samples
result_popMat <- replicate(nsim, build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                                              survival_rate = S, survival_rate_sd = SD_survival_rate,
                                              recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                                              init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                                              stochastic = stochastic))


# Checking structure
str(result_popMat)

# indexing of each matrix; first 10 slices
result_popMat[c(1, 7, 13, 19, 25, 31, 37, 43, 49, 55)]

# construct the for loop to store the nsim values 

# Selecting matrix every 6 values
selected_indices <- seq(1, length(result_popMat), by = 6)

# Checking the size of matrix 
nrow_mat <- nrow(result_popMat[[selected_indices[1]]])
ncol_mat <- ncol(result_popMat[[selected_indices[1]]])

# Creating empty array
nsim_array <- array(NA, dim = c(nrow_mat, ncol_mat, length(selected_indices)))

# Filling every matrix in each slice
for (k in seq_along(selected_indices)) {
  nsim_array[,,k] <- result_popMat[[selected_indices[k]]]
}

# Checking the result
print(dim(nsim_array))    # array size 2*2*10
print(nsim_array[,,1:10])  # slice 1 to 10
print(nsim_array[,,100])   # slice 100
print(nsim_array[,,1000])   # slice 1000


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

#---------------------#
# Creating the for loop to calculate lambda with 1000 simulations
# Creating vector to store lambda
lambda_values <- numeric(dim(nsim_array)[3])

# Using for loop to calculate lambda from each matrix
for (k in 1:dim(nsim_array)[3]) {
  mat <- nsim_array[,,k]
  lambda_values[k] <- lambda(mat)   
}

# Checking lambda result
print(lambda_values)

# mean
mean_lambda <- mean(lambda_values)

# standard deviation
sd_lambda <- sd(lambda_values)

# SE (standard error) = sd / sqrt(n)
se_lambda <- sd_lambda / sqrt(length(lambda_values))

# 95% CI
ci_lower <- mean_lambda - 1.96 * se_lambda
ci_upper <- mean_lambda + 1.96 * se_lambda

# creating density plot
plot(density(lambda_values),
     main = "Density Plot of Lambda with Mean & 95% CI",
     xlab = "Lambda",
     col = "blue",
     lwd = 2)

# adding mean
abline(v = mean_lambda, col = "red", lwd = 2, lty = 2)

# adding 95% CI
abline(v = ci_lower, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper, col = "darkgreen", lwd = 2, lty = 3)

legend("topright",
       legend = c("Density", "Mean", "95% CI"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3),
       lwd = 2)

#---------------------#
# Creating the for loop to run matrix element sensitivities with 1000 simulations
# Creating list to store sensitivity matrices
sens_results <- vector("list", dim(nsim_array)[3])

for (k in 1:dim(nsim_array)[3]) {
  mat <- nsim_array[,,k]
  sens_results[[k]] <- sensitivity(mat)
}

# Checking slice 1
print(sens_results[[1]])

# sensitivity of element (1,1) in every slice
s11 <- sapply(sens_results, function(x) x[1,1])
print(s11)

(s12 <- sapply(sens_results, function(x) x[1,2]))
(s21 <- sapply(sens_results, function(x) x[2,1]))
(s22 <- sapply(sens_results, function(x) x[2,2]))

#---------------------#
# Creating the for loop to run matrix element elasticities with 1000 simulations
# Creating list to store elasticity
elas_results <- vector("list", dim(nsim_array)[3])

for (k in 1:dim(nsim_array)[3]) {
  mat <- nsim_array[,,k]
  elas_results[[k]] <- elasticity(mat)
}

# # Checking slice 1
print(elas_results[[1]])

# elasticity of element (1,1) in every slice
e11 <- sapply(elas_results, function(x) x[1,1])
print(e11)

(e12 <- sapply(elas_results, function(x) x[1,2]))
(e21 <- sapply(elas_results, function(x) x[2,1]))
(e22 <- sapply(elas_results, function(x) x[2,2]))

## ----------------------------------------------------------------------------------------

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
