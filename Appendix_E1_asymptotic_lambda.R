# Appendix E1 - Asymptotic sensitivity analysis code: Matrix element sensitivities and elasticities 
# Last check: 1 June 2026

rm(list=ls(all=TRUE))		# remove all objects

# package require
library(popbio)         
library(truncnorm)    

# Set seed
mySeed <- 90
set.seed(mySeed)

# The parameters were estimated using survey data collected every two years (2019, 2021 & 2023), 
# based on previous research (see Phosri et al., 2025)
# Parameters and uncertainty from open population Spatial Capture-Recapture

# Initial population proportion
init_adultProp_mean <- 0.89 # Mean of Proportion of initial population that is adult; 3 years
init_adultProp_SD <- 0.05 # SD of Proportion of initial population that is adult

# Growth rate (lambda)
lambda <- 1.16
SE_growth_rate <- 0.12
SD_growth_rate <- SE_growth_rate

# Recruitment (f)
f <- 0.76
SE_recruitment <- 0.12
SD_recruitment <- SE_recruitment

# Survival rate (S)
S <- 0.49
SE_survival_rate <- 0.091
SD_survival_rate <- SE_survival_rate
#-----------------------------------------------------------------------------------------------------
# Building projection matrices

build_popMat <- function(growth_rate, growth_rate_sd,
                         survival_rate, survival_rate_sd,
                         recruitment_rate, recruitment_rate_sd,
                         init_adultProp, init_adultProp_SD,
                         stochastic = TRUE) {
  
  # Sample/set demographic parameters
  
  if(stochastic){
  
    # True survival (2 years interval)
    S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) 
    
    # Recruitment (2 years interval)
    f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) 
    
    # Population growth rate
    lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) 
    
    # Initial adult Proportion
    init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
    
  }else{
        S <- survival_rate
        f <- recruitment_rate
        lambda <- growth_rate
        init_adultProp <- init_adultProp_mean
      }
  
        E = abs(lambda - S - f) # Emigration rate
        S1yr <- S / sqrt(lambda)
        f1yr <- f / sqrt(lambda)
        E1yr <- E / sqrt(lambda)
  
  # Build projections matrix
  # (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
  
        A <- matrix(NA, nrow = 2, ncol = 2)
        A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
        A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
        A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
        A[1, 2] <- f1yr # Adults producing juveniles within 1 year
  
  # Calculate "per adult recruitment rate" for time-step
        f1yr_ad <- f1yr /  init_adultProp
        A[1, 2] <- f1yr_ad
  
  
  return(list(A = A,
              S1yr = S1yr,
              f1yr = f1yr,
              f1yr_ad = f1yr_ad,
              E1yr = E1yr,
              init_adultProp = init_adultProp))
}

# Set simulation parameters
nsim <- 10000

# Build population projection matrix for mean values
popMat <- build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                       survival_rate = S, survival_rate_sd = SD_survival_rate,
                       recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                       init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                       stochastic = FALSE)

# Checking the structure
str(popMat)
A <- popMat$A
str(A)

# Assemble nsim replicate matrices
nsim_VRs <- list()
nsim_array <- array(dim = c(nrow(A), ncol(A), nsim))

for(i in 1:nsim){
              result_popMat <- build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                                survival_rate = S, survival_rate_sd = SD_survival_rate,
                                recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                                init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                                stochastic = TRUE)
  
              nsim_array[,,i] <- result_popMat$A
              result_popMat$A <- NULL
              nsim_VRs <- append(nsim_VRs, list(result_popMat))
}

# Checking the result
print(dim(nsim_array))    # array size 2*2*10000

#-----------------------------------------------------------------------------------------------------
# Asymptotic analysis
#-----------------------------------------------------------------------------------------------------
# Asymptotic population growth rate 
# Creating the for loop to calculate lambda with 10000 simulations
# Creating vector to store lambda
lambda_values <- numeric(nsim)

# Using for loop to calculate lambda from each matrix
for(k in 1:dim(nsim_array)[3]) {
        mat <- nsim_array[,,k]
        lambda_values[k] <- popbio::lambda(mat)   
}

# Checking lambda result
View(lambda_values)

# mean value
(mean_lambda <- mean(lambda_values))

# median value
(median_lambda <- median(lambda_values))

# standard deviation
(sd_lambda <- sd(lambda_values))

# 95% CI
(ci_lower <- quantile(lambda_values, prob = 0.025)) 
(ci_upper <- quantile(lambda_values, prob = 0.975)) 

# Creating density plot for asymptotic population growth rate
plot(density(lambda_values),
     main = "",
     xlab = "Lambda",
     col = "blue",
     lwd = 2)

    # adding median
    abline(v = median_lambda, col = "red", lwd = 2, lty = 3)

    # adding 95% CI
    abline(v = ci_lower, col = "gray", lwd = 2, lty = 2)
    abline(v = ci_upper, col = "gray", lwd = 2, lty = 2)

    legend("topright",
      legend = c("Density", "Median", "95% CI"),
      col = c("blue", "red", "gray"),
      lty = c(1, 3, 2), lwd = 2)

#-----------------------------------------------------------------------------------------------------
# Matrix element sensitivities and elasticities
#-----------------------------------------------------------------------------------------------------
# Creating the for loop to run matrix element sensitivities with 10000 simulations
# Creating list to store sensitivity matrices
sens_results <- vector("list", nsim)
    
for (k in 1:nsim) {
            mat <- nsim_array[,,k]
            sens_results[[k]] <- sensitivity(mat)
}

# Checking slice 1
print(sens_results[[1]])
    
# sensitivity of element (1,1) in every slice
s11 <- sapply(sens_results, function(x) x[1,1])
View(s11)
hist(s11, main = "")
med <- median(s11)
abline(v = med, col = "red", lwd = 2, lty = 3)

# Sensitivities of element (1,2), (2,1) and (2,2) in every slice    
s12 <- sapply(sens_results, function(x) x[1,2])
View(s12)
s21 <- sapply(sens_results, function(x) x[2,1])
View(s21)
s22 <- sapply(sens_results, function(x) x[2,2])
View(s22)

# Matrix element elasticities 
# Creating the for loop to run matrix element elasticities with 10000 simulations
# Creating list to store elasticity
elas_results <- vector("list", nsim)

for (k in 1:nsim) {
            mat <- nsim_array[,,k]
            elas_results[[k]] <- elasticity(mat)
}

# # Checking slice 1
print(elas_results[[1]])

# elasticity of element (1,1) in every slice
e11 <- sapply(elas_results, function(x) x[1,1])
View(e11)

# elasticities of element (1,2), (2,1) and (2,2) in every slice
e12 <- sapply(elas_results, function(x) x[1,2])
View(e12)
e21 <- sapply(elas_results, function(x) x[2,1])
View(e21)
e22 <- sapply(elas_results, function(x) x[2,2])
View(e22)

# need object "appendix_E1_R" for Appendix E2
save.image(file="appendix_E1.RData")
#----------------------------------------- END NOT RUN------------------------------------------------------------
