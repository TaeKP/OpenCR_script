
library(popbio)

## Set seed
mySeed <- 90
set.seed(mySeed)

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

## Build population projection matrix for mean values
popMat <- build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                       survival_rate = S, survival_rate_sd = SD_survival_rate,
                       recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                       init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                       stochastic = FALSE)

A <- popMat$A

# Checking the structure
str(popMat)
str(A)

## Assemble nsim replicate matrices
result_popMat <- replicate(nsim, build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                                           survival_rate = S, survival_rate_sd = SD_survival_rate,
                                           recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                                           init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                                           stochastic = stochastic))

nsim_VRs <- list()
nsim_array <- array(dim = c(nrow(A), ncol(A), nsim))


for(i in 1:nsim){
  
  result_popMat <- build_popMat(growth_rate = lambda, growth_rate_sd = SD_growth_rate,
                                survival_rate = S, survival_rate_sd = SD_survival_rate,
                                recruitment_rate = f, recruitment_rate_sd = SD_recruitment,
                                init_adultProp = init_adultProp_mean, init_adultProp_SD = init_adultProp_SD,
                                stochastic = stochastic)
  
  nsim_array[,,i] <- result_popMat$A
  
  result_popMat$A <- NULL
  nsim_VRs <- append(nsim_VRs, list(result_popMat))
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
lambda_values <- numeric(nsim)

# Using for loop to calculate lambda from each matrix
for (k in 1:dim(nsim_array)[3]) {
  mat <- nsim_array[,,k]
  lambda_values[k] <- popbio::lambda(mat)   
}

# Checking lambda result
print(lambda_values)

# mean
mean_lambda <- mean(lambda_values)

# median
median_lambda <- median(lambda_values)

# standard deviation
sd_lambda <- sd(lambda_values)

# 95% CI
ci_lower <- quantile(lambda_values, prob = 0.025) 
ci_upper <- quantile(lambda_values, prob = 0.975) 

# creating density plot
plot(density(lambda_values),
     main = "Density Plot of Lambda with Mean & 95% CI",
     xlab = "Lambda",
     col = "blue",
     lwd = 2)

# adding mean
abline(v = mean_lambda, col = "red", lwd = 2, lty = 2)

# adding median
abline(v = median_lambda, col = "red", lwd = 2, lty = 3)

# adding 95% CI
abline(v = ci_lower, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper, col = "darkgreen", lwd = 2, lty = 3)

legend("topright",
       legend = c("Density", "Mean", "Median", "95% CI"),
       col = c("blue", "red", "red", "darkgreen"),
       lty = c(1, 2, 3, 3),
       lwd = 2)

#---------------------#
# Creating the for loop to run matrix element sensitivities with 1000 simulations
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
print(s11)

(s12 <- sapply(sens_results, function(x) x[1,2]))
(s21 <- sapply(sens_results, function(x) x[2,1]))
(s22 <- sapply(sens_results, function(x) x[2,2]))

#---------------------#
# Creating the for loop to run matrix element elasticities with 1000 simulations
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

elas_S1yr <- sens_S1yr * (popMat$S1yr/lambda)
elas_E1yr <- sens_E1yr * (popMat$E1yr/lambda)
elas_f1yr <- sens_f1yr * (popMat$f1yr_ad/lambda)

## ----------------------------------------------------------------------------------------
## Vital rate sensitivities for 1000 samples

## 1.1 Vital rate sensitivity for survival
# Using "lapply" function for list object, then obtained the sens_results (1000 samples) 
# multiply by the derivative of matrix A for survival (S1yr) (1 matrix) and using sum function to calculate each value in each slice

# Multiply each matrix in the list by derivA_S1yr
# get 1000 values of sens_S1yr_sims from sum function
sens_S1yr_sims <- lapply(sens_results, function(x) sum(derivA_S1yr %*% x))

# This is sensitivity of survival from 1000 samples
sens_S1yr_sims

## 1.2 Estimating the elasticity for the Survival
# Building the object for the S1yr values
result_S1 <- c()
for(i in 1:nsim){
  result_S1 <- c(result_S1, nsim_VRs[[i]]$S1yr)
}

# Checking the structure of result_S1 for list object include 1000 values
str(result_S1)

# Initialize empty list to store results
elas_S1yr_sims <- vector("list", length = 1000)

# Loop through positions
for (i in seq_along(sens_S1yr_sims)) {
  elas_S1yr_sims[[i]] <- sens_S1yr_sims[[i]] * (result_S1[i] / lambda_values[i])
}

# This is elasticity of survival from 1000 samples
elas_S1yr_sims

# Plot the density of sensitivity of survival
# Flatten into a numeric vector
vec_sens_S1yr_sims <- unlist(sens_S1yr_sims)

library(ggplot2)

# make dataframe
df_vec_sens_S1yr_sims <- data.frame(value = vec_sens_S1yr_sims)

ggplot(df_vec_sens_S1yr_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_sens_S1yr_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of sensitivity of survival")

# mean
mean_sens_S1yr_sims <- mean(df_vec_sens_S1yr_sims$value)
mean_sens_S1yr_sims

# standard deviation
sd_sens_S1yr_sims <- sd(df_vec_sens_S1yr_sims$value)
sd_sens_S1yr_sims

# SE (standard error) = sd / sqrt(n)
se_sens_S1yr_sims <- sd_sens_S1yr_sims / sqrt(length(df_vec_sens_S1yr_sims$value))
se_sens_S1yr_sims

# 95% CI
ci_lower_sens_S1yr_sims <- mean_sens_S1yr_sims - 1.96 * se_sens_S1yr_sims
ci_lower_sens_S1yr_sims
ci_upper_sens_S1yr_sims <- mean_sens_S1yr_sims + 1.96 * se_sens_S1yr_sims
ci_upper_sens_S1yr_sims

# Plot the density of elasticity of survival
# Flatten into a numeric vector
vec_elas_S1yr_sims <- unlist(elas_S1yr_sims)

# make dataframe
df_vec_elas_S1yr_sims <- data.frame(value = vec_elas_S1yr_sims)

ggplot(df_vec_elas_S1yr_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_elas_S1yr_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of elasticity of survival")

# mean
mean_elas_S1yr_sims <- mean(df_vec_elas_S1yr_sims$value)
mean_elas_S1yr_sims

# standard deviation
sd_elas_S1yr_sims <- sd(df_vec_elas_S1yr_sims$value)
sd_elas_S1yr_sims

# SE (standard error) = sd / sqrt(n)
se_elas_S1yr_sims <- sd_elas_S1yr_sims / sqrt(length(df_vec_elas_S1yr_sims$value))
se_elas_S1yr_sims

# 95% CI
ci_lower_elas_S1yr_sims <- mean_elas_S1yr_sims - 1.96 * se_elas_S1yr_sims
ci_lower_elas_S1yr_sims
ci_upper_elas_S1yr_sims <- mean_elas_S1yr_sims + 1.96 * se_elas_S1yr_sims
ci_upper_elas_S1yr_sims

## -----------------------------------------
## 2.1 Vital rate sensitivities for emigration
# Using "lapply" function for list object, then obtained the sens_results (1000 samples) 
# multiply by the derivative of matrix A for emigration (E1yr) (1 matrix)

# Multiply each matrix in the list by derivA_E1yr
# get 1000 values of sens_E1yr_sims from sum function
sens_E1yr_sims <- lapply(sens_results, function(x) sum(derivA_E1yr %*% x))

# This is sensitivity of emigration from 1000 samples
sens_E1yr_sims

## 2.2 Estimating the elasticity for the Emigration
# Building the object for the E1yr values
result_E1 <- c()
for(i in 1:nsim){
  result_E1 <- c(result_E1, nsim_VRs[[i]]$E1yr)
}

# Checking the structure of result_E1 for list object include 1000 values
str(result_E1)

# Initialize empty list to store results
elas_E1yr_sims <- vector("list", length = 1000)

# Loop through positions
for (i in seq_along(sens_E1yr_sims)) {
  elas_E1yr_sims[[i]] <- sens_E1yr_sims[[i]] * (result_E1[[i]] / lambda_values[i])
}

# This is elasticity of emigration from 1000 samples
elas_E1yr_sims

# Plot the density of sensitivity of emigration
# Flatten into a numeric vector
vec_sens_E1yr_sims <- unlist(sens_E1yr_sims)

# make dataframe
df_vec_sens_E1yr_sims <- data.frame(value = vec_sens_E1yr_sims)

ggplot(df_vec_sens_E1yr_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_sens_E1yr_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of sensitivity of emigration")

# mean
mean_sens_E1yr_sims <- mean(df_vec_sens_E1yr_sims$value)
mean_sens_E1yr_sims

# standard deviation
sd_sens_E1yr_sims <- sd(df_vec_sens_E1yr_sims$value)
sd_sens_E1yr_sims

# SE (standard error) = sd / sqrt(n)
se_sens_E1yr_sims <- sd_sens_E1yr_sims / sqrt(length(df_vec_sens_E1yr_sims$value))
se_sens_E1yr_sims

# 95% CI
ci_lower_sens_E1yr_sims <- mean_sens_E1yr_sims - 1.96 * se_sens_E1yr_sims
ci_lower_sens_E1yr_sims
ci_upper_sens_E1yr_sims <- mean_sens_E1yr_sims + 1.96 * se_sens_E1yr_sims
ci_upper_sens_E1yr_sims

# Plot the density of elasticity of emigration
# Flatten into a numeric vector
vec_elas_E1yr_sims <- unlist(elas_E1yr_sims)

# make dataframe
df_vec_elas_E1yr_sims <- data.frame(value = vec_elas_E1yr_sims)

ggplot(df_vec_elas_E1yr_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_elas_E1yr_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of elasticity of emigration")

# mean
mean_elas_E1yr_sims <- mean(df_vec_elas_E1yr_sims$value)
mean_elas_E1yr_sims

# standard deviation
sd_elas_E1yr_sims <- sd(df_vec_elas_E1yr_sims$value)
sd_elas_E1yr_sims

# SE (standard error) = sd / sqrt(n)
se_elas_E1yr_sims <- sd_elas_E1yr_sims / sqrt(length(df_vec_elas_E1yr_sims$value))
se_elas_E1yr_sims

# 95% CI
ci_lower_elas_E1yr_sims <- mean_elas_E1yr_sims - 1.96 * se_elas_E1yr_sims
ci_lower_elas_E1yr_sims
ci_upper_elas_E1yr_sims <- mean_elas_E1yr_sims + 1.96 * se_elas_E1yr_sims
ci_upper_elas_E1yr_sims

## -----------------------------------------
## 3.1 Vital rate sensitivities for recruitment
# Using "lapply" function for list object, then obtained the sens_results (1000 samples) 
# multiply by the derivative of matrix A for recruitment (f1yr_ad) (1 matrix)

# Multiply each matrix in the list by derivA_f1yr
sens_f_samples <- lapply(sens_results, function(x) derivA_f1yr %*% x)
sens_f_samples

# get 1000 values which extract only element [1,2] for adult recruitment
sens_f1yr_sims <- sapply(sens_f_samples, function(x) x[1, 2])

length(sens_f1yr_sims)  # must be 1000

# This is sensitivity of recruitment from 1000 samples
sens_f1yr_sims

## 3.2 Estimating the elasticity for the Recruitment
# Building the object for the f1yr_ad values
result_f1 <- c()
for(i in 1:nsim){
  result_f1 <- c(result_S1, nsim_VRs[[i]]$f1yr)
}

# Checking the structure of result_f1 for list object include 1000 values
str(result_f1)

# Initialize empty list to store results
elas_f1yr_sims <- vector("list", length = 1000)

# Loop through positions
for (i in seq_along(sens_f1yr_sims)) {
  elas_f1yr_sims[[i]] <- sens_f1yr_sims[[i]] * (result_f1[[i]] / lambda_values[i])
}

# This is elasticity of recruitment from 1000 samples
elas_f1yr_sims

# Plot the density of sensitivity of recruitment
# Flatten into a numeric vector
vec_sens_f1yr_sims <- unlist(sens_f1yr_sims)

# make dataframe
df_vec_sens_f1yr_sims <- data.frame(value = vec_sens_f1yr_sims)

ggplot(df_vec_sens_f1yr_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_sens_f1yr_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of sensitivity of recruitment")

# mean
mean_sens_f1yr_sims <- mean(df_vec_sens_f1yr_sims$value)
mean_sens_f1yr_sims

# standard deviation
sd_sens_f1yr_sims <- sd(df_vec_sens_f1yr_sims$value)
sd_sens_f1yr_sims

# SE (standard error) = sd / sqrt(n)
se_sens_f1yr_sims <- sd_sens_f1yr_sims / sqrt(length(df_vec_sens_f1yr_sims$value))
se_sens_f1yr_sims

# 95% CI
ci_lower_sens_f1yr_sims <- mean_sens_f1yr_sims - 1.96 * se_sens_f1yr_sims
ci_lower_sens_f1yr_sims
ci_upper_sens_f1yr_sims <- mean_sens_f1yr_sims + 1.96 * se_sens_f1yr_sims
ci_upper_sens_f1yr_sims

# Plot the density of elasticity of recruitment
# Flatten into a numeric vector
vec_elas_f1yr_sims <- unlist(elas_f1yr_sims)

# make dataframe
df_vec_elas_f1yr_sims <- data.frame(value = vec_elas_f1yr_sims)

ggplot(df_vec_elas_f1yr_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_elas_f1yr_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of elasticity of recruitment")

# mean
mean_elas_f1yr_sims <- mean(df_vec_elas_f1yr_sims$value)
mean_elas_f1yr_sims

# standard deviation
sd_elas_f1yr_sims <- sd(df_vec_elas_f1yr_sims$value)
sd_elas_f1yr_sims

# SE (standard error) = sd / sqrt(n)
se_elas_f1yr_sims <- sd_elas_f1yr_sims / sqrt(length(df_vec_elas_f1yr_sims$value))
se_elas_f1yr_sims

# 95% CI
ci_lower_elas_f1yr_sims <- mean_elas_f1yr_sims - 1.96 * se_elas_f1yr_sims
ci_lower_elas_f1yr_sims
ci_upper_elas_f1yr_sims <- mean_elas_f1yr_sims + 1.96 * se_elas_f1yr_sims
ci_upper_elas_f1yr_sims

## compare the 3 sensitivity models
# makes data frame for 3 sensitivity models
df_all_sens_models <- data.frame(
  value = c(vec_sens_S1yr_sims, vec_sens_E1yr_sims, vec_sens_f1yr_sims),
  model = rep(c("Model S", "Model E", "Model f"), each = 1000)
)

df_all_sens_models

# calculate the mean value for each model
library(dplyr)
means_sens <- df_all_sens_models %>%
  group_by(model) %>%
  summarise(mean_val = mean(value), .groups = "drop")

# Plot densities together with each mean line
ggplot(df_all_sens_models, aes(x = value, colour = model, fill = model)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  xlim(-3, 3) +
  scale_colour_grey(start = 0.2, end = 0.8) +   # line colours from dark to light grey
  scale_fill_grey(start = 0.8, end = 0.2) +     # fill shades reversed for contrast
  theme_minimal(base_size = 13) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_vline(data = means_sens, aes(xintercept = mean_val, linetype = model),
             colour = "black", linewidth = 0.1)+
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Density comparison of 3 sensitivity models",
       x = "Value", y = "Density")

# makes data frame for 3 elasticity models
df_all_elas_models <- data.frame(
  value = c(vec_elas_S1yr_sims, vec_elas_E1yr_sims, vec_elas_f1yr_sims),
  model = rep(c("Model S", "Model E", "Model f"), each = 1000)
)

df_all_elas_models

# calculate the mean value for each model
means_elas <- df_all_elas_models %>%
  group_by(model) %>%
  summarise(mean_val = mean(value), .groups = "drop")

# Plot densities together with each mean line
ggplot(df_all_elas_models, aes(x = value, colour = model, fill = model)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  xlim(-3, 3) +
  scale_colour_grey(start = 0.2, end = 0.8) +   # line colours from dark to light grey
  scale_fill_grey(start = 0.8, end = 0.2) +     # fill shades reversed for contrast
  theme_minimal(base_size = 13) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_vline(data = means_elas, aes(xintercept = mean_val, linetype = model),
             colour = "black", linewidth = 0.1)+
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Density comparison of 3 elasticity models",
       x = "Value", y = "Density")

## ----------------------------------------------------------------------------------------
## Set up the scenarios with the matrix perturbation
# Set up the value for perturbation
# which reduce 5% of each vital rates (juvenile survival/adult survival/adult recruitment)
per_fact <- 0.95
per_fact

# 4.1 Scenarios 1. Reduce 5% Juvenile survival
# Derivative matrix for S1yr_juv
derivA_S1yr_juv <- matrix(NA, nrow = 2, ncol = 2)
derivA_S1yr_juv[1, 1] <- 0 
derivA_S1yr_juv[2, 1] <- 1 * per_fact
derivA_S1yr_juv[2, 2] <- 1 
derivA_S1yr_juv[1, 2] <- 0

derivA_S1yr_juv

# Using "lapply" function for list object, then obtained the sens_results (1000 samples) 
# multiply by the derivative of matrix A for Juvenile survival (S1yr_juv) (1 matrix)

sens_perJS_sims <- lapply(sens_results, function(x) sum(derivA_S1yr_juv %*% x))

# This is sensitivity of Juvenile survival from 1000 samples
sens_perJS_sims

# 4.2 Estimating the elasticity for the perturbation Juvenile survival
# Initialize empty list to store results
elas_perJS_sims <- vector("list", length = 1000)

# Loop through positions
for (i in seq_along(sens_perJS_sims)) {
  elas_perJS_sims[[i]] <- sens_perJS_sims[[i]] * (result_S1[[i]] / lambda_values[i])
}

# This is elasticity of perturbation Juvenile survival from 1000 samples
elas_perJS_sims

# Plot the density of sensitivity of perturbation Juvenile survival
# Flatten into a numeric vector
vec_sens_perJS_sims <- unlist(sens_perJS_sims)

# make dataframe
df_vec_sens_perJS_sims <- data.frame(value = vec_sens_perJS_sims)

ggplot(df_vec_sens_perJS_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_sens_perJS_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of sensitivity of perturbation Juvenile survival")

# mean
mean_sens_perJS_sims <- mean(df_vec_sens_perJS_sims$value)
mean_sens_perJS_sims

# standard deviation
sd_sens_perJS_sims <- sd(df_vec_sens_perJS_sims$value)
sd_sens_perJS_sims

# SE (standard error) = sd / sqrt(n)
se_sens_perJS_sims <- sd_sens_perJS_sims / sqrt(length(df_vec_sens_perJS_sims$value))
se_sens_perJS_sims

# 95% CI
ci_lower_sens_perJS_sims <- mean_sens_perJS_sims - 1.96 * se_sens_perJS_sims
ci_lower_sens_perJS_sims
ci_upper_sens_perJS_sims <- mean_sens_perJS_sims + 1.96 * se_sens_perJS_sims
ci_upper_sens_perJS_sims

# Plot the density of elasticity of perturbation Juvenile survival
# Flatten into a numeric vector
vec_elas_perJS_sims <- unlist(elas_perJS_sims)

# make dataframe
df_vec_elas_perJS_sims <- data.frame(value = vec_elas_perJS_sims)

ggplot(df_vec_elas_perJS_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_elas_perJS_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of elasticity of perturbation Juvenile survival")

# mean
mean_elas_perJS_sims <- mean(df_vec_elas_perJS_sims$value)
mean_elas_perJS_sims

# standard deviation
sd_elas_perJS_sims <- sd(df_vec_elas_perJS_sims$value)
sd_elas_perJS_sims

# SE (standard error) = sd / sqrt(n)
se_elas_perJS_sims <- sd_elas_perJS_sims / sqrt(length(df_vec_elas_perJS_sims$value))
se_elas_perJS_sims

# 95% CI
ci_lower_elas_perJS_sims <- mean_elas_perJS_sims - 1.96 * se_elas_perJS_sims
ci_lower_elas_perJS_sims
ci_upper_elas_perJS_sims <- mean_elas_perJS_sims + 1.96 * se_elas_perJS_sims
ci_upper_elas_perJS_sims

## -----------------------------------------
# 5.1 Scenarios 2. Reduce 5% Adult survival
# Derivative matrix for S1yr_ad
derivA_S1yr_ad <- matrix(NA, nrow = 2, ncol = 2)
derivA_S1yr_ad[1, 1] <- 0 
derivA_S1yr_ad[2, 1] <- 1 
derivA_S1yr_ad[2, 2] <- 1 * per_fact
derivA_S1yr_ad[1, 2] <- 0

derivA_S1yr_ad

# Using "lapply" function for list object, then obtained the sens_results (1000 samples) 
# multiply by the derivative of matrix A for Adult survival (S1yr_ad) (1 matrix)

sens_perAS_sims <- lapply(sens_results, function(x) sum(derivA_S1yr_ad %*% x))

# This is sensitivity of adult survival from 1000 samples
sens_perAS_sims

# 5.2 Estimating the elasticity for the perturbation Adult survival
# Initialize empty list to store results
elas_perAS_sims <- vector("list", length = 1000)

# Loop through positions
for (i in seq_along(sens_perAS_sims)) {
  elas_perAS_sims[[i]] <- sens_perAS_sims[[i]] * (result_S1[[i]] / lambda_values[i])
}

# This is elasticity of perturbation Adult survival from 1000 samples
elas_perAS_sims

# Plot the density of sensitivity of perturbation Adult survival
# Flatten into a numeric vector
vec_sens_perAS_sims <- unlist(sens_perAS_sims)

# make dataframe
df_vec_sens_perAS_sims <- data.frame(value = vec_sens_perAS_sims)

ggplot(df_vec_sens_perAS_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_sens_perAS_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of sensitivity of perturbation Adult survival")

# mean
mean_sens_perAS_sims <- mean(df_vec_sens_perAS_sims$value)
mean_sens_perAS_sims

# standard deviation
sd_sens_perAS_sims <- sd(df_vec_sens_perAS_sims$value)
sd_sens_perAS_sims

# SE (standard error) = sd / sqrt(n)
se_sens_perAS_sims <- sd_sens_perAS_sims / sqrt(length(df_vec_sens_perAS_sims$value))
se_sens_perAS_sims

# 95% CI
ci_lower_sens_perAS_sims <- mean_sens_perAS_sims - 1.96 * se_sens_perAS_sims
ci_lower_sens_perAS_sims
ci_upper_sens_perAS_sims <- mean_sens_perAS_sims + 1.96 * se_sens_perAS_sims
ci_upper_sens_perAS_sims

# Plot the density of elasticity of perturbation Adult survival
# Flatten into a numeric vector
vec_elas_perAS_sims <- unlist(elas_perAS_sims)

# make dataframe
df_vec_elas_perAS_sims <- data.frame(value = vec_elas_perAS_sims)

ggplot(df_vec_elas_perAS_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_elas_perAS_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of elasticity of perturbation Adult survival")

# mean
mean_elas_perAS_sims <- mean(df_vec_elas_perAS_sims$value)
mean_elas_perAS_sims

# standard deviation
sd_elas_perAS_sims <- sd(df_vec_elas_perAS_sims$value)
sd_elas_perAS_sims

# SE (standard error) = sd / sqrt(n)
se_elas_perAS_sims <- sd_elas_perAS_sims / sqrt(length(df_vec_elas_perAS_sims$value))
se_elas_perAS_sims

# 95% CI
ci_lower_elas_perAS_sims <- mean_elas_perAS_sims - 1.96 * se_elas_perAS_sims
ci_lower_elas_perAS_sims
ci_upper_elas_perAS_sims <- mean_elas_perAS_sims + 1.96 * se_elas_perAS_sims
ci_upper_elas_perAS_sims

## -----------------------------------------
# 6.1 Scenarios 3. Reduce 5% Adult recruitment
# Derivative matrix for f1yr_ad
derivA_f1yr_ad <- matrix(NA, nrow = 2, ncol = 2)
derivA_f1yr_ad[1, 1] <- 0
derivA_f1yr_ad[2, 1] <- 0 
derivA_f1yr_ad[2, 2] <- 0 
derivA_f1yr_ad[1, 2] <- 1 * per_fact

derivA_f1yr_ad

# Using "lapply" function for list object, then obtained the sens_results (1000 samples) 
# multiply by the derivative of matrix A for Adult recruitment (f1yr_ad) (1 matrix)

sens_perf_samples <- lapply(sens_results, function(x) derivA_f1yr_ad %*% x)
sens_perf_samples

# get 1000 values which extract only element [1,2] for perturbed adult recruitment
sens_perf_sims <- sapply(sens_perf_samples, function(x) x[1, 2])

length(sens_perf_sims)  # must be 1000

# This is sensitivity of perturbation adult recruitment from 1000 samples
sens_perf_sims

# 6.2 Estimating the elasticity for the perturbation Adult recruitment
# Initialize empty list to store results
elas_perf_sims <- vector("list", length = 1000)

# Loop through positions
for (i in seq_along(sens_perf_sims)) {
  elas_perf_sims[[i]] <- sens_perf_sims[[i]] * (result_f1[[i]] / lambda_values[i])
}

# This is elasticity of perturbation Adult recruitment from 1000 samples
elas_perf_sims

# Plot the density of sensitivity of perturbation Adult recruitment
# Flatten into a numeric vector
vec_sens_perf_sims <- unlist(sens_perf_sims)

# make dataframe
df_vec_sens_perf_sims <- data.frame(value = vec_sens_perf_sims)

ggplot(df_vec_sens_perf_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_sens_perf_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of sensitivity of perturbation Adult recruitment")

# mean
mean_sens_perf_sims <- mean(df_vec_sens_perf_sims$value)
mean_sens_perf_sims

# standard deviation
sd_sens_perf_sims <- sd(df_vec_sens_perf_sims$value)
sd_sens_perf_sims

# SE (standard error) = sd / sqrt(n)
se_sens_perf_sims <- sd_sens_perf_sims / sqrt(length(df_vec_sens_perf_sims$value))
se_sens_perf_sims

# 95% CI
ci_lower_sens_perf_sims <- mean_sens_perf_sims - 1.96 * se_sens_perf_sims
ci_lower_sens_perf_sims
ci_upper_sens_perf_sims <- mean_sens_perf_sims + 1.96 * se_sens_perf_sims
ci_upper_sens_perf_sims

# Plot the density of elasticity of perturbation Adult recruitment
# Flatten into a numeric vector
vec_elas_perf_sims <- unlist(elas_perf_sims)

# make dataframe
df_vec_elas_perf_sims <- data.frame(value = vec_elas_perf_sims)

ggplot(df_vec_elas_perf_sims, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
  xlim(-5, 5) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = mean(df_vec_elas_perf_sims$value), colour = "red", linetype = "dashed") +
  labs(title = "Density plot of elasticity of perturbation Adult recruitment")

# mean
mean_elas_perf_sims <- mean(df_vec_elas_perf_sims$value)
mean_elas_perf_sims

# standard deviation
sd_elas_perf_sims <- sd(df_vec_elas_perf_sims$value)
sd_elas_perf_sims

# SE (standard error) = sd / sqrt(n)
se_elas_perf_sims <- sd_elas_perf_sims / sqrt(length(df_vec_elas_perf_sims$value))
se_elas_perf_sims

# 95% CI
ci_lower_elas_perf_sims <- mean_elas_perf_sims - 1.96 * se_elas_perf_sims
ci_lower_elas_perf_sims
ci_upper_elas_perf_sims <- mean_elas_perf_sims + 1.96 * se_elas_perf_sims
ci_upper_elas_perf_sims

## Makes data frame for 3 perturbation sensitivity models
df_all_per_sens_models <- data.frame(
  value = c(vec_sens_perJS_sims, vec_sens_perAS_sims, vec_sens_perf_sims),
  model = rep(c("Model perJS", "Model perAS", "Model perf"), each = 1000)
)

df_all_per_sens_models

# calculate mean value of each model
means_per_sens <- df_all_per_sens_models %>%
  group_by(model) %>%
  summarise(mean_val = mean(value), .groups = "drop")

# Plot densities together with each mean line
ggplot(df_all_per_sens_models, aes(x = value, colour = model, fill = model)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  xlim(-3, 3) +
  scale_colour_grey(start = 0.2, end = 0.8) +   # line colours from dark to light grey
  scale_fill_grey(start = 0.8, end = 0.2) +     # fill shades reversed for contrast
  theme_minimal(base_size = 13) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_vline(data = means_per_sens, aes(xintercept = mean_val, linetype = model),
             colour = "black", linewidth = 0.1)+
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Density comparison of 3 perturbation sensitivity models",
       x = "Value", y = "Density")

## Makes data frame for 3 perturbation elasticity models
df_all_per_elas_models <- data.frame(
  value = c(vec_elas_perJS_sims, vec_elas_perAS_sims, vec_elas_perf_sims),
  model = rep(c("Model perJS", "Model perAS", "Model perf"), each = 1000)
)

df_all_per_elas_models

# calculate mean value of each model
means_per_elas <- df_all_per_elas_models %>%
  group_by(model) %>%
  summarise(mean_val = mean(value), .groups = "drop")

# Plot densities together with each mean line
ggplot(df_all_per_elas_models, aes(x = value, colour = model, fill = model)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  xlim(-3, 3) +
  scale_colour_grey(start = 0.2, end = 0.8) +   # line colours from dark to light grey
  scale_fill_grey(start = 0.8, end = 0.2) +     # fill shades reversed for contrast
  theme_minimal(base_size = 13) +
  theme(
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_vline(data = means_per_elas, aes(xintercept = mean_val, linetype = model),
             colour = "black", linewidth = 0.1)+
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Density comparison of 3 perturbation elasticity models",
       x = "Value", y = "Density")

## Make the table of vital rate sensitivity and elasticity for baseline models
sens_base_set <- data.frame(
  Sensitivity = c("Survival", "Emigration", "Recruitment"),
  Mean_value = c(mean_sens_S1yr_sims, mean_sens_E1yr_sims, mean_sens_f1yr_sims),
  ci_lower = c(ci_lower_sens_S1yr_sims, ci_lower_sens_E1yr_sims,     ci_lower_sens_f1yr_sims),
  ci_upper = c(ci_upper_sens_S1yr_sims, ci_upper_sens_E1yr_sims, ci_upper_sens_f1yr_sims)
)

sens_base_set

elas_base_set <- data.frame(
  Elasticity = c("Survival", "Emigration", "Recruitment"),
  Mean_value = c(mean_elas_S1yr_sims, mean_elas_E1yr_sims, mean_elas_f1yr_sims),
  ci_lower = c(ci_lower_elas_S1yr_sims, ci_lower_elas_E1yr_sims,     ci_lower_elas_f1yr_sims),
  ci_upper = c(ci_upper_elas_S1yr_sims, ci_upper_elas_E1yr_sims, ci_upper_elas_f1yr_sims)
)

elas_base_set

## Make the table of vital rate sensitivity and elasticity for perturbation models
sens_per_set <- data.frame(
  Sensitivity = c("Juvenile Survival", "Adult Survival", "Adult Recruitment"),
  Mean_value = c(mean_sens_perJS_sims, mean_sens_perAS_sims, mean_sens_perf_sims),
  ci_lower = c(ci_lower_sens_perJS_sims, ci_lower_sens_perAS_sims,     ci_lower_sens_perf_sims),
  ci_upper = c(ci_upper_sens_perJS_sims, ci_upper_sens_perAS_sims, ci_upper_sens_perf_sims)
)

sens_per_set

elas_per_set <- data.frame(
  Elasticity = c("Juvenile Survival", "Adult Survival", "Adult Recruitment"),
  Mean_value = c(mean_elas_perJS_sims, mean_elas_perAS_sims, mean_elas_perf_sims),
  ci_lower = c(ci_lower_elas_perJS_sims, ci_lower_elas_perAS_sims,     ci_lower_elas_perf_sims),
  ci_upper = c(ci_upper_elas_perJS_sims, ci_upper_elas_perAS_sims, ci_upper_elas_perf_sims)
)

elas_per_set
