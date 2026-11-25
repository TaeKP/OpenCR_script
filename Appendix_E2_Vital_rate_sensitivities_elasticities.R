# Appendix E2 - Asymptotic sensitivity analysis code: Vital rate sensitivities and elasticities 
# Last check: 1 June 2026

rm(list=ls(all=TRUE))		# remove all objects

# package require
library(ggplot2)
library(dplyr)

# Set seed
mySeed <- 90
set.seed(mySeed)

#-----------------------------------------------------------------------------------------------------
# Vital rate sensitivities and elasticities
#-----------------------------------------------------------------------------------------------------
# Derivative of lambda with respect to each vital rate
# potential effect a small absolute change in vital rate could have on population growth rate
# *Obtained the lambda value from Appendix E1

# load objects from Appendix E1
load("appendix_E1.RData")

# Projection matrix
#A[1, 1] <- 0 
#A[2, 1] <- S1yr - E1yr 
#A[2, 2] <- S1yr - E1yr 
#A[1, 2] <- f1yr_ad

# Derivative matrix for survival (S1yr)
derivA_S1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_S1yr[1, 1] <- 0 
derivA_S1yr[2, 1] <- 1 
derivA_S1yr[2, 2] <- 1 
derivA_S1yr[1, 2] <- 0

# Derivative matrix for emigration (E1yr)
derivA_E1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_E1yr[1, 1] <- 0 
derivA_E1yr[2, 1] <- -1 
derivA_E1yr[2, 2] <- -1
derivA_E1yr[1, 2] <- 0

# Derivative matrix for recruitment (f1yr_ad)
derivA_f1yr <- matrix(NA, nrow = 2, ncol = 2)
derivA_f1yr[1, 1] <- 0 
derivA_f1yr[2, 1] <- 0 
derivA_f1yr[2, 2] <- 0 
derivA_f1yr[1, 2] <- 1

#-----------------------------------------------------------------------------------------------------
# Model_1: Vital rate sensitivity and elasticity for survival
#-----------------------------------------------------------------------------------------------------
## 1.1 Vital rate sensitivity for survival
# Using "lapply" function for list object, then obtained the sens_results (10000 samples) 
# multiply by the derivative of matrix A for survival (S1yr) (1 matrix) 
# and using sum function to calculate each value in each slice

# Multiply each matrix in the list by derivA_S1yr
# get 10000 values of sens_S1yr_sims from sum function

sens_S1yr_sims <- lapply(sens_results, function(x) sum(derivA_S1yr %*% x))

# Sensitivity of survival from 10000 samples
View(sens_S1yr_sims)

## 1.2 Estimating the elasticity for survival
# building the object for the S1yr values
result_S1 <- c()
for(i in 1:nsim){
  result_S1 <- c(result_S1, nsim_VRs[[i]]$S1yr)
}

# Checking the structure of result_S1 for list object include 10000 values
str(result_S1)

# Initialize empty list to store results
elas_S1yr_sims <- vector("list", length = 10000)

# Loop through positions
for (i in seq_along(sens_S1yr_sims)) {
  elas_S1yr_sims[[i]] <- sens_S1yr_sims[[i]] * (result_S1[i] / lambda_values[i])
}

# Elasticity of survival from 10000 samples
View(elas_S1yr_sims)
#----------------------------------------- END of vital rate sensitivity for survival --------------------------------

## Plot the density of sensitivity of survival
# Flatten into a numeric vector
vec_sens_S1yr_sims <- unlist(sens_S1yr_sims)

# Make dataframe
df_vec_sens_S1yr_sims <- data.frame(value = vec_sens_S1yr_sims)

ggplot(df_vec_sens_S1yr_sims, aes(x = value)) +
      geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
      xlim(-5, 5) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = median(df_vec_sens_S1yr_sims$value), colour = "red", linetype = "dashed") +
      labs(title = "Sensitivity of pop. growth rate (survival)")

# Calculating statistic results
# median
(median_sens_S1yr_sims <- median(df_vec_sens_S1yr_sims$value))

# 95% CI
(ci_lower_sens_S1yr_sims <- quantile(df_vec_sens_S1yr_sims$value, probs = 0.025))
(ci_upper_sens_S1yr_sims <- quantile(df_vec_sens_S1yr_sims$value, probs = 0.975))

## Plot the density of elasticity of survival
# Flatten into a numeric vector
vec_elas_S1yr_sims <- unlist(elas_S1yr_sims)

# make dataframe
df_vec_elas_S1yr_sims <- data.frame(value = vec_elas_S1yr_sims)

ggplot(df_vec_elas_S1yr_sims, aes(x = value)) +
      geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
      xlim(-5, 5) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = median(df_vec_elas_S1yr_sims$value), colour = "red", linetype = "dashed") +
      labs(title = "Elasticity of pop. growth rate (survival)")

# Calculating statistic results
# median
(median_elas_S1yr_sims <- median(df_vec_elas_S1yr_sims$value))

# 95% CI
(ci_lower_elas_S1yr_sims <- quantile(df_vec_elas_S1yr_sims$value, probs = 0.025))
(ci_upper_elas_S1yr_sims <- quantile(df_vec_elas_S1yr_sims$value, probs = 0.975))

#----------------------------------------- END of vital rate sensitivity for survival plot ---------------------------
#-----------------------------------------------------------------------------------------------------
# Model_2: Vital rate sensitivity and elasticity for emigration
#-----------------------------------------------------------------------------------------------------
## 2.1 Vital rate sensitivities for emigration
# Using "lapply" function for list object, then obtained the sens_results (10000 samples) 
# multiply by the derivative of matrix A for emigration (E1yr) (1 matrix)

# Multiply each matrix in the list by derivA_E1yr
# get 10000 values of sens_E1yr_sims from sum function
sens_E1yr_sims <- lapply(sens_results, function(x) sum(derivA_E1yr %*% x))

# Sensitivity of emigration from 10000 samples
View(sens_E1yr_sims)

## 2.2 Estimating the elasticity for emigration
# Building the object for the E1yr values
result_E1 <- c()
for(i in 1:nsim){
  result_E1 <- c(result_E1, nsim_VRs[[i]]$E1yr)
}

# Checking the structure of result_E1 for list object include 10000 values
str(result_E1)

# Initialize empty list to store results
elas_E1yr_sims <- vector("list", length = 10000)

# Loop through positions
for (i in seq_along(sens_E1yr_sims)) {
  elas_E1yr_sims[[i]] <- sens_E1yr_sims[[i]] * (result_E1[[i]] / lambda_values[i])
}

# Elasticity of emigration from 10000 samples
View(elas_E1yr_sims)
#----------------------------------------- END of vital rate sensitivity for emigration ------------------------------

## Plot the density of sensitivity of emigration
# Flatten into a numeric vector
vec_sens_E1yr_sims <- unlist(sens_E1yr_sims)

# make dataframe
df_vec_sens_E1yr_sims <- data.frame(value = vec_sens_E1yr_sims)

ggplot(df_vec_sens_E1yr_sims, aes(x = value)) +
      geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
      xlim(-5, 5) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = median(df_vec_sens_E1yr_sims$value), colour = "red", linetype = "dashed") +
      labs(title = "Sensitivity of pop. growth rate (emigration)")

# Calculating statistic results
# median
(median_sens_E1yr_sims <- median(df_vec_sens_E1yr_sims$value))

# 95% CI
(ci_lower_sens_E1yr_sims <- quantile(df_vec_sens_E1yr_sims$value, probs = 0.025))
(ci_upper_sens_E1yr_sims <- quantile(df_vec_sens_E1yr_sims$value, probs = 0.975))

## Plot the density of elasticity of emigration
# Flatten into a numeric vector
vec_elas_E1yr_sims <- unlist(elas_E1yr_sims)

# make dataframe
df_vec_elas_E1yr_sims <- data.frame(value = vec_elas_E1yr_sims)

ggplot(df_vec_elas_E1yr_sims, aes(x = value)) +
      geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
      xlim(-5, 5) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = median(df_vec_elas_E1yr_sims$value), colour = "red", linetype = "dashed") +
      labs(title = "Elasticity of pop. growth rate (emigration)")

# Calculating statistic results
# median
(median_elas_E1yr_sims <- median(df_vec_elas_E1yr_sims$value))

# 95% CI
(ci_lower_elas_E1yr_sims <- quantile(df_vec_elas_E1yr_sims$value, probs = 0.025))
(ci_upper_elas_E1yr_sims <- quantile(df_vec_elas_E1yr_sims$value, probs = 0.975))

#----------------------------------------- END of vital rate sensitivity for emigration plot -------------------------
#-----------------------------------------------------------------------------------------------------
# Model_3: Vital rate sensitivity and elasticity for recruitment
#-----------------------------------------------------------------------------------------------------
## 3.1 Vital rate sensitivities for recruitment
# Using "lapply" function for list object, then obtained the sens_results (10000 samples) 
# multiply by the derivative of matrix A for recruitment (f1yr_ad) (1 matrix)

# Multiply each matrix in the list by derivA_f1yr
sens_f_samples <- lapply(sens_results, function(x) derivA_f1yr %*% x)
View(sens_f_samples)

# get 10000 values which extract only element [1,2] for adult recruitment
sens_f1yr_sims <- sapply(sens_f_samples, function(x) x[1, 2])

length(sens_f1yr_sims)  # must be 10000

# Sensitivity of recruitment from 10000 samples
View(sens_f1yr_sims)

## 3.2 Estimating the elasticity for recruitment
# Building the object for the f1yr_ad values
result_f1 <- c()
for(i in 1:nsim){
  result_f1 <- c(result_f1, nsim_VRs[[i]]$f1yr)
}

# Checking the structure of result_f1 for list object include 10000 values
str(result_f1)

# Initialize empty list to store results
elas_f1yr_sims <- vector("list", length = 10000)

# Loop through positions
for (i in seq_along(sens_f1yr_sims)) {
  elas_f1yr_sims[[i]] <- sens_f1yr_sims[[i]] * (result_f1[[i]] / lambda_values[i])
}

# Elasticity of recruitment from 10000 samples
View(elas_f1yr_sims)
#----------------------------------------- END of vital rate sensitivity for recruitment -----------------------------

## Plot the density of sensitivity of recruitment
# Flatten into a numeric vector
vec_sens_f1yr_sims <- unlist(sens_f1yr_sims)

# make dataframe
df_vec_sens_f1yr_sims <- data.frame(value = vec_sens_f1yr_sims)

ggplot(df_vec_sens_f1yr_sims, aes(x = value)) +
      geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
      xlim(-5, 5) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = median(df_vec_sens_f1yr_sims$value), colour = "red", linetype = "dashed") +
      labs(title = "Sensitivity of pop. growth rate (recruitment)")

# Calculating statistic results
# median
(median_sens_f1yr_sims <- median(df_vec_sens_f1yr_sims$value))

# 95% CI
(ci_lower_sens_f1yr_sims <- quantile(df_vec_sens_f1yr_sims$value, probs = 0.025))
(ci_upper_sens_f1yr_sims <- quantile(df_vec_sens_f1yr_sims$value, probs = 0.975))

## Plot the density of elasticity of recruitment
# Flatten into a numeric vector
vec_elas_f1yr_sims <- unlist(elas_f1yr_sims)

# make dataframe
df_vec_elas_f1yr_sims <- data.frame(value = vec_elas_f1yr_sims)

ggplot(df_vec_elas_f1yr_sims, aes(x = value)) +
      geom_density(fill = "skyblue", alpha = 0.4, colour = "blue") +
      xlim(-5, 5) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = median(df_vec_elas_f1yr_sims$value), colour = "red", linetype = "dashed") +
      labs(title = "Elasticity of pop. growth rate (recruitment)")

# Calculating statistic results
# median
(median_elas_f1yr_sims <- median(df_vec_elas_f1yr_sims$value))

# 95% CI
(ci_lower_elas_f1yr_sims <- quantile(df_vec_elas_f1yr_sims$value, probs = 0.025))
(ci_upper_elas_f1yr_sims <- quantile(df_vec_elas_f1yr_sims$value, probs = 0.975))
#----------------------------------------- END of vital rate sensitivity for recruitment plot ------------------------

## Compare the 3 sensitivity models
# makes data frame for 3 sensitivity models
df_all_sens_models <- data.frame(
                      value = c(vec_sens_S1yr_sims, vec_sens_E1yr_sims, vec_sens_f1yr_sims),
                      Model = rep(c("Survival", "Emigration", "Recruitment"), each = 10000))

View(df_all_sens_models)

# Calculate the median value for each model
median_sens <- df_all_sens_models %>%
              group_by(Model) %>%
              summarise(median_val = median(value), .groups = "drop")

# Plot densities together with each median line
ggplot(df_all_sens_models, aes(x = value, colour = Model, fill = Model)) +
      geom_density(alpha = 0.3, linewidth = 0.5) + xlim(-3, 3) +
      scale_colour_manual(values = c("Survival" = "black", 
                                    "Emigration" = "dimgray", 
                                    "Recruitment"= "darkgray")) +
      scale_fill_manual(values = c("Survival" = "gray10", 
                                  "Emigration" = "gray90", 
                                  "Recruitment"= "gray50")) + 
      theme_minimal(base_size = 13) +
      theme(axis.line = element_line(color = "black")) +
      geom_hline(yintercept = 0, colour = "black") +
      geom_vline(data = median_sens, aes(xintercept = median_val, linetype = Model), 
                 colour = "black", linewidth = 0.3) +
      scale_y_continuous(expand = c(0,0)) +
      labs(title = "Comparison of sensitivity models", x = "Value", y = "Density")
#----------------------------------------- END of 3 sensitivity models comparison ------------------------------------

# Makes data frame for 3 elasticity models
df_all_elas_models <- data.frame(
                      value = c(vec_elas_S1yr_sims, vec_elas_E1yr_sims, vec_elas_f1yr_sims),
                      Model = rep(c("Survival", "Emigration", "Recruitment"), each = 10000))

View(df_all_elas_models)

# Calculate the median value for each model
median_elas <- df_all_elas_models %>%
              group_by(Model) %>%
              summarise(median_val = median(value), .groups = "drop")

# Plot densities together with each median line
ggplot(df_all_elas_models, aes(x = value, colour = Model, fill = Model)) +
      geom_density(alpha = 0.3, linewidth = 0.5) + 
      xlim(-3, 3) +
      scale_colour_manual(values = c("Survival" = "black", 
                                    "Emigration" = "dimgray", 
                                    "Recruitment"= "darkgray")) + 
      scale_fill_manual(values = c("Survival" = "gray10", 
                                  "Emigration" = "gray90", 
                                  "Recruitment"= "gray50")) +
      theme_minimal(base_size = 13) +
      theme(axis.line = element_line(color = "black")) + 
      geom_vline(data = median_elas, aes(xintercept = median_val, linetype = Model),
             colour = "black", linewidth = 0.3) + 
      geom_hline(yintercept = 0, colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      labs(title = "Comparison of elasticity models", x = "Value", y = "Density")

#----------------------------------------- END of 3 elasticity models comparison -------------------------------------

