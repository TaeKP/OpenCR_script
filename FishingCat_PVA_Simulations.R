
# VITAL RATE & POPULATION PARAMETERS FROM SCR #
#---------------------------------------------#

S <- 0.49 # True survival (2 yrs)
f <- 0.76 # Recruitment (2 yrs)
lambda <- 1.16 # Population growth rate
initN <- 81 # Initial population size
#init_adultProp <- 0.4 # Proportion of initial population that is adult
init_adultProp_mean <- 0.89 # Mean of Proportion of initial population that is adult; 3 years
init_adultProp_SD <- 0.05 # SD of Proportion of initial population that is adult

E = abs(lambda - S - f) # Emigration rate
pR = 1 - E # Probability of remaining in study area
  
phi = S - E # Apparent survival  

#-------------------------------------------------------------------------------
### Implement uncertainty on all parameters in simplest model ###
# Originally code from Dr.Matt using "popbio" package
#-------------------------------------------------------------------------------
library(popbio)

# obtaining the results from OpenCR including mean, SE, 95%CI
# calculating the uncertainty (SD) from SE and sample size (n; number of observation)
# using the new sample size, due to the open model was used the number of minimum/maximum individual
# in the data set to calculate mean and SE, so in this case we used 34 or 61 to be a sample size.
# individual identification = 34(2019); 55(2021); 61(2023)
# follows SD = SE * sqrt(n)

# Initial population
# mean = 81
SE_pop <- 10
n <- 61

## we get: 
#(SD_pop <- SE_pop*sqrt(n))
SD_pop <- SE_pop

## Growth rate (lambda)
# mean = 1.16
SE_growth_rate <- 0.12

# we get: 
#(SD_growth_rate <- SE_growth_rate*sqrt(n))
SD_growth_rate <- SE_growth_rate

## Recruitment (f)
# mean = 0.76
SE_recruitment <- 0.12

# we get: 
#(SD_recruitment <- SE_recruitment*sqrt(n))
SD_recruitment <- SE_recruitment

## Survival rate
# mean = 0.49
SE_survival_rate <- 0.091

# we get: 
#(SD_survival_rate <- SE_survival_rate*sqrt(n))
SD_survival_rate <- SE_survival_rate

#-------------------------------------------------------------------------------
# Parameters and uncertainty that should be in simplest model

initial_population <- initN     # Initial population size
initialN_sd <- SD_pop            # obtained from the uncertainty calculation

growth_rate <- lambda          # growth rate (lambda)
growth_rate_sd <- SD_growth_rate

recruitment_rate <- f          # recruitment rate (f)
recruitment_rate_sd <- SD_recruitment

survival_rate <- S        # True survival rate
survival_rate_sd <- SD_survival_rate

carrying_capacity <- 140     # Carrying capacity of the environment; based on the suitable habitat and FC's home range size in KSRY
n_years <- 10                  # Number of years to simulate
simulations <- 1000          # Number of simulation runs; test

test_contantValues <- FALSE
if(test_contantValues){
  initialN_sd <- growth_rate_sd <- recruitment_rate_sd <- survival_rate_sd <- 0
}

# 2-YEAR POPULATION MODEL (Population growth rate only) #
#-------------------------------------------------------#

#-------------------------------------------------------------------------------
# Function to simulate growth_rate with stochasticity

source("Function_pva_simulation_2years_lambda.R")
#-------------------------------------------------------------------------------
# Running the simulation multiple times
results <- replicate(simulations, pva_simulation(initial_population, growth_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability <- mean(results[n_years,] == 0))



# 2-YEAR POPULATION MODEL (Vital rates) #
#---------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

N <- c(initN, rep(NA, n_years - 1)) # due to error in 1-year model; need to add object "N" outside the simulation function

pva_simulation_2yrs <- function(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years) {
  
  N <- c(initN, rep(NA, n_years - 1))
  Surv <- Rec <- Emi <- rep(NA, n_years)
  
  #N[1] <- truncnorm::rtruncnorm(1, mean = initN, sd = initialN_sd, a = 0)  # in case we need the random number of initN
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  E = abs(lambda - S - f) # Emigration rate
  
  for(t in 3:length(N)){
    
    N[t] <- N[t-2] * (S + f - E)
    Surv[t] <- N[t-2]*S
    Rec[t] <- N[t-2]*f
    Emi[t] <- N[t-2]*E
    
    #N[t] <- ifelse(N[t] > carrying_capacity, carrying_capacity, N[t])
    N[t] <- ifelse(N[t] > carrying_capacity, N[t-2], N[t])
    
    for(j in 3:length(N)){
      
      if (N[j-2] < 1) { # Extinction event
        N[j-2] <- 0
      }
      break
    }
  }
  return(N)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
results_2yrs <- replicate(simulations, pva_simulation_2yrs(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_2yrs <- mean(results_2yrs[n_years-1,] == 0))


#-------------------------------------------------------------------------------

# 1-YEAR POPULATION MODEL (vital rates) #
#---------------------------------------#

# Derivation:
# N[t+1] = N[t] * sqrt(lambda) 
#        = N[t] * sqrt(S + f - E)
# 
# N[t+2] = N[t]*lambda 
#        = N[t] * (S + f - E)
#        = N[t+1] * sqrt(lambda)
#        = N[t] * sqrt(S + f - E)
# 
# 
# Surv[t+2] = N[t] * S
#           = N[t+1] * S1yr
#           = (Surv[t+1] + Rec[t+1] - Emi[t+1]) * S1yr
#           = (N[t] * S1yr + N[t] * f1yr + N[t] * E1yr) * S1yr
#           = N[t] * (S1yr + f1yr + E1yr) * S1yr
# 
# --> S = (S1yr + f1yr + E1yr) * S1yr
# 
# 
# Rec[t+2] = N[t] * f
#          = N[t+1] * f1yr
# 
# --> f = (S1yr + f1yr + E1yr) * f1yr
# 
# 
# Emi[t+2] = N[t] * E
#          = N[t+1] * E1yr
# 
# --> E = (S1yr + f1yr + E1yr) * E1yr
# 
# 
# Assuming (S1yr + f1yr + E1yr) = sqrt(S + f - E) = sqrt(lamnda)
# 
# --> X = sqrt(lambda) * X1yr


#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_1yr <- function(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years) {
  
  N_alt <- c(initN, rep(NA, n_years - 1))
  Surv_alt <- Rec_alt <- Emi_alt <- rep(NA, n_years)
  
  #N_alt[1] <- truncnorm::rtruncnorm(1, mean = initN, sd = initialN_sd, a = 0)  # in case we need the random number of initN
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  E = abs(lambda - S - f) # Emigration rate
  
  for(t in 2:length(N)){
    
    S1yr <- S / sqrt(lambda)
    f1yr <- f / sqrt(lambda)
    E1yr <- E / sqrt(lambda)
    
    N_alt[t] <- N_alt[t-1] * (S1yr + f1yr - E1yr)
    Surv_alt[t] <- N_alt[t-1]*S1yr
    Rec_alt[t] <- N_alt[t-1]*f1yr
    Emi_alt[t] <- N_alt[t-1]*E1yr
    
    #N_alt[t] <- ifelse(N_alt[t] > carrying_capacity, carrying_capacity, N_alt[t])
    N_alt[t] <- ifelse(N_alt[t] > carrying_capacity, N_alt[t-1], N_alt[t])
    
    if (N_alt[t] < 1) { # Extinction event
      N_alt[t] <- 0
      break
    }
  }
  return(N_alt)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
results_1yr <- replicate(simulations, pva_simulation_1yr(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_1yr <- mean(results_1yr[n_years,] == 0))


# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_age_str <- function(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  for(t in 2:length(N)){
    
    # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
    A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
    A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
    A[1, 2] <- f1yr # Adults producing juveniles within 1 year
    
    # Calculate "per adult recruitment rate" for time-step
    current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
    f1yr_ad <- f1yr /  current_adultProp
    A[1, 2] <- f1yr_ad
    
    # Project
    N_mat2[, t] <- A %*% N_mat2[, t-1]
    
    if(sum(N_mat2[, t]) > carrying_capacity){
      N_mat2[, t] <- N_mat2[, t-1]
    }
    
    if (N_mat2[t] < 1) { # Extinction event
      N_mat2[t] <- 0
      break
    }
  }
  return(N_mat2)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str <- replicate(simulations, pva_simulation_age_str(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years)) )


# VISUAL COMPARISON #
#-------------------#

library(magrittr)
library(ggplot2)

# Write results as data frames

sim_1 <- reshape2::melt(results) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "2-year, lambda only",
                Year = (Year*2)-1) %>%
  dplyr::filter(Year <= n_years)

sim_2 <- reshape2::melt(results_2yrs) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "2-year, vital rates") %>%
  dplyr::filter((Year %% 2) != 0) # Drop even years that are "skipped" by model

sim_3 <- reshape2::melt(results_1yr) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates")

sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary <- rbind(sim_1, sim_2, sim_3, sim_4) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 1	
## Decreasing the survival rate from 0.49 to 0.29 (Fixed other values)
#S <- 0.49 # True survival (2 yrs)
S <- 0.29 # True survival (2 yrs) # scenario 1
survival_rate <- S        # True survival rate
#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".
source("Function_pva_simulation_age_str.R")


#-------------------------------------------------------------------------------
# Running the simulation multiple times

# Baseline scenario
(results_age_str_s1 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity,
                                                                     n_years)) )

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s1 <- reshape2::melt(apply(results_age_str_s1, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 1")

## Combine results and summarise
simSummary_s1 <- rbind(sim_s1) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s1, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 2	
## Decreasing the recruitment rate from 0.76 to 0.56 (Fixed other values)
S <- 0.49 # True survival (2 yrs)
#f <- 0.76 # Recruitment (2 yrs)  # mean value
f <- 0.56 # Recruitment (2 yrs)   # scenario 2
recruitment_rate <- f          # recruitment rate (f)
survival_rate <- S        # True survival rate
#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_age_str_s2 <- function(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  for(t in 2:length(N)){
    
    # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
    A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
    A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
    A[1, 2] <- f1yr # Adults producing juveniles within 1 year
    
    # Calculate "per adult recruitment rate" for time-step
    current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
    f1yr_ad <- f1yr /  current_adultProp
    A[1, 2] <- f1yr_ad
    
    # Project
    N_mat2[, t] <- A %*% N_mat2[, t-1]
    
    if(sum(N_mat2[, t]) > carrying_capacity){
      N_mat2[, t] <- N_mat2[, t-1]
    }
    
    if (N_mat2[t] < 1) { # Extinction event
      N_mat2[t] <- 0
      break
    }
  }
  return(N_mat2)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str_s2 <- replicate(simulations, pva_simulation_age_str(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years)) )

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s2 <- reshape2::melt(apply(results_age_str_s2, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 2")

## Combine results and summarise
simSummary_s2 <- rbind(sim_s2) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s2, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 3	
## Decreasing the growth rate from 1.16 to 0.96 (Fixed other values)
f <- 0.76 # Recruitment (2 yrs)   # mean value
#lambda <- 1.16                   # mean value
lambda <- 0.96 # Population growth rate  # scenario 3
growth_rate <- lambda          # growth rate (lambda)
recruitment_rate <- f          # recruitment rate (f)
#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_age_str_s3 <- function(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  for(t in 2:length(N)){
    
    # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
    A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
    A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
    A[1, 2] <- f1yr # Adults producing juveniles within 1 year
    
    # Calculate "per adult recruitment rate" for time-step
    current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
    f1yr_ad <- f1yr /  current_adultProp
    A[1, 2] <- f1yr_ad
    
    # Project
    N_mat2[, t] <- A %*% N_mat2[, t-1]
    
    if(sum(N_mat2[, t]) > carrying_capacity){
      N_mat2[, t] <- N_mat2[, t-1]
    }
    
    if (N_mat2[t] < 1) { # Extinction event
      N_mat2[t] <- 0
      break
    }
  }
  return(N_mat2)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str_s3 <- replicate(simulations, pva_simulation_age_str(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years)) )

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s3 <- reshape2::melt(apply(results_age_str_s3, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 3")

## Combine results and summarise
simSummary_s3 <- rbind(sim_s3) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s3, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 4	
## Decreasing the adult proportion from 0.92 to 0.72 (Fixed other values)
lambda <- 1.16                    # mean value
#init_adultProp_mean <- 0.92      # mean value
init_adultProp_mean <- 0.72  # scenario 4
growth_rate <- lambda          # growth rate (lambda)
#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_age_str_s4 <- function(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  for(t in 2:length(N)){
    
    # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
    A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
    A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
    A[1, 2] <- f1yr # Adults producing juveniles within 1 year
    
    # Calculate "per adult recruitment rate" for time-step
    current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
    f1yr_ad <- f1yr /  current_adultProp
    A[1, 2] <- f1yr_ad
    
    # Project
    N_mat2[, t] <- A %*% N_mat2[, t-1]
    
    if(sum(N_mat2[, t]) > carrying_capacity){
      N_mat2[, t] <- N_mat2[, t-1]
    }
    
    if (N_mat2[t] < 1) { # Extinction event
      N_mat2[t] <- 0
      break
    }
  }
  return(N_mat2)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str_s4 <- replicate(simulations, pva_simulation_age_str(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years)) )

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s4 <- reshape2::melt(apply(results_age_str_s4, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 4")

## Combine results and summarise
simSummary_s4 <- rbind(sim_s4) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s4, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
# VISUAL COMPARISON #
# 4 scenarios (1-4)
#-------------------#

## Combine results and summarise
simSummary_4models <- rbind(sim_s1, sim_s2, sim_s3, sim_s4) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_4models, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()
#-------------------------------------------------------------------------------

## Increasing vital rate due to high management context and low treats in the area. 
## Which increase 20% of each vital rate
#-------------------------------------------------------------------------------
## Scenario 5	
## Increasing the survival rate

#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

(results_age_str_s5 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity,
                                                                     pertFac.S = 1.2,
                                                                     n_years)) )

#-------------------------------------------------------------------------------
# Running the simulation multiple times

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s5 <- reshape2::melt(apply(results_age_str_s5, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 5")

## Combine results and summarise
simSummary_s5 <- rbind(sim_s5) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s5, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 6	
## Increasing the recruitment rate from 0.76 to 0.96 (Fixed other values)
S <- 0.49 # True survival (2 yrs)
#f <- 0.76 # Recruitment (2 yrs)  # mean value
f <- 0.96 # Recruitment (2 yrs)   # scenario 6
recruitment_rate <- f          # recruitment rate (f)
survival_rate <- S        # True survival rate
#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_age_str_s6 <- function(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  for(t in 2:length(N)){
    
    # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
    A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
    A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
    A[1, 2] <- f1yr # Adults producing juveniles within 1 year
    
    # Calculate "per adult recruitment rate" for time-step
    current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
    f1yr_ad <- f1yr /  current_adultProp
    A[1, 2] <- f1yr_ad
    
    # Project
    N_mat2[, t] <- A %*% N_mat2[, t-1]
    
    if(sum(N_mat2[, t]) > carrying_capacity){
      N_mat2[, t] <- N_mat2[, t-1]
    }
    
    if (N_mat2[t] < 1) { # Extinction event
      N_mat2[t] <- 0
      break
    }
  }
  return(N_mat2)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str_s6 <- replicate(simulations, pva_simulation_age_str(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years)) )

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s6 <- reshape2::melt(apply(results_age_str_s6, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 6")

## Combine results and summarise
simSummary_s6 <- rbind(sim_s6) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s6, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 7	
## Increasing the growth rate from 1.16 to 1.36 (Fixed other values)
f <- 0.76 # Recruitment (2 yrs)   # mean value
#lambda <- 1.16                   # mean value
lambda <- 1.36 # Population growth rate  # scenario 7
growth_rate <- lambda          # growth rate (lambda)
recruitment_rate <- f          # recruitment rate (f)
#-------------------------------------------------------------------------------
# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_age_str_s7 <- function(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  for(t in 2:length(N)){
    
    # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
    A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
    A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
    A[1, 2] <- f1yr # Adults producing juveniles within 1 year
    
    # Calculate "per adult recruitment rate" for time-step
    current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
    f1yr_ad <- f1yr /  current_adultProp
    A[1, 2] <- f1yr_ad
    
    # Project
    N_mat2[, t] <- A %*% N_mat2[, t-1]
    
    if(sum(N_mat2[, t]) > carrying_capacity){
      N_mat2[, t] <- N_mat2[, t-1]
    }
    
    if (N_mat2[t] < 1) { # Extinction event
      N_mat2[t] <- 0
      break
    }
  }
  return(N_mat2)
}

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str_s7 <- replicate(simulations, pva_simulation_age_str(initN, growth_rate, survival_rate, recruitment_rate, init_adultProp, carrying_capacity, n_years)) )

#-------------------------------------------------------------------------------

# Write results as data frames

sim_s7 <- reshape2::melt(apply(results_age_str_s7, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "scenario 7")

## Combine results and summarise
simSummary_s7 <- rbind(sim_s7) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s7, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
# VISUAL COMPARISON #
# 3 scenarios (5-7)
#-------------------#

## Combine results and summarise
simSummary_3models <- rbind(sim_s5, sim_s6, sim_s7) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_3models, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()
