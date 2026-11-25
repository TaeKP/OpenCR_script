# Appendix C1 - Population viability analysis: 
# Last check: 29 September 2026

rm(list=ls(all=TRUE))		# remove all objects

# required packages
# library(truncnorm)    

# The parameters were estimated using survey data collected every two years (2019, 2021 & 2023), 
# based on previous research (see Phosri et al., 2025)

# Parameters from open population Spatial Capture-Recapture
S <- 0.49 					# True survival (2 years interval)
f <- 0.76 					# Recruitment (2 years interval)
lambda <- 1.16 				# Population growth rate
initN <- 66 				# Initial population size; based on 275 km2
init_adultProp_mean <- 0.89 # Mean of Proportion of initial population that is adult; 3 years
init_adultProp_SD <- 0.05 	# SD of Proportion of initial population that is adult
E = abs(lambda - S - f) 	# Emigration rate
pR = 1 - E 					# Probability of remaining in study area
phi = S - E 				# Apparent survival  
#-------------------------------------------------------------------------------

# Implement uncertainty on all parameters
# Initial population
SD_pop <- SE_pop <- 8

# Growth rate (lambda)
SD_growth_rate <- SE_growth_rate <- 0.12

# Recruitment (f)
SD_recruitment <- SE_recruitment <- 0.12

# Survival rate
SD_survival_rate <- SE_survival_rate <- 0.091

#-------------------------------------------------------------------------------

# Parameters and uncertainty
initial_population <- initN     	# Initial population size
initialN_sd <- SD_pop            	# Obtained from the uncertainty calculation

growth_rate <- lambda          		# growth rate (lambda)
growth_rate_sd <- SD_growth_rate

recruitment_rate <- f          		# recruitment rate (f)
recruitment_rate_sd <- SD_recruitment

survival_rate <- S        			# True survival rate
survival_rate_sd <- SD_survival_rate

carrying_capacity <- 116     		# based on the suitable habitat and FC's home range size in KSRY
n_years <- 10                  		# Number of years to simulate
simulations <- 10000          		# Number of simulations

#-----------------------------------------------------------------------------------------------------
# Model_C1: Population growth rate only model for 2 years
#-----------------------------------------------------------------------------------------------------
# Function to simulate PVA: Growth rate + stochasticity

pva_simulation_G2ys <- function(initial_population, growth_rate, carrying_capacity, n_years) {
  population <- numeric(n_years)
  population[1] <- initial_population
  
  for (year in 2:n_years) {
    # Adding randomness
    stochastic_growthrate <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0) 
    population[year] <- population[year - 1] * stochastic_growthrate
    population[year] <- ifelse(population[year] > carrying_capacity, population[year-1], population[year])
    
    # Extinction event
    if (population[year] < 1) { 	
      population[year] <- 0
      break
    }
  }
  return(population)
}

# Run simulation
results_MC1 <- replicate(simulations, pva_simulation_G2ys(initial_population, growth_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_MC1 <- mean(results_MC1[n_years,] == 0))

#----------------------------------------- END of model C1 -------------------------------------------

#-----------------------------------------------------------------------------------------------------
# Model_C2: Incorporating vital rates model (Survival & Recruitment) for 2 years
#-----------------------------------------------------------------------------------------------------
# Function to simulate PVA: Growth rate + Survival + Recruitment + stochasticity

pva_simulation_GSR2ys <- function(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years) {
  N <- c(initN, rep(NA, n_years - 1)) # add N at time t1
  Surv <- Rec <- Emi <- rep(NA, n_years)
  
  # True survival (2 years interval)
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) 
  
  # Recruitment (2 years interval) 
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) 
  
  # Population growth rate  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) 
  
  # Emigration rate
  E = abs(lambda - S - f) 
  
  for(t in 3:length(N)){    
    N[t] <- N[t-2] * (S + f - E)
    Surv[t] <- N[t-2]*S
    Rec[t] <- N[t-2]*f
    Emi[t] <- N[t-2]*E
    N[t] <- ifelse(N[t] > carrying_capacity, N[t-2], N[t])
    
    # Extinction event
    for(j in 3:length(N)){
      if (N[j-2] < 1) { 
        N[j-2] <- 0
      }
      break
    }
  }
  return(N)
}

# Run simulation 
results_MC2 <- replicate(simulations, pva_simulation_GSR2ys(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_MC2 <- mean(results_MC2[n_years-1,] == 0))

#----------------------------------------- END of model C2 -------------------------------------------

# See Appendix C2 for the parameter derivative from the 2-year time step to develop the 1-year time step

#-----------------------------------------------------------------------------------------------------
# Model_C3: Incorporating vital rates model for 1 year
#-----------------------------------------------------------------------------------------------------
# Function to simulate PVA: Growth rate + Survival + Recruitment + stochasticity, 
# equivalent to the function "pva_simulation_GSR2ys()", but Model_C3 for annual vital rates

pva_simulation_GSR1yr <- function(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years) {
  N <- c(initN, rep(NA, n_years - 1))  # add N at time t1
  N_alt <- c(initN, rep(NA, n_years - 1))
  Surv_alt <- Rec_alt <- Emi_alt <- rep(NA, n_years)  
  
  # True survival (2 years interval)  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) 
  
  # Recruitment (2 years interval)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) 
  
  # Population growth rate
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) 
  
  # Emigration rate
  E = abs(lambda - S - f) 
  
  for(t in 2:length(N)){
    S1yr <- S / sqrt(lambda)
    f1yr <- f / sqrt(lambda)
    E1yr <- E / sqrt(lambda)
    N_alt[t] <- N_alt[t-1] * (S1yr + f1yr - E1yr)
    Surv_alt[t] <- N_alt[t-1]*S1yr
    Rec_alt[t] <- N_alt[t-1]*f1yr
    Emi_alt[t] <- N_alt[t-1]*E1yr
    N_alt[t] <- ifelse(N_alt[t] > carrying_capacity, N_alt[t-1], N_alt[t])
    
    # Extinction event
    if (N_alt[t] < 1) { 
      N_alt[t] <- 0
      break
    }
  }
  return(N_alt)
}

# Run simulation
results_MC3 <- replicate(simulations, pva_simulation_GSR1yr(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_MC3 <- mean(results_MC3[n_years,] == 0))
#----------------------------------------- END NOT RUN------------------------------------------------------------
