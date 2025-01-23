
# VITAL RATE & POPULATION PARAMETERS FROM SCR #
#---------------------------------------------#

S <- 0.49 # True survival (2 yrs)
f <- 0.76 # Recruitment (2 yrs)
lambda <- 1.16 # Population growth rate
initN <- 81 # Initial population size
init_adultProp <- 0.4 # Proportion of initial population that is adult

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

initial_population <- 81     # Initial population size
initialN_sd <- SD_pop            # obtained from the uncertainty calculation

growth_rate <- 1.16          # growth rate (lambda)
growth_rate_sd <- SD_growth_rate

recruitment_rate <- 0.76          # recruitment rate (f)
recruitment_rate_sd <- SD_recruitment

survival_rate <- 0.49        # True survival rate
survival_rate_sd <- SD_survival_rate

carrying_capacity <- 140     # Carrying capacity of the environment; based on the suitable habitat and FC's home range size in KSRY
years <- 10                  # Number of years to simulate
simulations <- 1000          # Number of simulation runs; test
#-------------------------------------------------------------------------------
# Function to simulate growth_rate with stochasticity

pva_simulation <- function(initial_population, growth_rate, carrying_capacity, years) {
  population <- numeric(years)
  population[1] <- truncnorm::rtruncnorm(1, mean = initial_population, sd = initialN_sd, a = 0) # we probably want this lognormal (or truncated at 0)
  
  for (year in 2:years) {
    
    stochastic_growthrate <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0) # Adding randomness
    
    population[year] <- population[year - 1] * stochastic_growthrate
    
    population[year] <- ifelse(population[year] > carrying_capacity, carrying_capacity, population[year])
    
    if (population[year] < 1) { # Extinction event
      population[year] <- 0
      break
    }
  }
  return(population)
}
#-------------------------------------------------------------------------------
# Running the simulation multiple times
results <- replicate(simulations, pva_simulation(initial_population, growth_rate, carrying_capacity, years))
#results <- replicate(simulations, pva_simulation(initial_population, survival_rate, carrying_capacity, years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability <- mean(results[years,] == 0))

# Plotting results
matplot(1:years, results, type = "l", lty = 1, col = rgb(0.1, 0.1, 0.1, 0.1),
        xlab = "Year", ylab = "Population Size", main = "PVA Simulation")
abline(h = carrying_capacity, col = "red", lty = 2)

cat("Estimated Probability of Extinction:", extinction_probability, "\n")

#-------------------------------------------------------------------------------
# Now the uncertainty was implemented in the growth_rate in the simple model.
#-------------------------------------------------------------------------------

# 2-YEAR POPULATION MODEL #
#-------------------------#

## implement uncertainty to all parameters for the 2-YEAR POPULATION MODEL
S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0) # Population growth rate
initN <- initial_population # Initial population size; based on open model (2023)
E = abs(lambda - S - f) # Emigration rate
#-------------------------------------------------------------------------------

n_years <- 10
N <- c(initN, rep(NA, n_years - 1))
Surv <- Rec <- Emi <- rep(NA, n_years)

for(t in 3:length(N)){
  
  N[t] <- N[t-2] * (S + f - E)
  Surv[t] <- N[t-2]*S
  Rec[t] <- N[t-2]*f
  Emi[t] <- N[t-2]*E
}

cbind(N, Surv, Rec, Emi, Surv + Rec - Emi)

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".

pva_simulation_2yrs <- function(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years) {
  
  N <- c(initN, rep(NA, n_years - 1))
  Surv <- Rec <- Emi <- rep(NA, n_years)
  
  #N[1] <- truncnorm::rtruncnorm(1, mean = initN, sd = initialN_sd, a = 0)  # in case we need the random number of initN
  
  for(t in 3:length(N)){
    
    S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
    f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
    lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0) # Population growth rate
    E = abs(lambda - S - f) # Emigration rate
    
    N[t] <- N[t-2] * (S + f - E)
    Surv[t] <- N[t-2]*S
    Rec[t] <- N[t-2]*f
    Emi[t] <- N[t-2]*E
    
    N[t] <- ifelse(N[t] > carrying_capacity, carrying_capacity, N[t])
    
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

View(results_2yrs)

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_2yr <- mean(results_2yrs[n_years-1,] == 0))

#-------------------------------------------------------------------------------
## Alternative way to simulate data
# adding the loop of simulation for 2 year model
sim_2yrs <- sapply(1:simulations, function(x) {
  
  n_years <- 10
  N <- c(initN, rep(NA, n_years - 1))
  Surv <- Rec <- Emi <- rep(NA, n_years)
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) 
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) 
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0) 
  initN <- initial_population 
  E = abs(lambda - S - f) 
  
  for(t in 3:length(N)){
    
    N[t] <- N[t-2] * (S + f - E)
    Surv[t] <- N[t-2]*S
    Rec[t] <- N[t-2]*f
    Emi[t] <- N[t-2]*E
  }
  N
})

View(sim_2yrs)

#-------------------------------------------------------------------------------

# 1-YEAR POPULATION MODEL #
#-------------------------#

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
S1yr <- S / sqrt(lambda)
f1yr <- f / sqrt(lambda)
E1yr <- E / sqrt(lambda)


N_alt <- c(initN, rep(NA, n_years - 1))
Surv_alt <- Rec_alt <- Emi_alt <- rep(NA, n_years)


for(t in 2:length(N)){
  
  N_alt[t] <- N_alt[t-1] * (S1yr + f1yr - E1yr)
  Surv_alt[t] <- N_alt[t-1]*S1yr
  Rec_alt[t] <- N_alt[t-1]*f1yr
  Emi_alt[t] <- N_alt[t-1]*E1yr
}

cbind(N, Surv, Rec, Emi, Surv + Rec - Emi)
cbind(N_alt, Surv_alt, Rec_alt, Emi_alt, Surv_alt + Rec_alt - Emi_alt)

# --> These are now equivalent 2-year and 1-year interval models
#-------------------------------------------------------------------------------
# adding the loop of simulation for 1 year model
sim_1yr <- sapply(1:simulations, function(x) {
  
  S1yr <- S / sqrt(lambda)
  f1yr <- f / sqrt(lambda)
  E1yr <- E / sqrt(lambda)
  
  N_alt <- c(initN, rep(NA, n_years - 1))
  Surv_alt <- Rec_alt <- Emi_alt <- rep(NA, n_years)
  
  for(t in 2:length(N)){
    
    N_alt[t] <- N_alt[t-1] * (S1yr + f1yr - E1yr)
    Surv_alt[t] <- N_alt[t-1]*S1yr
    Rec_alt[t] <- N_alt[t-1]*f1yr
    Emi_alt[t] <- N_alt[t-1]*E1yr
  }
  N_alt
})

View(sim_1yr)
#-------------------------------------------------------------------------------

# 1-YEAR POPULATION WITH 2 AGE CLASSES #
#--------------------------------------#

# Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
A <- matrix(NA, nrow = 2, ncol = 2)

A[1, 1] <- 0 # Juveniles producing juveniles within 1 year
A[2, 1] <- S1yr - E1yr # Juvenile becoming adults within 1 year
A[2, 2] <- S1yr - E1yr # Adults remaining alive (& in study area) within 1 year
A[1, 2] <- f1yr # Adults producing juveniles within 1 year

N_mat <- matrix(NA, nrow = 2, ncol = n_years)
N_mat[,1] <- c(1 - init_adultProp, init_adultProp)*initN

for(t in 2:length(N)){
  N_mat[, t] <- A %*% N_mat[, t-1]
}

N_mat
cbind(colSums(N_mat), N_alt) 

# --> These are not equivalent
# --> The reason for that is that f1yr (and f) are defined "per capita" for
#     the whole population, not just for the adults. If we set A[1, 1] = f1yr
#     as well, we obtain an equivalent model. But that is not really what we
#     want to do here. 
# --> Instead, we will try to "re-calibrate" f1yr to correspond to recruitment
#     rate per adult as opposed to recruitment rate per individual

N_mat <- matrix(NA, nrow = 2, ncol = n_years)
N_mat[,1] <- c(1 - init_adultProp, init_adultProp)*initN

for(t in 2:length(N)){
  
  # Calculate "per adult recruitment rate" for time-step
  f1yr_ad <- f1yr / init_adultProp 
  A[1, 2] <- f1yr_ad
  
  # Project
  N_mat[, t] <- A %*% N_mat[, t-1]
}

N_mat
cbind(colSums(N_mat), N_alt) 

# --> Now, the first time step projects identical to the non-structured model
# --> After that, it is different again though. The reason for that is that
#     population structure (juvenile:adult ratio) does change, but we calculated
#     f1yr_ad only for the initial population structure. 
# --> Population structure changes because the initial population does not 
#     correspond to the stable age structure. 
# --> This means that if we keep "adjusting" f1yr_ad until reaching stable
#     age strucutre, we should get the correct result.

N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN

for(t in 2:length(N)){
  
  # Calculate "per adult recruitment rate" for time-step
  current_adultProp <- N_mat2[2, t-1] / sum(N_mat2[, t-1])
  f1yr_ad <- f1yr /  current_adultProp
  A[1, 2] <- f1yr_ad
  
  # Project
  N_mat2[, t] <- A %*% N_mat2[, t-1]
}

N_mat2
cbind(colSums(N_mat2), N_alt) 
# --> There we go, same as the non-structured model now

# Further thoughts: 
# -> We should probably compare PVAs done with i) two-year unstructured model, 
#    ii) one-year unstructured model, and iii) one-year age-structured model
# -> For iii), it would be best to start PVAs with a relatively "realistic"
#    population age structure. Can we look at the camera trap data and calculate
#    a proportion of juveniles (<1 year old) vs. adults?
# -> Regarding partitioning of reproduction: This is a bit trickier to do than
#    I first thought, mainly because of the way we have to infer annual recruitment. 
#    There may be way to "approximate" it in a Nimble model, but I am not sure
#    whether that makes sense. Something like this (code won't run, it's just a
#    doodle):


# Population projection
for(t in 2:Tmax){
  
  # New juveniles entering the population
  N[1, t] ~ dpois(N[2, t-1]*RecRate[t])
  
  # Juveniles and adults surviving (and remaining in population)
  N[2, t] ~ dbin(S1yr - E1yr, sum(N[1:2, t-1]))
}

# Informative priors for vital rates
S1yr <- S / sqrt(lambda)
f1yr <- f / sqrt(lambda)
E1yr <- E / sqrt(lambda)

E <- abs(lambda - S - f)

S ~ dnorm(mean = 0.49, sd = 0.091)
f ~ dnorm(mean = 0.76, sd = 0.12)
lambda ~ dnorm(mean = 1.16, sd = 0.12)

# Data likelihood for litter size data
kitten_count[x] ~ dpois(LS)
LS ~ dunif(0, 8)

# Data likelihood for breeding probability data
female_rep[x] ~ dbern(pB)
pB ~ dunif(0, 1)

# Link recruitment rate to vital rates / f1yr
for(t in 1:Tmax){
  log(RecRate[t]) ~ dnorm(mean = log(pB * LS * 0.5), sd = sigma_Res)
  
  RecRate[t] <- f1yr / (N[2, t] / sum(N[1:2, t]))
}
sigma_Res ~ dunif(0, 5)
# The ISSUE: This won't work, because we cannot have the same variable on the 
# left-hand side twice (multiple node definitions are not allowed).

# One way for approaching with might be to use another version of the dataset. 
# Let's call it "kitten_per_adult". It would be kitten counts on all pictures/
# individuals irrespective of whether i) adult is male/female and ii) any kittens 
# are present. The expectation there would be pB * LS, but we will link it to f1yr.

for(t in 1:Tmax){
  log(RecRate[t]) ~ dnorm(mean = log(0.5 * pB * LS), sd = sigma_Res)
  
  f1yr_ad[t] <- f1yr / (N[2, t] / sum(N[1:2, t]))
}
sigma_Res ~ dunif(0, 5)

kitten_per_adult[x] ~ dpois(mean(f1yr_ad[1:Tmax]))

# This solution is not ideal either because we have to depart from the original
# structure of the population model, and the link between f1yr and pB*LS is only
# implicit. 

# Another way might be to play around with the order of the implied relationship
# f1yr_ad[t] ~ 0.5 * pB * LS and use that to replace the non-informative priors
# for LS or pB. 

pB <- (2 * mean(f1yr_ad[1:Tmax]))/LS

# This will only work if pB does not appear in the population model though, i.e.
# if we define RecRate[t] <- f1yr_ad[t].
# This is again suboptimal because it restricts the information flow between 
# population dynamics and pB / LS. 

# All in all, none of the solutions for decomposing recruitment rate into different
# fecundity parameters post-hoc seems particularly promising. 
# I think this will have to be done in a fully integrated Bayesian model later
# (i.e. a model jointly analysing SCR and reproduction data).


