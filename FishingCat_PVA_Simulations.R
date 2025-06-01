
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

# Function to simulate growth_rate with stochasticity

source("Function_pva_simulation_2years_vital_rates.R")
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

source("Function_pva_simulation_1year_vital_rates.R")

#-------------------------------------------------------------------------------
# Running the simulation multiple times
results_1yr <- replicate(simulations, pva_simulation_1yr(initN, growth_rate, survival_rate, recruitment_rate, carrying_capacity, n_years))

# Calculate the probability of extinction by the end of the simulation period
(extinction_probability_1yr <- mean(results_1yr[n_years,] == 0))


# 1-YEAR POPULATION WITH 2 AGE CLASSES (Vital rates, age structure) #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------------------
## Function to simulate the stochasticity, equivalent to the function "pva simulation()".
# Baseline scenario

source("Function_pva_simulation_age_str.R")

#-------------------------------------------------------------------------------
# Running the simulation multiple times
(results_age_str <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                  growth_rate, growth_rate_sd,
                                                                  survival_rate, survival_rate_sd,
                                                                  recruitment_rate, recruitment_rate_sd,
                                                                  init_adultProp, init_adultProp_SD,
                                                                  carrying_capacity, 
                                                                  pertFac.S = 1, pertFac.f = 1, pertFac.E = 1,
                                                                  n_years)) )

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
## Function to simulate the stochasticity.
## Calling the baseline pva simulation function
source("Function_pva_simulation_age_str.R")

#-------------------------------------------------------------------------------
## Scenario 1	
## 10% Decreasing the survival rate 
## perturbation factor "S" at 0.9

# Running the simulation multiple times
(results_age_str_s1 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 0.9, pertFac.f = 1, pertFac.E = 1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s1 <- reshape2::melt(apply(results_age_str_s1, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "10% survival decreased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s1 <- rbind(sim_s1, sim_4) %>%
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
## 20% Decreasing the survival rate 
## perturbation factor "S" at 0.8

# Running the simulation multiple times
(results_age_str_s2 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 0.8, pertFac.f = 1, pertFac.E = 1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s2 <- reshape2::melt(apply(results_age_str_s2, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "20% survival decreased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s2 <- rbind(sim_s2, sim_4) %>%
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
## 30% Decreasing the survival rate 
## perturbation factor "S" at 0.7

# Running the simulation multiple times
(results_age_str_s3 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 0.7, pertFac.f = 1, pertFac.E = 1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s3 <- reshape2::melt(apply(results_age_str_s3, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "30% survival decreased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s3 <- rbind(sim_s3, sim_4) %>%
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

## Comparing 3 decreasing survival scenarios with baseline scenario
## Combine results and summarise
simSummary_s1t3 <- rbind(sim_s1, sim_s2, sim_s3, sim_4) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s1t3, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 4	
## 10% Decreasing the recruitment rate 
## perturbation factor "f" at 0.9

# Running the simulation multiple times
(results_age_str_s4 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 1, pertFac.f = 0.9, pertFac.E = 1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s4 <- reshape2::melt(apply(results_age_str_s4, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "10% recruitment decreased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s4 <- rbind(sim_s4, sim_4) %>%
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
## Scenario 5	
## 20% Decreasing the recruitment rate 
## perturbation factor "f" at 0.8

# Running the simulation multiple times
(results_age_str_s5 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 1, pertFac.f = 0.8, pertFac.E = 1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s5 <- reshape2::melt(apply(results_age_str_s5, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "20% recruitment decreased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s5 <- rbind(sim_s5, sim_4) %>%
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
## 30% Decreasing the recruitment rate 
## perturbation factor "f" at 0.7

# Running the simulation multiple times
(results_age_str_s6 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 1, pertFac.f = 0.7, pertFac.E = 1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s6 <- reshape2::melt(apply(results_age_str_s6, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "30% recruitment decreased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s6 <- rbind(sim_s6, sim_4) %>%
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

## Comparing 3 decreasing recruitment scenarios with baseline scenario
## Combine results and summarise
simSummary_s4t6 <- rbind(sim_s4, sim_s5, sim_s6, sim_4) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s4t6, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------
## Scenario 7
## Due to the high threat and poor management could cause animals emigrate or move out to some better area
## Therefore, the emigration rate must be increased
## 10% increasing the emigration rate 
## perturbation factor "E" at 1.1

# Running the simulation multiple times
(results_age_str_s7 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 1, pertFac.f = 1, pertFac.E = 1.1,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s7 <- reshape2::melt(apply(results_age_str_s7, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "10% emigration increased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s7 <- rbind(sim_s7, sim_4) %>%
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
## Scenario 8
## Due to the high threat and poor management could cause animals emigrate or move out to some better area
## Therefore, the emigration rate must be increased
## 20% increasing the emigration rate 
## perturbation factor "E" at 1.2

# Running the simulation multiple times
(results_age_str_s8 <- replicate(simulations, pva_simulation_age_str(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 1, pertFac.f = 1, pertFac.E = 1.2,
                                                                     n_years)) )

# Write results as data frames
# compare with baseline scenario

sim_s8 <- reshape2::melt(apply(results_age_str_s8, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "20% emigration increased")

#baseline scenario
sim_4 <- reshape2::melt(apply(results_age_str, c(2, 3), sum)) %>%
  dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
  dplyr::mutate(Model = "1-year, vital rates & age structure")

## Combine results and summarise
simSummary_s8 <- rbind(sim_s8, sim_4) %>%
  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
  dplyr::group_by(Model, Year) %>%
  dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

## Plot
ggplot(simSummary_s8, aes(x = Year, group = Model)) + 
  geom_line(aes(y = median_N, color = Model)) + 
  geom_ribbon(aes(ymin = lCI_N, ymax = uCI_N, fill = Model), alpha = 0.2) + 
  xlim(1, n_years-1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw()

#-------------------------------------------------------------------------------