# Appendix D - Population viability analysis: 
# Last check: 29 September 2026

# required packages
 library(magrittr)
 library(ggplot2)
 library(truncnorm)    
 library(dplyr)        
 library(reshape2)    

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

carrying_capacity <- 116     		# based on the suitable habitat and FC's home range size in KSRY and habitat improvement ; assigned < 83, 93, 116 individuals; for varying CC
n_years <- 10                  		# Number of years to simulate
simulations <- 10000          		# Number of simulations

#-----------------------------------------------------------------------------------------------------
# Model_1: Incorporating vital rates model for 1 year and age‑structure (Baseline scenario)
#-----------------------------------------------------------------------------------------------------
# Function to simulate PVA: Growth rate + Survival + Recruitment + age‑structure + stochasticity, 
# equivalent to the function "pva_simulation_GSR1ys()" in Appendix C2, but Model_1 for annual vital rates and age‑structure

# Example for the factor perturbation of survival rate (S)
# pertFac.S = 1 --> vital rate remains unchanged
# pertFac.S < 1 (but >= 0) --> vital rate decreases proportionally
# pertFac.S > 1 --> Vital rate increases

pva_simulation_GSR1yr_age <- function(initN, 
                                   growth_rate, growth_rate_sd,
                                   survival_rate, survival_rate_sd,
                                   recruitment_rate, recruitment_rate_sd,
                                   init_adultProp, init_adultProp_SD,
                                   carrying_capacity, 
                                   pertFac.S = 1, pertFac.f = 1, pertFac.E = 1,
                                   n_years) {
            
            N <- c(initN, rep(NA, n_years - 1))  # add N at time t1
            
            # True survival (2 years interval)
            S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) 
  
            # Recruitment (2 years interval)
            f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) 
  
            # Population growth rate
            lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) 
  
            # Emigration rate
            E = abs(lambda - S - f) 
  
            # Initial adult Proportion
            init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
            N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
            N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
            S1yr <- (S / sqrt(lambda)) * pertFac.S 
            f1yr <- (f / sqrt(lambda)) * pertFac.f
            E1yr <- (E / sqrt(lambda)) * pertFac.E
  
            for(t in 2:length(N)){
    
              # Projection matrix (1 = juvenile, < 1 year old; 2 = adult, > 1 year old)
              A <- matrix(NA, nrow = 2, ncol = 2)
              A[1, 1] <- 0 				        # Juveniles producing juveniles within 1 year
              A[2, 1] <- S1yr - E1yr 			# Juvenile becoming adults within 1 year
              A[2, 2] <- S1yr - E1yr 			# Adults remaining alive (& in study area) within 1 year
              A[1, 2] <- f1yr 				    # Adults producing juveniles within 1 year
    
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

# Running the simulation
results_M1 <- replicate(simulations, pva_simulation_GSR1yr_age(initN, 
                                                                  growth_rate, growth_rate_sd,
                                                                  survival_rate, survival_rate_sd,
                                                                  recruitment_rate, recruitment_rate_sd,
                                                                  init_adultProp, init_adultProp_SD,
                                                                  carrying_capacity, 
                                                                  pertFac.S = 1, pertFac.f = 1, pertFac.E = 1,
                                                                  n_years)) 

# Calculate the probability of extinction by the end of the simulation period
# Extract population sizes at the final year across all samples
final_pop <- results_M1[, n_years, ]   # age classes × samples

# Sum across age classes to get total population per sample
final_total <- colSums(final_pop, na.rm = TRUE)

# Extinction indicator: 1 if total population == 0
extinct <- final_total == 0     # assigned < 50 ; for quasi-extinction 

# Probability of extinction
extinction_probability_age_str <- mean(extinct)
extinction_probability_age_str

#----------------------------------------- END of model 1 ------------------------------------------------------------
# Plot the baseline scenario, which obtain the data from sim_1

sim_1 <- reshape2::melt(apply(results_M1, c(2, 3), sum)) %>%
        dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
        dplyr::mutate(Model = "1-year, vital rates and age structure")

## Combine results and summarize from sim_1
simSummary_base <- rbind(sim_1) %>%
                  dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
                  dplyr::group_by(Model, Year) %>%
                  dplyr::summarise(mean_N = mean(PopSize),
                    median_N = median(PopSize),
                    sd_N = sd(PopSize),
                    lCI_N = quantile(PopSize, probs = 0.025),
                    uCI_N = quantile(PopSize, probs = 0.975),
                    .groups = "keep") 

# Black line color and Blue median line Plot
ggplot() +
        # spaghetti plot (all simulations in black)
        geom_line(data = sim_1,
            aes(x = Year, y = PopSize, group = SimNo),
            color = "black", alpha = 0.1) +
  
        # median line (blue)
        geom_line(data = simSummary_base,
            aes(x = Year, y = median_N),
            color = "blue", linewidth = 1) +
  
        # custom x-axis ticks      
        scale_x_continuous(breaks = seq(0, n_years, by = 2), expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +   
  
        # horizontal line for carrying capacity
        geom_hline(yintercept = carrying_capacity, linetype = "dashed", color = "red", linewidth = 1) +
  
        theme_minimal() +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(color = "black")) +
  
        labs(title = "a) Baseline",
            x = "Year", y = "Population size")

#----------------------------------------- END of baseline scenario plot ---------------------------------------------

#-----------------------------------------------------------------------------------------------------
# Model_2: Incorporating vital rates model for 1 year and age‑structure (Management scenario 1)
#-----------------------------------------------------------------------------------------------------
# Obtained the function to simulate PVA from model_1 "pva_simulation_GSR1yr_age", 
# but Model_2 perturbation factor "S" at 0.95, which is 5% survival rate decreased

# Running the simulation
results_M2 <- replicate(simulations, pva_simulation_GSR1yr_age(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 0.95, pertFac.f = 1, pertFac.E = 1,
                                                                     n_years)) 

# Calculate the probability of extinction by the end of the simulation period
# Extract population sizes at the final year across all samples
final_pop_s1 <- results_M2[, n_years, ]   

# Sum across age classes to get total population per sample
final_total_s1 <- colSums(final_pop_s1, na.rm = TRUE)

# Extinction indicator: 1 if total population == 0
extinct_s1 <- final_total_s1 == 0   # assigned < 50 ; for quasi-extinction 

# Probability of extinction
extinction_probability_age_str_s1 <- mean(extinct_s1)
extinction_probability_age_str_s1

#----------------------------------------- END of model 2 ------------------------------------------------------------

# Plot the management scenario 1
# Write results as data frames
sim_2 <- reshape2::melt(apply(results_M2, c(2, 3), sum)) %>%
          dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
          dplyr::mutate(Model = "5% survival decreased")

## Combine results and summarize from sim_2
simSummary_s1 <- rbind(sim_2) %>%
                dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
                dplyr::group_by(Model, Year) %>%
                dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

# Black line color and Blue median line plot
ggplot() +
        # spaghetti plot (all simulations in black)
        geom_line(data = sim_2,
            aes(x = Year, y = PopSize, group = SimNo),
            color = "black", alpha = 0.1) +
        
        # median line (blue)
        geom_line(data = simSummary_s1,
            aes(x = Year, y = median_N),
            color = "blue", linewidth = 1) +
        
        # custom x-axis ticks      
        scale_x_continuous(breaks = seq(0, n_years, by = 2), expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
  
        # horizontal line for carrying capacity
        geom_hline(yintercept = carrying_capacity, linetype = "dashed", color = "red", linewidth = 1) +
  
        theme_minimal() +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(color = "black")) +
  
        labs(title = "b) Survival rate decreased",
            x = "Year", y = "Population size")

#----------------------------------------- END of management scenario 1 plot -----------------------------------------

#-----------------------------------------------------------------------------------------------------
# Model_3: Incorporating vital rates model for 1 year and age‑structure (Management scenario 2)
#-----------------------------------------------------------------------------------------------------
# Obtained the function to simulate PVA from model_1 "pva_simulation_GSR1yr_age", 
# but Model_3 perturbation factor "f" at 0.95, which is 5% recruitment rate decreased

# Running the simulation
results_M3 <- replicate(simulations, pva_simulation_GSR1yr_age(initN, 
                                                                     growth_rate, growth_rate_sd,
                                                                     survival_rate, survival_rate_sd,
                                                                     recruitment_rate, recruitment_rate_sd,
                                                                     init_adultProp, init_adultProp_SD,
                                                                     carrying_capacity, 
                                                                     pertFac.S = 1, pertFac.f = 0.95, pertFac.E = 1,
                                                                     n_years)) 

# Calculate the probability of extinction by the end of the simulation period
# Extract population sizes at the final year across all samples
final_pop_s2 <- results_M3[, n_years, ]   

# Sum across age classes to get total population per sample
final_total_s2 <- colSums(final_pop_s2, na.rm = TRUE)

# Extinction indicator: 1 if total population == 0
extinct_s2 <- final_total_s2 == 0   # assigned < 50 ; for quasi-extinction 

# Probability of extinction
extinction_probability_age_str_s2 <- mean(extinct_s2)
extinction_probability_age_str_s2

#----------------------------------------- END of model 3 ------------------------------------------------------------

# Plot the management scenario 2
# Write results as data frames
sim_3 <- reshape2::melt(apply(results_M3, c(2, 3), sum)) %>%
          dplyr::rename(Year = Var1, SimNo = Var2, PopSize = value) %>%
          dplyr::mutate(Model = "5% recruitment decreased")

## Combine results and summarize from sim_3
simSummary_s2 <- rbind(sim_3) %>%
                dplyr::mutate(PopSize = ifelse(is.na(PopSize), 0, PopSize)) %>%
                dplyr::group_by(Model, Year) %>%
                dplyr::summarise(mean_N = mean(PopSize),
                   median_N = median(PopSize),
                   sd_N = sd(PopSize),
                   lCI_N = quantile(PopSize, probs = 0.025),
                   uCI_N = quantile(PopSize, probs = 0.975),
                   .groups = "keep") 

# Black line color and Blue median line plot
ggplot() +
        # spaghetti plot (all simulations in black)
        geom_line(data = sim_3,
            aes(x = Year, y = PopSize, group = SimNo),
            color = "black", alpha = 0.1) +
      
        # median line (blue)
        geom_line(data = simSummary_s2,
            aes(x = Year, y = median_N),
            color = "blue", linewidth = 1) +
        
        # custom x-axis ticks
        scale_x_continuous(breaks = seq(0, n_years, by = 2), expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +  
  
        # horizontal line for carrying capacity
        geom_hline(yintercept = carrying_capacity, linetype = "dashed", color = "red", linewidth = 1) +
  
        theme_minimal() +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(color = "black")) +
  
        labs(title = "c) Recruitment rate decreased",
            x = "Year", y = "Population size")

# Create an extinction risk table
extinction_set <- data.frame(
                  Scenario = c("Baseline", "5% survival decreased", "5% recruitment decreased"),
                  Extinction_probability = c(extinction_probability_age_str, 
                             extinction_probability_age_str_s1, 
                             extinction_probability_age_str_s2))
extinction_set
#----------------------------------------- END NOT RUN------------------------------------------------------------