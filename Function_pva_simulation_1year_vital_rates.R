
# Function to simulate growth_rate with stochasticity
# 1-YEAR POPULATION MODEL (vital rates)

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