
# Function to simulate growth_rate with stochasticity
# 2-YEAR POPULATION MODEL (vital rates)

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