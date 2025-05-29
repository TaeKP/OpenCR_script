
# pertFac.X = 1 --> vital rate remains unchanged
# pertFac.X < 1 (but >= 0) --> vital rate decreases proportionally
# pertFac.X > 1 --> Vital rate increases

pva_simulation_age_str <- function(initN, 
                                   growth_rate, growth_rate_sd,
                                   survival_rate, survival_rate_sd,
                                   recruitment_rate, recruitment_rate_sd,
                                   init_adultProp, init_adultProp_SD,
                                   carrying_capacity, 
                                   pertFac.S, pertFac.f, pertFac.E,
                                   n_years) {
  
  S <- truncnorm::rtruncnorm(1, mean = survival_rate, sd = survival_rate_sd, a = 0, b = 1) # True survival (2 yrs)
  f <- truncnorm::rtruncnorm(1, mean = recruitment_rate, sd = recruitment_rate_sd, a = 0) # Recruitment (2 yrs)
  
  lambda <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0, b = S + f) # Population growth rate
  
  E = abs(lambda - S - f) # Emigration rate
  
  init_adultProp <- truncnorm::rtruncnorm(1, mean = init_adultProp_mean, sd = init_adultProp_SD, a = 0, b = 1)
  
  N_mat2 <- matrix(NA, nrow = 2, ncol = n_years)
  N_mat2[,1] <- c(1 - init_adultProp, init_adultProp)*initN
  
  S1yr <- (S / sqrt(lambda))*pertFac.S 
  f1yr <- (f / sqrt(lambda))*pertFac.f
  E1yr <- (E / sqrt(lambda))*pertFac.E
  
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
