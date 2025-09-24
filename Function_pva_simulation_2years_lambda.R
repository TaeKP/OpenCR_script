
# Function to simulate growth_rate with stochasticity
# 2-YEAR POPULATION MODEL (Population growth rate only)

pva_simulation <- function(initial_population, growth_rate, carrying_capacity, n_years) {
  population <- numeric(n_years)
  #population[1] <- truncnorm::rtruncnorm(1, mean = initial_population, sd = initialN_sd, a = 0) # we probably want this lognormal (or truncated at 0)
  population[1] <- initial_population
  
  for (year in 2:n_years) {
    
    stochastic_growthrate <- truncnorm::rtruncnorm(1, mean = growth_rate, sd = growth_rate_sd, a = 0) # Adding randomness
    
    population[year] <- population[year - 1] * stochastic_growthrate
    
    #population[year] <- ifelse(population[year] > carrying_capacity, carrying_capacity, population[year])
    population[year] <- ifelse(population[year] > carrying_capacity, population[year-1], population[year])
    
    if (population[year] < 1) { # Extinction event
      population[year] <- 0
      break
    }
  }
  return(population)
}