# burden functions--------------------------------------------------------------

# 1. calculate incidence -------------------------------------------------------

calc_incidence <- function(FOI, l_lim, u_lim){
  if (FOI < 0) {
    stop("FOI must be greater than zero") # debugging
  }
  
  get_incidence <- function(infec_prob){
    (infec_prob) / (u_lim - l_lim)
  } # nested function: this function recognize infec_prob = p_I 
  
  p_I <- exp(-FOI * l_lim) - exp(-FOI * u_lim) # probability of infection over a specific range
  
  incidence <- get_incidence(p_I) # passing incidence into nested function
  # converts incidence per unit range 
  return(incidence)
}

#2. incidence to total infection numbers ---------------------------------------

incidence_to_numbs <- function(incidence, n_j){
  incidence * n_j # n_j is total population per pixel j 
}

calc_infections <- function(FOI,
                            n_j,
                            u_lim,
                            l_lim){
  incids <- calc_incidence(FOI, u_lim, l_lim)
  
  infection_numbers_j <- lapply(incids, incidence_to_numbs, n_j)
  
  tot_infection       <- vapply(infection_numbers_j, sum, numeric(1))
  
  sum(tot_infection)
}

# 3. estimate total infections and age-stratified burden -----------------------

# function with returning only original df
age_strat_burden_org <- function(age_vec_df, df, n_age_groups, l_lim, u_lim) {
  
  for(i in 1:nrow(df)){
    
    #1. fetch the age proportion 
    start_idx <- (i - 1) * n_age_groups + 1
    
    end_idx   <- i * n_age_groups
    
    age_vec_current <- age_vec_df[start_idx:end_idx, 'prop']
    
    #2. breakdown total population in each row 
    
    age_str_pop <- df$pop[i] * age_vec_current
    
    #3. incidence in each specific age band 
    
    incidence_rates <- mapply(calc_incidence, FOI = rep(df$FOI[i], n_age_groups), l_lim, u_lim)   # calculate incidence rates for each age band
    
    #4. infections per age band
    
    infection <- mapply(incidence_to_numbs, incidence = incidence_rates, n_j = age_str_pop)  # calculates the number of infections for each age band by multiplying the incidence rates by the corresponding populations in age_str_pop
    
    #5. total sum of infections per row 
    
    total_infection[i] <- sum(infection)  # sums up the infections for all age bands for the ith row
  }
  
  # Add it to the original dataframe
  
  df$total_infection <- total_infection
  return(df)
}

# function with returning both age-stratified burden and original df
# matching two dfs by country code in a wide format dataframe 
age_strat_burden <- function(age_vec_df, df) {
  
  # define age band
  
  n_age_groups <- 20
  
  l_lim <- seq(0, 95, length.out = n_age_groups)
  
  u_lim <- seq(5, 100, length.out = n_age_groups)
  
  # initialise a matrix to store infections for each age band for each row of df
  
  infections_per_age_band <- matrix(nrow = nrow(df), ncol = n_age_groups)
  
  total_infection <- numeric(nrow(df))
  
  for(i in 1:nrow(df)){
    
    #1. match the country 
    
    country_row <- match(df$iso3[i], age_vec_df$iso3)
    
    if(!is.na(country_row)) {
      age_vec_current <- age_vec_df[country_row, 4:23]
      
      #2. breakdown total population in each row 
      
      age_str_pop <- df$pop[i] * age_vec_current
      
      #3. incidence in each specific age band 
      
      incidence_rates <- mapply(calc_incidence, FOI = rep(df$layer[i], n_age_groups), l_lim, u_lim)   # calculate incidence rates for each age band
      
      #4. infections per age band
      
      infection <- mapply(incidence_to_numbs, incidence = incidence_rates, n_j = age_str_pop)  # calculates the number of infections for each age band by multiplying the incidence rates by the corresponding populations in age_str_pop
      
      #5. store the infections per age band in the matrix 
      
      infections_per_age_band[i, ] <- infection
      
      #6. total sum of infections per row 
      
      total_infection[i] <- sum(infection)  # sums up the infections for all age bands for the ith row
      
    } else {
      warning(paste("No matching country found in age_vec_df for"), df$iso3[i])
    }
  }
  
  # Add it to the original dataframe
  
  df$total_infection <- total_infection
  
  return(list(
    updated_df         = df,
    infection_per_band = infections_per_age_band
  ))
}
