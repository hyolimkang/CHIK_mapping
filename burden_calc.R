
load("rasterDf.RData")

pop_count <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Worldpop_population_2020_5k.tif"
pop_count <- raster(pop_count)
plot(pop_count)

pop_count  <- projectRaster(pop_count, elev, method="bilinear")
plot(pop_count)

pop_spdf <- as(pop_count, "SpatialPixelsDataFrame")
pop_df <- as.data.frame(pop_count, xy = TRUE, na.rm = T)


foi_comb <- merge(raster_df, pop_df, by = c("x", "y"))
foi_comb$tot_pop <- foi_comb$Worldpop_population_2020_5k * 5000

foi_comb_rast <- foi_comb[,c(1,2,5)]


# incidence calculations https://github.com/lorecatta/drep/blob/master/R/calculating_R0.R

calc_incidence <- function(FOI, l_lim, u_lim){
  if (FOI <=0) {
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
