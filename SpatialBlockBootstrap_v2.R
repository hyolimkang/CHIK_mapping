
FOI <- list()
FOI_sf_list <- list()
grid_blocks_list <- list()
block_data_list <- list()
split_list <- list()
train_list <- list()
test_list <- list()
rf_model <- list()
rf_model_logit <- list()
merge_list <- list()
ensemble_list <- list()
ensemble_org <- list()
prediction <- list()
prediction_logit <- list()
prediction_logit_org <- list()

load("mcmc_tot.RData")
# 1. create 100 random FOI dataset and convert to sf ----------------------------------------
for(i in 1:100){
  
  current_dataset <- data.frame(FOI = numeric(76), lat = numeric(76), long = numeric(76))
  
  for(k in 1:nrow(chik_foi)){
    non_zero_values <- mcmc_tot[[k]][, 2][mcmc_tot[[k]][, 2] != 0]
    # randomly sample logit scale FOI from each distribution
    current_dataset$FOI[k] <- sample(non_zero_values, 1)
    current_dataset$logitFOI[k] <- log(current_dataset$FOI[k] / (1 - current_dataset$FOI[k]))
    # Assign corresponding latitude and longitude
    current_dataset$lat[k] <- mcmc_tot[[k]]$lat[1]
    current_dataset$long[k] <- mcmc_tot[[k]]$long[1]
    current_dataset$study_no[k] <- mcmc_tot[[k]]$study_no[1]
    current_dataset$ID[k]    <- mcmc_tot[[k]]$ID[1]
  }
  
  FOI[[i]] <- current_dataset
  
  # create blocks 
  FOI_sf <- st_as_sf(FOI[[i]], coords = c("long", "lat"))
  data_extent <- st_bbox(FOI_sf)
  FOI_sf_list[[i]] <- FOI_sf
}

FOI_df <- do.call(rbind, FOI_sf_list)
  
#2. Merge each of 100 FOI dataset with covariate dataset for a complete data ---------------

for(k in 1:length(FOI_sf_list)) {
  # Convert sf objects to regular data frames
  df1 <- as.data.frame(FOI_sf_list[[k]])
  
  # Merge the data frames based on "ID"
  merged_df <- left_join(df1, p_covs, by = "ID")
  
  merge_sf <- st_as_sf(merged_df, coords = c("Longitude", "Latitude"), crs = st_crs(FOI_sf_list[[k]]))
  
  merge_list[[k]] <- merge_sf
}

#3. Create grid blocks for each 100 FOI+covariate dataset (500x500km: 50 blocks each) -------------

for(i in 1:100){
  cell_size_deg <- 500 / 111
  
  # Create a grid of specified block size
  grid_blocks <- st_make_grid(merge_list[[i]], cellsize = c(cell_size_deg, cell_size_deg))
  grid_blocks_sf <- st_as_sf(grid_blocks)
  # Assign a unique block ID to each block
  grid_blocks_sf$blockID <- seq_len(nrow(grid_blocks_sf))
  grid_blocks_list[[i]] <- grid_blocks_sf
  
  # Extract longitude and latitude from geometry
  merge_list[[i]]$long <- st_coordinates(merge_list[[i]])[, 1]
  merge_list[[i]]$lat <- st_coordinates(merge_list[[i]])[, 2]
  
  # Perform a spatial join to identify which blocks contain data points
  blocks_with_data <- st_join(grid_blocks_list[[i]], merge_list[[i]], join = st_intersects)
  blocks_with_data <- blocks_with_data[!is.na(blocks_with_data$logitFOI.x), ]
  block_data_list[[i]] <- blocks_with_data
}

#4. For each 100 FOIdataset, randomly sample (N-1) train blocks (where N = total N blocks) and leave one block out for test -------- 

for(i in 1:100){
  split <- split(block_data_list[[i]], block_data_list[[i]]$blockID)
  split_list[[i]] <- split
  sampled_indices <- sample(1:length(split_list[[i]]), length(split_list[[i]]) - 1, replace = FALSE)
  non_samp <- setdiff(1:length(split_list[[i]]), sampled_indices)
  train_set <- do.call(rbind, split_list[[i]][sampled_indices])
  train_list[[i]] <- train_set ## seems like almost 99% of data are being used for training
  test_set <- split_list[[i]][non_samp]
  test_list[[i]] <- do.call(rbind, test_set)
}

# Parallel setup ---------------------------------------------------------------
# run each time
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

#5. Run RF model using logit FOI for each 100 FOI dataset-----------------------

start_time <- Sys.time()
rf_mod_logit <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  randomForest(logitFOI.x ~ Tsuit + PRCP + CHIKRisk + Albo + Aegyp + GDP + NDVI,
                                 data = train_list[[i]],
                                 importance = TRUE)
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs

# 6. Predict FOI for 1 held-out block in each of 100 FOI dataset ---------------

start_time <- Sys.time()
prediction_logit_org <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  logistic(predict(rf_mod_logit[[i]], newdata  = test_list[[i]], type = "response"))
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins


start_time <- Sys.time()
prediction_logit <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  predict(rf_mod_logit[[i]], newdata  = test_list[[i]], type = "response")
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

# 7. Make predictions using RF models trained for each of 100 FOI dataset-------

start_time <- Sys.time()
ensemble_list <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  logistic(predict(rf_mod_logit[[i]], newdata  = pred.data, type = "response"))
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 49.66074 mins

# 8. prediction dataframe (logit to original transform) ---------------------------

## add logit to original FOI to the test_list
test_list <- lapply(test_list, function(df) {
  adjusted_FOI <- logistic(df$logitFOI.x)
  df$FOI_org <- adjusted_FOI
  return(df)
})


combined_results_rf <- list()

for (i in 1:length(prediction_logit_org)) {
  # Create a data frame from the predictions
  pred_df <- data.frame(ID = test_list[[i]]$ID, Predicted_FOI = prediction_logit_org[[i]], Predicted_logit = prediction_logit[[i]])
  
  test_list[[i]]$test_ID <- i 
  
  test_df <- test_list[[i]][, c("ID", "FOI.x", "logitFOI.x", "x", "blockID", "region", "test_ID")]
  
  combined_df <- merge(test_df, pred_df, by = "ID")
  
  combined_df$RMSE <- sqrt(mean((combined_df$FOI.x - combined_df$Predicted_FOI)^2))
  
  combined_df$RMSE_logit <- sqrt(mean((combined_df$logitFOI.x - combined_df$Predicted_logit)^2))
  
  # Combine with the rest of the results
  combined_results_rf[[i]] <- combined_df
}

combined_df_rf <- do.call(rbind, combined_results_rf)

rmse_summary_rf <- combined_df_rf %>% group_by(blockID) %>% 
  summarise(
    meanRMSE = mean(RMSE, na.rm = TRUE),
    meanRMSE_logit =  mean(RMSE_logit, na.rm = TRUE),
    .groups = 'drop'
  )

mean_rmse <- mean(rmse_summary_rf$meanRMSE, na.rm = TRUE)

# 9. plots ---------------------------------------------------------------------

ggplot(data = rmse_summary_rf)+
  geom_sf(aes(fill = meanRMSE)) +
  scale_fill_viridis_c()+
  geom_sf(data = p_covs_sf, color = 'red', size = 2)+
  geom_sf_text(data = blocks_with_data, aes(label = blockID), size = 3, check_overlap = TRUE)+
  labs(title = "Spatial Blocks with Mean RMSE; RF", fill = "Mean RMSE") +
  theme_minimal() +
  theme(legend.position = "right")

x_lim <- c(0, 0.15)
y_lim <- c(0, 0.15)

# logitFOI -> original transformed
ggplot(combined_df_rf, aes(x = FOI.x, y = Predicted_FOI, color = factor(test_ID))) +
  geom_point() +  # Add points
  geom_point(alpha = 0.5) +
  labs(x = "Actual", y = "Fitted", title = "FOI, RF") +  # Labels and title
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  coord_fixed(ratio = 1) +  # Equal aspect ratio
  xlim(y_lim) +
  ylim(y_lim) +
  theme_bw()

# logitFOI 
x_lim <- c(-14, -1)
y_lim <- c(-14, -1)

ggplot(combined_df_rf, aes(x = logitFOI.x, y = Predicted_logit, color = factor(test_ID))) +
  geom_point() +  # Add points
  geom_point(alpha = 0.5) +
  labs(x = "Actual", y = "Fitted", title = "FOI, RF (logit)") +  # Labels and title
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  coord_fixed(ratio = 1) +  # Equal aspect ratio
  xlim(y_lim) +
  ylim(y_lim) +
  #geom_text(aes(label = ID), vjust = -0.5, hjust = 0.5) +
  theme_bw()

plot_list <- list()

for (i in 1:length(combined_results_rf)) {
  # Create each ggplot object and store in the list
  p <- ggplot(combined_results_rf[[i]], aes(x = FOI.x, y = Predicted_FOI)) +
    geom_point() +
    labs(x = "Actual", y = "Fitted", title = paste("FOI, RF", i)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_fixed(ratio = 1) +
    xlim(y_lim) +
    ylim(y_lim) +
    theme_bw()
  plot_list[[i]] <- p
}

pdf("combined_plots.pdf", width = 8, height = 8)

# Combine and print or save the plots in batches of 9
for (j in seq(1, length(plot_list), by = 9)) {
  # Combine 9 plots or the remainder if less than 9
  combined_plot <- wrap_plots(plot_list[j:min(j+8, length(plot_list))], ncol = 3)
  
  # If you want to display the plot in an R environment
  print(combined_plot)
}

dev.off()


# each of 100 logit FOI model fitted to the entire data (logit scale)
plot_list <- list()

pdf("combined_rf_alldata.pdf", width = 8, height = 8)
par(mar = c(2, 2, 2, 2))

for(j in seq(1, length(rf_mod_logit), by = 10)) {
  
  # Set up a multi-plot layout: 3 rows, 3 column
  par(mfrow = c(3, 3))
  
  # Loop through each chunk of 10 models
  for(i in j:min(j+9, length(rf_mod_logit))) {
    
    # Predict and plot
    p.rf <- logistic(predict(rf_mod_logit[[i]], newdata=p_covs, type = "response"))
    plot(p.rf ~ p_covs$FOI, asp=1, pch=20, xlab="actual", ylab="fitted", 
         main=paste("FOI, Random Forest", i),
         xaxs="i", yaxs="i", xlim=c(0,0.08),  ylim=c(0,0.08))
    grid()
    abline(h=0, v=0, col="gray")
    abline(a=0, b=1)
    
    text(p.rf, p_covs$logitFOI, labels=p_covs$ID, cex=0.7, pos=4)
  }
}


pdf("combined_rf_alldata_logit.pdf", width = 8, height = 8)
par(mar = c(2, 2, 2, 2))

for(j in seq(1, length(rf_mod_logit), by = 10)) {
  
  # Set up a multi-plot layout: 3 rows, 3 column
  par(mfrow = c(3, 3))
  
  # Loop through each chunk of 10 models
  for(i in j:min(j+9, length(rf_mod_logit))) {
    
    # Predict and plot
    p.rf <- predict(rf_mod_logit[[i]], newdata=p_covs, type = "response")
    plot(p.rf ~ p_covs$logitFOI, asp=1, pch=20, xlab="actual", ylab="fitted", 
         main=paste("FOI, Random Forest", i),
         xaxs="i", yaxs="i")
    grid()
    abline(h=0, v=0, col="gray")
    abline(a=0, b=1)
    
    text(p.rf, p_covs$logitFOI, labels=p_covs$ID, cex=0.7, pos=4)
  }
}




dev.off()

# 10. Ensemble map -------------------------------------------------------------
template = tsuit

mapRasCon <- function(pred_vec) { 
  output <- template
  values(output) = as.numeric(pred_vec)
  return(output)
}

ensemble_raster <- list()
for(i in 1:length(ensemble_list)){
  mapresult          <- mapRasCon(ensemble_list[[i]])
  ensemble_raster[[i]] <- mapresult 
}

bootstrap_stack <- stack(ensemble_raster)
avg_boogstrap_map <- calc(bootstrap_stack, fun = median)

boot_CHIK_rf <- plotRaster(avg_boogstrap_map)+scale_fill_viridis(option = "rocket", direction = -1)
ggsave(filename = paste0("Map_figs/CHIK_boot_rf",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), boot_CHIK_rf, height=6, width=12, dpi=900)



