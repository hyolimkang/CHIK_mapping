library(randomForest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(exactextractr)
library(raster)
library(snowfall)
library(viridisLite)
library(gbm)
library(spatialRF)
#library(MODIStsp)
library(sf)
library(geodata)
library(ecospat)
library(NHANES)
library(dplyr)
library(class)
library(vctrs)
library(caret)
library(blockCV)
library(mgcv)
setwd("D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping") # window
source("Raster/fixNAs_adj.R")
source("plotRaster.R")
options(scipen = 999)
chik_foi<- read.csv("chik_foi.csv")

#1. load covariates data -------------------------------------------------------
load("Raster/covariates.RData")

#2. RF model (basic: train 100%) -----------------------------------------------
p_covs$FOI <- as.numeric(as.character(p_covs$FOI))
sapply(p_covs, class)
rf.mod <- randomForest(FOI ~ Temp + PRCP + Pop_dens + GDP + Albo + Aegyp + NDVI +Elev,
                       data = p_covs,
                       importance = TRUE)

pred.data <- data.frame(Temp        = as.vector(temp),
                        PRCP        = as.vector(precip), 
                        Pop_dens    = as.vector(pop_dens),
                        GDP         = as.vector(gdp_nat),
                        Albo        = as.vector(albo),
                        Aegyp       = as.vector(aegyp),
                        NDVI        = as.vector(ndvi),
                        Elev        = as.vector(elev))

rf.preds <- predict(rf.mod, pred.data, type = "response")

# transform all negative values to zero
rf.preds <- pmax(predict(rf.mod, pred.data, type = "response"), 0)


#3. GAM model ------------------------------------------------------------------
gam.mod <- gam(FOI ~ s(Temp) + s(PRCP) + s(Elev) + s(Pop_dens) + s(GDP) + s(Albo) + s(Aegyp) + s(NDVI),
               data = p_covs,
               method = "REML")
pred.data <- data.frame(Temp        = as.vector(temp),
                        PRCP        = as.vector(precip),
                        Elev        = as.vector(elev),
                        Pop_dens    = as.vector(pop_dens),
                        GDP         = as.vector(gdp_nat),
                        Albo        = as.vector(albo),
                        Aegyp       = as.vector(aegyp),
                        NDVI        = as.vector(ndvi))

plot(gam.mod, pages = 1, all.terms = T)

#4. brt model ------------------------------------------------------------------
brt.mod <- gbm(FOI ~ Temp + PRCP + Pop_dens + GDP + Albo + Aegyp + NDVI,
               data = p_covs,
               distribution = "gaussian", # FOI is continuous
               n.trees = 5000,            # Number of boosting trees (can be tuned)
               interaction.depth = 4,     # Depth of tree (can be tuned)
               shrinkage = 0.01,          # Learning rate (can be tuned)
               bag.fraction = 0.5)  

brt.preds <- predict(brt.mod, pred.data) 


#4. variable importance --------------------------------------------------------
# important variables
png(filename = paste0("varImp",
                      gsub("-", "_", Sys.Date()),
                      ".png"),
    width = 7, height = 5, units = "in", res = 300)
varImpPlot(rf.mod)
dev.off()
imp <-importance(rf.mod)
impvar <-rownames(imp)[order(imp[,1],decreasing=TRUE)]

# collect plots
plotList <- list()
for(i in seq_along(impvar)){ 
  plotList[[i]] <- partialPlot(rf.mod, p_covs, impvar[i],
                               xlab=impvar[i], 
                               main=paste("Partial Dependence on", impvar[i]),
                               plot = FALSE)
}
# plotting partial dependencies of surveillance model
png(filename = paste0("Chikmod_varPartial",
                      gsub("-", "_", Sys.Date()),
                      ".png"),
    width = 14, height = 10, units = "in", res = 300)
par(mfrow=c(4,3))
for(i in 1:length(plotList)){
  plot(plotList[[i]], type = "l",
       xlab = impvar[i],
       ylab = "Relative contribution to FOI")
  abline(h = 0, lty = 2)
}
dev.off()


#4.make map --------------------------------------------------------------------

template = temp

mapRasCon <- function(pred_vec) { 
  output <- template
  values(output) = as.numeric(pred_vec)
  return(output)
}

chikmap.rf <- mapRasCon(rf.preds) 
writeRaster(chikmap.rf, filename=paste0("Map_figs/raster", gsub("-", "_", Sys.Date()),".tif"), overwrite=TRUE)
gc()
rfmodel_plot <- plotRaster(chikmap.rf)
ggsave(filename = paste0("Map_figs/chik_map_rf",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), rfmodel_plot, height = 6, width = 12, dpi=900)


chikmap.gam <- mapRasCon(gam.preds) 
writeRaster(chikmap.gam, filename=paste0("Map_figs/raster", gsub("-", "_", Sys.Date()),".tif"), overwrite=TRUE)
gc()
gammodel_plot <- plotRaster(chikmap.gam)
ggsave(filename = paste0("Map_figs/chik_map_gam",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), gammodel_plot, height = 6, width = 12, dpi=900)

chikmap.brt <- mapRasCon(brt.preds) 
writeRaster(chikmap.brt, filename=paste0("Map_figs/raster", gsub("-", "_", Sys.Date()),".tif"), overwrite=TRUE)
gc()
brtmodel_plot <- plotRaster(chikmap.brt)
ggsave(filename = paste0("Map_figs/chik_map_brt",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), brtmodel_plot, height = 6, width = 12, dpi=900)


#5. masking map ----------------------------------------------------------------
# convert temperature suitablity to 0-1 index
# following: https://www.nature.com/articles/s41564-019-0476-8 
temp = temp / max(as.vector(temp), na.rm = T)


# temperature suitaiblity masking - 0.04 = 2 weeks of suitability = 1 serial interval
mask_ts <- (temp > 0.04)
# aedes aegypti masking
# from Kraemer paper: https://www.nature.com/articles/s41564-019-0376-y
# Aegypti thresholds = 0.13, 0.41, 0.70
# albopictus thresholds = 0.16, 0.36, 0.53
mask_ae = aegyp > 0.41
mask_al = albo > 0.36
mask_all = (mask_ts * (mask_ae + mask_al)) > 0 # areas with ae or al and suitable for dengue

CHIK_range_mask <- chikmap.rf * mask_all

r_CHIK <- plotRaster(CHIK_range_mask)+scale_fill_viridis(option = "rocket", direction = -1)
ggsave(filename = paste0("Map_figs/CHIK_rangemap",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), r_CHIK, height=6, width=12, dpi=900)

#check
raster_spdf <- as(chikmap.rf, "SpatialPixelsDataFrame")
raster_df <- as.data.frame(raster_spdf)


#6. cross validation------------------------------------------------------------
#https://github.com/HannaMeyer/CAST?tab=readme-ov-file

# first test-train split the data using createDataPartition. 
# Here we are using 75% of the data for training.
index = createDataPartition(p_covs$FOI, p = 0.75, list = FALSE)
# create test data set
foi_train = p_covs[index, ] # used for training and cross-validation
# create train data set
foi_test = p_covs[-index,] # held out and used later for final model evaluation

# To illustrate tuning, we perform k-nearest neighbors with 5-fold cross-validation.
foi_rf_mod = train(
  FOI ~ Temp + PRCP + Elev + Pop_dens + GDP + Albo + Aegyp + NDVI,
  data = foi_train,
  preProcess = c("center", "scale"),
  method = "rf",
  trControl = trainControl(method = "cv", number = 10,
                           savePredictions = "final")
)
foi_rf_mod$results
rf.cv.result <- foi_rf_mod$resample
colMeans(rf.cv.result[,1:3])

foi_glm_mod = train(
  FOI ~ Temp + PRCP + Elev + Pop_dens + GDP + Albo + Aegyp + NDVI,
  data = foi_train,
  preProcess = c("center", "scale"),
  method = "glm",
  trControl = trainControl(method = "cv", number = 10)
)

results <- resamples(list(RandomForest = foi_rf_mod, GLM = foi_glm_mod))
summary(results)


# Step 1: Train the final model on the entire training dataset (foi_train)
final_rf_mod <- randomForest(FOI ~ Temp + PRCP + Elev + Pop_dens + GDP + Albo + Aegyp + NDVI, data = foi_train)
# Step 2: Make predictions on the test dataset (foi_test)
foi_test_predictions <- predict(final_rf_mod, foi_test)
actual_values <- foi_test$FOI  # Assuming the response variable is named 'FOI'
rmse <- sqrt(mean((actual_values - foi_test_predictions)^2))
mae <- mean(abs(actual_values - foi_test_predictions))

#7. k-fold iterations-----------------------------------------------------------
predictions <- list()
test_data <- list()
train_data <- list()

folds <- createFolds(foi_train$FOI, k = 10, list = TRUE)
for(i in seq_along(folds)){
  train_fold <- foi_train[-folds[[i]],] # excluding i-th fold 
  test_fold  <- foi_train[folds[[i]],] # i-th fold itself 
  
  # Each train_fold and test_fold pair is different in each iteration of the loop 
  # because folds[[i]] changes, ensuring that each fold is used exactly once for validation.
  
  
  # train the model
  rf_mod <- randomForest(FOI ~Temp + PRCP + Elev + Pop_dens + GDP + Albo + Aegyp + NDVI, data = train_fold)
  
  # make predictions
  predictions[[i]] <- predict(rf_mod, test_fold)
  
  #store test data
  test_data[[i]] <- test_fold
  
  #store train data
  train_data[[i]] <- train_fold
}

# Ensure folds is created from foi_train, not from the entire p_covs.
# Each train_fold is a combination of 9 folds, and each test_fold is the remaining 1 fold (assuming 10-fold CV).
# After training and validating in the loop, you can average the predictions or performance metrics across all folds to estimate model performance.
# The final step would be to evaluate your model on foi_test to see how it performs on unseen data.

# rf.mod for each test data (10 times)

rf.preds.cv <- list()

for(i in 1:10){
  result <- predict(rf.mod, test_data[[i]], type = "response")
  rf.preds.cv[[i]] <- result
}


chikmap.rf.cv <- mapRasCon(rf.preds.cv[[1]]) 
writeRaster(chikmap.rf.cv, filename=paste0("Map_figs/raster", gsub("-", "_", Sys.Date()),".tif"), overwrite=TRUE)
gc()
rfmodel_plot <- plotRaster(chikmap.rf)
ggsave(filename = paste0("Map_figs/chik_map_rf",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), rfmodel_plot, height = 6, width = 12, dpi=900)

#8. evaluate cross validation predictions --------------------------------------

refmodel <- randomForest(FOI ~ Temp + PRCP + Elev + Pop_dens + GDP + Albo + Aegyp + NDVI, 
                         data = foi_train,
                         importance = TRUE)
test.predictions <- predict(refmodel, foi_test)


cv_predictions <- foi_rf_mod$pred
spatial_index <- cv_predictions$rowIndex
prediction_df <- data.frame(spatial_index, cv_predictions$pred)

average_predictions <- prediction_df %>%
  group_by(spatial_index) %>%
  summarize(avg_prediction = mean(cv_predictions.pred, na.rm = TRUE))


rf.preds.cv <- predict(final_rf_mod, pred.data, type = "response")
chikmap.rf.cv <- mapRasCon(rf.preds.cv) 
writeRaster(chikmap.rf.cv, filename=paste0("Map_figs/raster", gsub("-", "_", Sys.Date()),".tif"), overwrite=TRUE)
gc()
rfmodel_plot_cv <- plotRaster(chikmap.rf.cv)
ggsave(filename = paste0("Map_figs/chik_map_rf_cv",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), rfmodel_plot_cv, height = 6, width = 12, dpi=900)

#9. block CV -------------------------------------------------------------------
chik_mask_df <- as(CHIK_range_mask, "SpatialPixelsDataFrame")
chik_mask_sf <- st_as_sf(chik_mask_df)
sb1 <- cv_spatial(x = chik_mask_sf,
                  column = "FOI", # the response column
                  k = 5, # number of folds
                  size = 350000, # size of the blocks in metres
                  selection = "random", # random blocks to fold
                  iteration = 50, # find evenly dispersed folds
                  biomod2 = T)
  
  
