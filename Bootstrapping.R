#https://github.com/giscience-fsu/sperrorest

library(sperrorest)
library(ranger)

num_samples <- 200
bootstrap_samp <- vector("list", num_samples)

for(i in 1:num_samples){
  sampled_indices <- sample(nrow(p_covs), nrow(p_covs), replace = T)
  bootstrap_samp[[i]] <- p_covs[sampled_indices, ]
}

# for each bootstrapped sample, create train-validation data (store in a list)

index <- list()
train <- list()
test  <- list()

for(i in 1:length(bootstrap_samp)){
  result = createDataPartition(bootstrap_samp[[i]]$FOI, p = 0.75, list = FALSE)
  index[[i]] <- result
  train[[i]] <- bootstrap_samp[[i]][index[[i]], ]
  test[[i]] <- bootstrap_samp[[i]][-index[[i]], ]
}

# fit the RF model to training dataset for each iteration of bootstrap

fit.train.model <- list()

for(i in 1:length(bootstrap_samp)){
  result <- ranger(FOI ~ Temp + PRCP + Elev + Pop_dens + GDP + Albo + Aegyp + NDVI,
                            data = train[[i]])
  fit.train.model[[i]] <- result
}

# use the model to predict the observation in the validation dataset

obs.validation <- list()

for(i in 1:length(bootstrap_samp)){
  result <- predict(fit.train.model[[i]],
                    data = test[[i]])
  
  predictions <- result$predictions
    
  df <- data.frame(dep.var     = test[[i]]$FOI,
                   predictions = predictions)
  
  obs.validation[[i]] <- df
}

# calculate sum of sqaured erros (sse)

sse <- list()

for(i in 1:length(bootstrap_samp)){
  
  result <- sum((obs.validation[[i]]$dep.var - obs.validation[[i]]$predictions)^2)
  sse[[i]] <- result
}

sse_val <- unlist(sse)
avg.sse <- mean(sse_val)
ssd_df <- as.data.frame(sse_val)

# 95%CI for predicted FOI

bootstrap.percentile <- lapply(obs.validation, function(sample) {
  
  predictions <- sample$predictions
  
  data.frame(
    Lower = quantile(predictions, probs = 0.025),
    Upper = quantile(predictions, probs = 0.975)
  )
})

combine.percentile <- do.call(rbind, bootstrap.percentile)

overall_lo <- quantile(combine.percentile$Lower, probs = 0.025)
overall_hi <- quantile(combine.percentile$Upper, probs = 0.975)


# block bootstrap package 
represampling_tile_bootstrap(
  data,
  coords = c("x", "y"),
  repetition = 1,
  nboot = -1,
  seed1 = NULL,
  oob = FALSE,
  ...
)

boostrap <- represampling_tile_bootstrap(data = chik_foi_all,
                             coords = c("long", "lat"),
                             repetition = 10,
                             nsplit = c(5, 5),
                             nboot = 1,
                             seed1 = NULL,
                             oob = FALSE)

train_data <- boostrap$training
test_data <- boostrap$testing


