#!/usr/bin/env Rscript
args_unparsed <- commandArgs(trailingOnly = TRUE)
args <- strsplit(args_unparsed, ",")
p <- as.integer(args[[1]])
spatial_basis <- as.numeric(args[[2]])
temporal_basis <- as.numeric(args[[3]])
ar_order <- as.integer(args[[4]])
aux_dim <- as.integer(args[[5]])
file_name <- paste0("pred_", args_unparsed[1], "_", args_unparsed[2], "_", args_unparsed[3], "_", args_unparsed[4], "_", args_unparsed[5], ".RData")

.libPaths(c("/projappl/project_2003994/project_rpackages", .libPaths()))

library(NonlinearBSS)
library(tensorflow)
library(keras3)

seed <- 29012025
temp_shift <- 80
n_s <- 58

data_full_train <- read.csv("/scratch/project_2003994/iVAEar1/athens_daily_train.csv")
head(data_full_train)
data_full_test <- read.csv("/scratch/project_2003994/iVAEar1/athens_daily_test.csv")
coords_time_train <- read.csv("/scratch/project_2003994/iVAEar1/athens_coords_train.csv")
head(coords_time_train)
coords_time_train <- coords_time_train[, 2:4]
coords_time_test <- read.csv("/scratch/project_2003994/iVAEar1/athens_coords_test.csv")
coords_time_test <- coords_time_test[, 2:4]
head(coords_time_train)

validation_time_points <- which(data_full_train$Date_numeric %in% 1091:1100)
data_full_validation <- data_full_train[validation_time_points, ]
coords_time_validation <- coords_time_train[validation_time_points, ]
data_train <- data_full_train[-validation_time_points, ]
coords_time_train <- coords_time_train[-validation_time_points, ]


variables <- c("Wind.Speed..U.", "Wind.Speed..V.", "Dewpoint.Temp", "Soil.Temp", "Temp", "Relative.Humidity", "PM10", "PM2.5", "NO2", "O3")
seas_var_names <- sapply(variables, function(name) paste0(name, "_seas"))
res_var_names <- sapply(variables, function(name) paste0(name, "_res"))

results <- data.frame(matrix(NA, nrow = 3, ncol = 6))
names(results) <- c("wMSE", "p", "spatial_basis", "temporal_basis", "ar_order", "method")

resiVAEar1 <- iVAEar_radial(as.matrix(data_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3] + temp_shift, p, n_s = n_s, 
    spatial_basis = spatial_basis, temporal_basis = temporal_basis, seasonal_period = 365, ar_order = ar_order, 
    aux_hidden_units = c(aux_dim), epochs = 40, initial_batch_size = 64, error_dist_sigma = 0.02, seed = seed,
    lr_start = 0.0001, lr_end = 0.0001)
resiVAEar1$week_component <- FALSE
unscaled <- FALSE
resiVAEar1$ar_order <- ar_order
st_trend_test_time <- predict_coords_to_IC(resiVAEar1, coords_time_validation[, 1:2], coords_time_validation[, 3] + temp_shift, unscaled = unscaled)
last_coords_time <- coords_time_train[which(coords_time_train[, 3] %in% ((1090 - (ar_order - 1)):1090)), ]
st_ar1_test_time <- predict_coords_to_IC_ar(resiVAEar1, last_coords_time[, 1:2], last_coords_time[, 3] + temp_shift, NULL, coords_time_validation[, 1:2], coords_time_validation[, 3] + temp_shift, NULL, unscaled = unscaled)
st_ar1_test_time <- st_ar1_test_time$preds

pred_obs_time_aux <- predict(resiVAEar1, st_ar1_test_time, IC_to_data = TRUE)
pred_obs_time_aux_mean <- predict(resiVAEar1, st_trend_test_time, IC_to_data = TRUE)

preds_mean <- pred_obs_time_aux_mean #+ data_full_test[, seas_var_names]
preds_aux <- pred_obs_time_aux #+ data_full_test[, seas_var_names]

mse_time_aux_mean <- apply((data_full_validation[, variables] - preds_mean)^2, 2, mean)
mse_time_aux <- apply((data_full_validation[, variables] - preds_aux)^2, 2, mean)
mse_time_aux_mean
mse_time_aux
sds <- apply(data_full_train[, res_var_names], 2, sd)
results[1, ] <- c(mean(mse_time_aux_mean/sds^2), p, toString(spatial_basis), toString(temporal_basis), ar_order, "iVAEar_trend")
results[2, ] <- c(mean(mse_time_aux/sds^2), p, toString(spatial_basis), toString(temporal_basis), ar_order, "iVAEar")

resiVAE_mean <- iVAE_radial_spatio_temporal(as.matrix(data_train[, variables]), as.matrix(coords_time_train[, 1:2]), 
    coords_time_train[, 3] + temp_shift, p, spatial_basis = spatial_basis, temporal_basis = temporal_basis, aux_hidden_units = c(aux_dim),
    seasonal_period = 365, epochs = 60, batch_size = 64, seed = seed)
resiVAE_mean$week_component <- FALSE

resiVAE_mean$max_season <- as.numeric(resiVAE_mean$max_season)
st_trend_test_time <- predict_coords_to_IC(resiVAE_mean, coords_time_validation[, 1:2], coords_time_validation[, 3] + temp_shift)
pred_obs_time_aux_mean <- predict(resiVAE_mean, st_trend_test_time, IC_to_data = TRUE)
preds_mean <- pred_obs_time_aux_mean #+ data_full_test[, seas_var_names]
mse_time_aux_mean2 <- apply((data_full_validation[, variables] - preds_mean)^2, 2, mean)
mse_time_aux_mean2
results[3, ] <- c(mean(mse_time_aux_mean2/sds^2), p, toString(spatial_basis), toString(temporal_basis), ar_order, "iVAE_trend")

save(results, file = paste0("athens_pred/", file_name))
results
