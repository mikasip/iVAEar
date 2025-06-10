library(NonlinearBSS)
library(tensorflow)
library(keras3)
library(lubridate)
library(ggplot2)
library(dplyr)
library(cutoffR)
library(sp)
library(spacetime)
library(sf)
library(gstat)
library(covatest)
library(forecast)
library(kernlab)

library(httpgd)
hgd()

#####################################
## HELPER FUNCTIONS #################
#####################################

fit_product_sum <- function(vv, res_var, verbose = FALSE) {
    # Temporal variogram estimation
    vv_t_temp <- data.frame(vv)
    vv_t_temp <- vv_t_temp %>% filter(spacelag == 0)
    vv_t <- data.frame(cbind(vv_t_temp$np, as.numeric(vv_t_temp$timelag), vv_t_temp$gamma))
    vv_t <- vv_t[-1, ]
    names(vv_t) <- c("np", "dist", "gamma")
    class(vv_t) <- c("gstatVariogram", "data.frame")
    vg_t <- fit.variogram(vv_t, vgm(0.3, c("Sph"), 20))
    if (verbose) {
        print(vg_t)
        plot(x = vv_t_temp$timelag, y = vv_t_temp$gamma, main = "Temporal variogram", pch = 19)
    }
    vgm_t <- vgm(
        psill = vg_t$psill[1], model = vg_t$model[1],
        range = vg_t$range[1], nugget = 0,
        kappa = vg_t$kappa[1]
    )
    tm_mod <- variogramLine(
        vgm_t,
        dist_vector = seq(0, 150)
    )
    if (verbose) lines(x = tm_mod[, 1], y = tm_mod[, 2], type = "l", col = "black")

    # Spatial variogram estimation
    vv_s_temp <- data.frame(vv)
    vv_s_temp <- vv_s_temp %>% filter(timelag == 0)
    vv_s <- data.frame(cbind(vv_s_temp$np, as.numeric(vv_s_temp$spacelag), vv_s_temp$gamma))
    vv_s <- vv_s[-1, ]
    names(vv_s) <- c("np", "dist", "gamma")
    class(vv_s) <- c("gstatVariogram", "data.frame")
    vg_s <- fit.variogram(vv_s, vgm(0.3, c("Sph"), 10000))
    if (verbose) {
        vg_s
        plot(x = vv_s_temp$spacelag, y = vv_s_temp$gamma, main = "Spatial variogram", pch = 19)
    }
    vgm_s <- vgm(
        psill = vg_s$psill[1], model = vg_s$model[1],
        range = vg_s$range[1], nugget = 0,
        kappa = vg_s$kappa[1]
    )
    s_mod <- variogramLine(
        vgm_s,
        dist_vector = seq(0, 20)
    )
    if (verbose) lines(x = s_mod[, 1], y = s_mod[, 2], type = "l", col = "black")

    # compute the k parameter for the PRODUCT-SUM MODEL
    sillT <- res_var - vg_t$psill[1] # c00-sillsp
    if (sillT < 0) sillT <- 0.01
    sillS <- res_var - vg_s$psill[1] # c00-sillt
    # k=(sillST-sillS-sillT)/(sillS*sillT)
    print(sillT)
    print(sillS)
    k <- abs((res_var - sillS - sillT) / (sillS * sillT))
    print(k)

    # fit the PRODUCT-SUM MODEL
    vgm_prodSum_model <- vgmST("productSum", space = vgm_s, time = vgm_t, k = k)

    vgm_prodSum_fit <- fit.StVariogram(vv, vgm_prodSum_model)

    return(list(vgm_s = vgm_s, vgm_t = vgm_t, vgm_st = vgm_prodSum_fit))
}

create_st_dataset <- function(data, var_names, length.out) {
    sp <- cbind(data$X_cent, data$Y_cent)
    sp <- unique(sp)
    sp.names <- unique(data$station_name)
    colnames(sp) <- c("x", "y")
    startDate <- as.Date(min(data$Date))
    sp2 <- sp::SpatialPoints(sp)
    row.names(sp2) <- sp.names
    time <- seq.Date(from = startDate, by = "day", length.out = length.out)
    ordered_df <- data[order(data[, "Date"]), ]
    stfdf <- STFDF(sp = sp2, time = time, data = data.frame("var" = ordered_df[, var_names]))
    return(stfdf)
}

variables <- c("Wind.Speed..U.", "Wind.Speed..V.", "Dewpoint.Temp", "Soil.Temp", "Temp", "Relative.Humidity", "PM10", "PM2.5", "NO2", "O3")
seas_var_names <- sapply(variables, function(name) paste0(name, "_seas"))
res_var_names <- sapply(variables, function(name) paste0(name, "_res"))

data_full_train <- read.csv("athens_daily_train.csv")
head(data_full_train)
data_full_test <- read.csv("athens_daily_test.csv")
coords_time_train <- read.csv("athens_coords_train.csv")
head(coords_time_train)
coords_time_train <- coords_time_train[, 2:4]
coords_time_test <- read.csv("athens_coords_test.csv")
coords_time_test <- coords_time_test[, 2:4]
head(coords_time_train)

n_s <- nrow(unique(coords_time_train[, 1:2]))

# Number of IC estimation:
# The results may vary slightly based on seed
resiVAEar1_p5 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3], 5, n_s = n_s, 
    epochs = 60, batch_size = 64)
resiVAEar1_p6 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3], 6, n_s = n_s, 
    epochs = 60, batch_size = 64)
resiVAEar1_p7 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3], 7, n_s = n_s, 
    epochs = 60, batch_size = 64)
resiVAEar1_p8 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3], 8, n_s = n_s, 
    epochs = 60, batch_size = 64)
resiVAEar1_p9 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3], 9, n_s = n_s, 
    epochs = 60, batch_size = 64)
resiVAEar1_p10 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3], 10, n_s = n_s, 
    epochs = 60, batch_size = 64)

pAICs <- c(-resiVAEar1_p5$elbo + 5,
          -resiVAEar1_p6$elbo + 6,
          -resiVAEar1_p7$elbo + 7,
          -resiVAEar1_p8$elbo + 8,
          -resiVAEar1_p9$elbo + 9,
          -resiVAEar1_p10$elbo + 10)
p <- which(pAICs == min(pAICs)) + 4


elbo_df <- data.frame(acc = c(as.numeric(resiVAEar1_p5$elbo), as.numeric(resiVAEar1_p6$elbo), as.numeric(resiVAEar1_p7$elbo), as.numeric(resiVAEar1_p8$elbo), as.numeric(resiVAEar1_p9$elbo), as.numeric(resiVAEar1_p10$elbo)), dim = c(5:10))
ggplot(elbo_df, aes(x = dim, y = acc)) +
  geom_line() +
  geom_point(aes(color = (dim != p))) +
  theme_minimal() +
  xlab("The latent dimension") +
  ylab("ELBO") +
  scale_x_continuous(breaks = seq(5, 9, by = 1)) +
  scale_color_manual(values = c("red", "black")) +
  theme(legend.position = "none")

# Parameters based on the validation results
p <- 9
spatial_basis <- c(2, 9)
temporal_basis <- c(9, 17)
seed <- 29012025
temp_shift <- 80
ar_order <- 2
n_s <- 58

sds <- apply(data_full_train[, res_var_names], 2, sd)

resiVAEar1 <- iVAEar_radial(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), coords_time_train[, 3] + temp_shift, p, n_s = n_s, 
    spatial_basis = spatial_basis, temporal_basis = temporal_basis, seasonal_period = 365, ar_order = ar_order,
    aux_hidden_units = c(16), epochs = 40, batch_size = 64, error_dist_sigma = 0.02, seed = seed,
    lr_start = 0.0001, lr_end = 0.0001)

last_coords_time <- coords_time_train[which(coords_time_train[, 3] %in% ((1100 - (ar_order - 1)):1100)), ]
st_ar_test_time <- predict_coords_to_IC_ar(resiVAEar1, last_coords_time[, 1:2], last_coords_time[, 3] + temp_shift, NULL, coords_time_test[, 1:2], coords_time_test[, 3] + temp_shift, NULL)
st_ar_test_time <- st_ar_test_time$preds
pred_obs_time_aux <- predict(resiVAEar1, st_ar_test_time, IC_to_data = TRUE)
mse_time_aux <- apply((data_full_test[, variables] - pred_obs_time_aux)^2, 2, mean)
mean(mse_time_aux/sds^2)

coords_time_shifted <- coords_time_train
coords_time_shifted[, 3] <- coords_time_shifted[, 3] + temp_shift
resiVAEar1_s <- iVAEar_segmentation(as.matrix(data_full_train[, variables]), as.matrix(coords_time_shifted), c(5000, 5000, 5), c(1, 1, 2), p, n_s = n_s, 
    ar_order = 3, time_dim = 3, seasonal_period = 365,
    aux_hidden_units = c(16), epochs = 40, batch_size = 64, error_dist_sigma = 0.02, seed = seed,
    lr_start = 0.0001, lr_end = 0.0001)

st_ar_test_time2 <- predict_coords_to_IC_ar(resiVAEar1_s, last_coords_time[, 1:2], last_coords_time[, 3] + temp_shift, NULL, coords_time_test[, 1:2], coords_time_test[, 3] + temp_shift, NULL)
st_ar_test_time2 <- st_ar_test_time2$preds
pred_obs_time_aux2 <- predict(resiVAEar1, st_ar_test_time2, IC_to_data = TRUE)
mse_time_aux2 <- apply((data_full_test[, variables] - pred_obs_time_aux2)^2, 2, mean)
mean(mse_time_aux2/sds^2)

# iVAE radial (without AR)
resiVAE_mean <- iVAE_radial_spatio_temporal(as.matrix(data_full_train[, variables]), as.matrix(coords_time_train[, 1:2]), 
    coords_time_train[, 3] + temp_shift, p, spatial_basis = spatial_basis, temporal_basis = temporal_basis, aux_hidden_units = c(16),
    seasonal_period = 365, epochs = 40, batch_size = 64, seed = seed, error_dist_sigma = 0.02, lr_start = 0.0001, lr_end = 0.0001)

resiVAE_mean$max_season <- as.numeric(resiVAE_mean$max_season)
resiVAE_mean$week_component <- FALSE
st_trend_test_time <- predict_coords_to_IC(resiVAE_mean, coords_time_test[, 1:2], coords_time_test[, 3] + temp_shift)
pred_obs_time_aux_mean <- predict(resiVAE_mean, st_trend_test_time, IC_to_data = TRUE)
mse_time_aux_mean <- apply((data_full_test[, variables] - pred_obs_time_aux_mean)^2, 2, mean)
mse_time_aux_mean
mean(mse_time_aux_mean/sds^2)

# ARIMA 
library(forecast)
preds_arima <- data_full_test[, c("station_name", "Date_numeric")]
for (var in res_var_names) {
    preds_arima[, var] <- NA
}

for (station in unique(data_full_train$station_name)) {
    for (variable in res_var_names) {
        station_idxs <- which(data_full_train$station_name == station)
        temp_data <- data_full_train[station_idxs, variable]
        model <- auto.arima(temp_data)
        temp_preds <- forecast:::predict.default(model, 24)$mean
        station_idxs_test <- which(preds_arima$station_name == station)
        preds_arima[station_idxs_test, variable] <- as.numeric(temp_preds)
    }   
}

preds_arima_seas <- preds_arima[, res_var_names] + data_full_test[, seas_var_names]
arima_mse <- apply((as.matrix(preds_arima_seas[, res_var_names]) - data_full_test[, variables])^2, 2, mean)
mean(arima_mse/sds^2)


# VARIMA

# Initialize predictions dataframe
preds_varima <- data_full_test[, c("station_name", "Date_numeric")]
for (var in res_var_names) {
    preds_varima[, var] <- NA
}

# Fit multivariate models by station
for (station in unique(data_full_train$station_name)) {
    station_idxs <- which(data_full_train$station_name == station)
    temp_data <- as.matrix(data_full_train[station_idxs, res_var_names])

    temp_integrated <- temp_data[-1, ]
    for(var in res_var_names) {
      temp_integrated[, var] <- diff(temp_data[, var], differences = 1)
    }
    lag_select <- vars::VARselect(temp_integrated[, res_var_names], lag.max = 8, type = "both")
    optimal_lag <- lag_select$selection["AIC(n)"]
    optimal_aic <- lag_select$criteria["AIC(n)", optimal_lag]
    tryCatch({
        model <- MTS::VARMA(temp_integrated[, res_var_names], p = 1, q = 1)
        if (is.null(model)) {
            stop("Model fitting failed, falling back to VAR")
        }
        if (model$aic < optimal_aic) {
            optimal_lag <- 1
            optimal_aic <- model$aic
            temp_preds <- MTS::VARMApred(model, h = 24)$pred
        } else {
            model <- vars::VAR(temp_integrated[, res_var_names], p = optimal_lag, type = "both")
            temp_preds <- predict(model, n.ahead = 24)$fcst
            temp_preds <- sapply(temp_preds, function(x) x[,"fcst"])
        }
    }, error = function(e) {
        model <- vars::VAR(temp_integrated[, res_var_names], p = optimal_lag, type = "both")
        temp_preds <- predict(model, n.ahead = 24)$fcst
        temp_preds <- sapply(temp_preds, function(x) x[,"fcst"])
    })
    # Apply integration if needed
    level_preds <- matrix(0, nrow = nrow(temp_preds), ncol = ncol(temp_preds))
    colnames(level_preds) <- res_var_names
    last_observed <- temp_data[nrow(temp_data), ]
    for(i in 1:length(res_var_names)) {
      var <- res_var_names[i]
      level_preds[, i] <- last_observed[i] + cumsum(temp_preds[, i])
    }
    # Store predictions
    station_idxs_test <- which(preds_varima$station_name == station)
    for (i in seq_along(res_var_names)) {
        preds_varima[station_idxs_test, res_var_names[i]] <- level_preds[, i]
    }
}

# KRIGING

for (var_name in res_var_names) {
    var_df <- data_full_train[, c("station_name", "X_cent", "Y_cent", "Date_numeric", "Date", var_name)]
    names(var_df) <- c("station_name", "X_cent", "Y_cent", "Date_numeric", "Date", "var")
    stfdf_data <- create_st_dataset(var_df, "var", n_t_train)
    vv_var <- variogramST(var ~ 1, stfdf_data, tlags = 0:20, cutoff = 30000)
    vv_var[1, 2:3] <- 0
    plot(vv_var, wireframe = TRUE, main = "var")
    var_res <- var(stfdf_data[, , "var"]@data[[1]], na.rm = TRUE)
    vg_obj <- fit_product_sum(vv_var, var_res, verbose = TRUE)
    
    pred_var_time <- krigeST(var ~ 1,
        data = stfdf_data,
        newdata = stfdf_validation_time, vg_obj$vgm_st, nmax = 40,
        stAni = 5000 / 7
    )

    pred_var_df_time <- as.data.frame(pred_var_time)
    pred_var_df_time <- pred_var_df_time[order(pred_var_df_time[, "sp.ID"]), ]
    data_full_test[, var_name] <- pred_var_df_time[, "var1.pred"]
}

pred_obs_time <- data_full_test[, res_var_names] + data_full_test[, seas_var_names]
pred_obs_time_scaled <- sweep(pred_obs_time, 2, sds_test_time, "/")

mse_pred_time <- apply((data_full_test[, ] - pred_obs_time)^2, 2, mean)
mse_pred_time
xtable(rbind(mse_time_aux, mse_pred_time, arima_mse))
mse_time_aux
mean(mse_pred_time/sds^2)
mean(mse_time_aux/sds^2)
