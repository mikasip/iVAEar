#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
N_sim <- as.integer(args[1])
n_layers <- as.integer(args[2])
nonlin <- args[3]
type <- "ar1_var_p6s6"
n_segments <- as.integer(args[4])
n_div <- as.integer(args[5])
start_seed <- as.integer(args[6])
reps_per_job <- as.integer(args[7])
add_trend <- as.integer(args[8])
trend_string <- ifelse(add_trend, "_", "_no_trend_")
filename <- paste0("res_", type, "_", n_layers, "_", nonlin, "_", n_segments, "_", n_div, trend_string, start_seed, ".RData")

.libPaths(c("/projappl/project_2003994/project_rpackages", .libPaths()))

library(gstat)
library(tensorflow)
library(keras)
library(RcppHungarian)
library(JADE)
library(SpaceTimeBSS)
library(fastICA)
library(NonlinearBSS)
library(keras3)
source("/scratch/project_2003994/STNICA/helpers.R")

c_matern <- function(h, range, nu) {
    if (h == 0) {
        return(1)
    }
    1 / (2^(nu - 1) * gamma(nu)) * (h / range)^nu * besselK(h / range, nu)
}

matern_covmat <- function(coords, range, nu, sigmas) {
    n_s <- dim(coords)[1]
    cov_mat <- matrix(NA, nrow = n_s, ncol = n_s)
    for (i in 1:n_s) {
        for (j in i:n_s) {
            h <- sqrt(sum((coords[i, ] - coords[j, ])^2))
            cov_val <- sigmas[i] * sigmas[j] * c_matern(h, range, nu)
            cov_mat[i, j] <- cov_val
            cov_mat[j, i] <- cov_val
        }
    }
    return(cov_mat)
}

trend <- function(coords_time, omega_x, omega_y, omega_t, phi, a) {
    trend <- apply(coords_time, 1, function(row) a * sin(omega_x*row[1] + omega_y*row[2] + omega_t*row[3] + phi))
    return(trend)
}

linear_trend <- function(coords_time, coef_x, coef_y, coef_t) {
    trend <- apply(coords_time, 1, function(row) row[1] * coef_x + row[2] * coef_y + row[3] * coef_t)
    return(trend)
}

n_s <- 100
n_t <- 500
p <- 6
s <- 6

params <- list(
    list(range = 0.2, nu = 0.5),
    list(range = 0.15, nu = 1),
    list(range = 0.1, nu = 0.25),
    list(range = 0.3, nu = 2),
    list(range = 0.05, nu = 0.75),
    list(range = 0.25, nu = 1.5)
)

shift_params <- list(
    list(range = 0.25, nu = 3),
    list(range = 0.15, nu = 2.5),
    list(range = 0.1, nu = 2),
    list(range = 0.3, nu = 1.5),
    list(range = 0.2, nu = 1),
    list(range = 0.05, nu = 0.5)
)

changing_AR1_coefs <- function(n, shift, scale) {
    coefs <- cos((2 * pi * c(1:n) * scale) / n - shift)
    return(coefs)
}

var_coefs <- function(n_t, coords, n_seg_time, n_seg_space) {
    var_mat <- matrix(NA, nrow = n_t, ncol = n_s)
    cent_coords <- coords[sample(1:nrow(coords), n_seg_space), ]
    cent_times <- c(1:n_t)[sample(1:n_t, n_seg_time)]
    labels_space <- unlist(apply(coords, 1, FUN = function(coord) {
        centered_points <- matrix(unlist(apply(cent_coords, 1,
            FUN = function(x) {
                x - coord
            }
        )), byrow = TRUE, nrow = n_seg_space)
        dists <- apply(centered_points, 1, FUN = function(x) {
            sqrt(x[1]^2 + x[2]^2)
        })
        return(which(dists == min(dists))[1])
    }))
    labels_time <- sapply(1:n_t, FUN = function(time) {
        dists <- abs(cent_times - time)
        return(which(dists == min(dists))[1])
    })
    for (k1 in 1:n_seg_space) {
        space_seg_inds <- which(labels_space == k1)
        for (k2 in 1:n_seg_time) {
            time_seg_inds <- which(labels_time == k2)
            var_mat[time_seg_inds, space_seg_inds] <- runif(1, 0.1, 5)
        }
    }
    return(var_mat)
}

epochs <- 60

res_df <- data.frame(matrix(ncol = 8, nrow = 3 * reps_per_job))
names(res_df) <- c("model", "MCC", "MD", "n_layers", "n_div", "seed", "nonlinearity", "model_layers")
ind <- 1
for (j in 1:reps_per_job) {
    seed <- start_seed + j * 10
    set.seed(seed)

    ar1_scale_params <- runif(p, 1, 10)

    coords <- matrix(runif(2 * n_s, 0, 1), ncol = 2)
    coords_time <- matrix(NA, nrow = n_s * n_t, ncol = 3)
    coords_time[, 1] <- rep(coords[, 1], n_t)
    coords_time[, 2] <- rep(coords[, 2], n_t)
    coords_time[, 3] <- rep(1:n_t, each = n_s)

    data <- array(NA, c(n_s * n_t, p))
    baseline_ar1s <- runif(p, 0.6, 0.99)
    for (i in 1:p) {
        shift_covmat <- matern_covmat(coords, shift_params[[i]]$range, shift_params[[i]]$nu, rep(1, n_s))
        shift_params_space <- t(rnorm(n_s, 0, 0.3)) %*% chol(shift_covmat)
        var_mat <- var_coefs(n_t, coords, 10, 5)
        cov_i <- matern_covmat(coords, params[[i]]$range, params[[i]]$nu, var_mat[1, ])
        st_data <- matrix(NA, ncol = n_s, nrow = n_t)
        st_data[1, ] <- rnorm(n_s) %*% chol(cov_i)
        t_seg <- 1
        baseline_ar1 <- baseline_ar1s[i]
        for (t_i in 2:n_t) {
            ar1_params <- cos((2 * pi * rep(t_i, n_s) * ar1_scale_params[i]) / n_t - shift_params_space)
            cov_i <- matern_covmat(coords, params[[i]]$range, params[[i]]$nu, var_mat[t_i, ])
            kernel_mat <- diag(c(ar1_params))
            sigma1 <- 1
            eps <- t(rnorm(n_s) %*% chol(sigma1^2 * cov_i))
            st_data[t_i, ] <- baseline_ar1 * kernel_mat %*% as.matrix(st_data[t_i - 1, ]) + eps
        }
        if (add_trend) {
            trend_i <- trend(coords_time, runif(1, 0.2, 4), runif(1, 0.2, 4), runif(1, 0.01, 0.1), runif(1, 0, 2*pi), runif(1, 0, 5)) + linear_trend(coords_time, runif(1, -3, 3), runif(1, -3, 3), runif(1, -0.01, 0.01))
        } else {
            trend_i <- 0
        }
        data[, i] <- as.vector(t(st_data)) + trend_i
    }

    data_scaled <- data / mean(sqrt(diag(var(data))))
    obs_data_obj <- mix_data_acyclic(data_scaled, s, 4, "elu")
    obs_data <- obs_data_obj$data
    graph1 <- obs_data_obj$edge_matrices[[1]]
    graph2 <- obs_data_obj$edge_matrices[[2]]

    
    resVAEar1_r <- iVAEar_radial(obs_data, coords_time[, 1:2],
        coords_time[, 3], p, n_s = n_s, ar_order = 1,
        hidden_units = c(128, 128, 128),
        aux_hidden_units = c(128, 128, 128), error_dist_sigma = 0.02, lr_start = 0.001,
        epochs = epochs, validation_split = 0, batch_size = 64, seed = seed
    )
    cormat2 <- cor(resVAEar1_r$IC, data)
    if (anyNA(cormat2)) {
        res_df[ind, ] <- c("iVAEar1_fail", 0, 1, 1, n_layers, n_div, seed, nonlin, 3)
    } else {
        MCCVAE <- absolute_mean_correlation(cormat2)
        MDVAE <- MD(cormat2, diag(p))
        res_df[ind, ] <- c("iVAEar1_r", MCCVAE, MDVAE, n_layers, n_div, seed, nonlin, 3)
    }
    ind <- ind + 1

    resVAEar1_r <- iVAEar_radial(obs_data, coords_time[, 1:2],
        coords_time[, 3], p, n_s = n_s, ar_order = 3,
        hidden_units = c(128, 128, 128),
        aux_hidden_units = c(128, 128, 128), error_dist_sigma = 0.02, lr_start = 0.001,
        epochs = epochs, validation_split = 0, batch_size = 64, seed = seed
    )
    cormat2 <- cor(resVAEar1_r$IC, data)
    if (anyNA(cormat2)) {
        res_df[ind, ] <- c("iVAEar3_fail", 0, 1, 1, n_layers, n_div, seed, nonlin, 3)
    } else {
        MCCVAE <- absolute_mean_correlation(cormat2)
        MDVAE <- MD(cormat2, diag(p))
        res_df[ind, ] <- c("iVAEar3_r", MCCVAE, MDVAE, n_layers, n_div, seed, nonlin, 3)
    }
    ind <- ind + 1

    resVAEar1_r <- iVAEar_radial(obs_data, coords_time[, 1:2],
        coords_time[, 3], p, n_s = n_s, ar_order = 5,
        hidden_units = c(128, 128, 128),
        aux_hidden_units = c(128, 128, 128), error_dist_sigma = 0.02, lr_start = 0.001,
        epochs = epochs, validation_split = 0, batch_size = 64, seed = seed
    )
    cormat2 <- cor(resVAEar1_r$IC, data)
    if (anyNA(cormat2)) {
        res_df[ind, ] <- c("iVAEar5_fail", 0, 1, 1, n_layers, n_div, seed, nonlin, 3)
    } else {
        MCCVAE <- absolute_mean_correlation(cormat2)
        MDVAE <- MD(cormat2, diag(p))
        res_df[ind, ] <- c("iVAEar5_r", MCCVAE, MDVAE, n_layers, n_div, seed, nonlin, 3)
    }
    ind <- ind + 1
}

res_test <- res_df

save(res_test, file = paste0("res_ar1_var_sim2/", filename))
res_test
