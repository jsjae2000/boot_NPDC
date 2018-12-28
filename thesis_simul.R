# 1. settings ----------------------------------------------------------------
## library
library(dplyr)
library(ggplot2)

## function
"%++%" <- function(x, y) paste0(x, y)

f1 <- function(x, y_1) {
  beta0 <- -1
  beta1 <- -.5
  beta2 <- .5
  
  out <- beta0 + beta1*x + beta2 * y_1
  out <- unname(out)
  return(out)
}

x_kern <- function(bw, x, X, type = "epane") {
  ### input type
  # bw : vector s.t. length(bw) == d | length(bw) == 1
  # x : vector s.t. length(x) == d
  # X : data frame s.t. nrow(X) == d
  # type : one of c("epane")
  # admit the case of common bandwidths
  if (length(bw) == 1) bw <- rep(bw, ncol(X))
  
  base <- sweep(X, 2, x) %>%
    sweep(2, FUN = "/", bw)
  (switch(type,
          epane = 3/4 * (1 - base**2)) %>%
      sweep(2, FUN = "/", bw) * (abs(base) <= 1)) %>%
    apply(1, prod) %>%
    as.double()
}

z_kern <- function(bw, z, Z) {
  ### input type
  # bw : vector s.t. length(bw) == d | length(bw) == 1
  # z : vector s.t. length(z) == d
  # Z : data frame s.t. nrow(Z) == d
  
  # ### toy example
  # Z <- data.frame(c(1,1,1,2,2,2,3,3,3), c(1,1,1,2,2,2,3,3,3))
  # z <- c(3,1)
  # bw <- c(2, 3)
  
  # admit the case of common bandwidths
  if (length(bw) == 1) bw <- rep(bw, ncol(Z))
  apply(Z, 1, function(Z_vec) prod(bw[z != Z_vec]))
}

cmn_preproc <- function(X, Y, d_lag) {
  ### input
  # X : data frame of covariates which allowed to contain continuous and discrete variables
  # Y : response vector(0/1)
  
  factor2double <- function(x) as.double(as.character(x))
  Y <- data.frame(resp__var__ = Y)
  ### check dimension of X and Y ###
  if (nrow(X) != nrow(Y)) stop("check dimension!")
  
  X_isfactor <- unlist(lapply(X, is.factor))
  conti_str <- names(which(!X_isfactor))
  disc_str <- names(which(X_isfactor))
  
  X <- mutate_if(X, is.factor, factor2double)
  df_lag <- cbind(Y, X)
  
  if (d_lag > 0) {
    for (i_lag in 1:d_lag) {
      df_lag[["resp__lagged__" %++% i_lag]] <- dplyr::lag(df_lag[["resp__var__"]], i_lag)  
    }
    df_lag <- df_lag[-(1:d_lag), ]
    disc_str <- union(disc_str, "resp__lagged__" %++% 1:d_lag)
  }
  
  attr(df_lag, "conti_str") <- conti_str
  attr(df_lag, "disc_str") <- disc_str
  
  # rare case
  if (any(duplicated(names(df_lag)))) stop("colnames in X should not contain c('resp__var__', 'resp__lagged__1',..., 'resp__lagged__d_lag'")
  return(df_lag)
}

cmn_defftns <- function(link_ftn, V_ftn) {
  if (link_ftn == "probit") {
    ginv <- pnorm
    gprime <- function(x) 1 / dnorm(qnorm(x))
    gprime2 <- function(x) (x * dnorm(qnorm(x))) * (dnorm(x)) / (dnorm(qnorm(x)))^2
  }
  if (V_ftn == "binomial") {
    Veval <- function(x) x * (1-x)
    Vprime <- function(x) 1 - 2*x
  }
  
  return(list(ginv = ginv,
              gprime = gprime,
              gprime2 = gprime2,
              Veval = Veval,
              Vprime = Vprime))
}

MNQLE_est <- function(X, Y, d_lag, point, bw, link_ftn = "probit", V_ftn = "binomial", method = "NR", init = NULL, tolerance = 1e-05, max_iter = 1000) {
  if (FALSE) {
    X = df_analysis[c("x_1")]
    Y = df_analysis[, "y"]
    d_lag <- 1
    
    point <- c(0.5, 1)
    bw = unname(c(n^{-1/5} * 1.06 * apply(X, 2, sd), n**-.4))
    
    link_ftn <- "probit"
    V_ftn <- "binomial"
    
    method = "NR"
    init = NULL
    tolerance = 1e-05
    max_iter = 1000
  }
  
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(link_ftn, V_ftn)
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  
  i_conti <- match(attr(df_lag, "conti_str"), names(df_lag)[-1])
  i_disc <- match(attr(df_lag, "disc_str"), names(df_lag)[-1])
  x_est <- point[i_conti]
  z_est <- point[i_disc]
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  return(MNQLE_optimize(x_est, z_est, x_bw, z_bw, ftns, df_lag, method, init, tolerance, max_iter))
}

MNQLE_bootstrap1 <- function(n_boot = 500, X, Y, estimates, d_lag, point, bw, link_ftn = "probit", V_ftn = "binomial", 
                             method = "NR", init = NULL, tolerance = 1e-05, max_iter = 1000, save_sample = TRUE) {
  if (FALSE) {
    estimates <- df_simul$est
    
    X = df_analysis[c("x_1")]
    Y = df_analysis[, "y"]
    d_lag <- 1
    
    point <- c(0.5, 1)
    bw = unname(c(n^{-1/5} * 1.06 * apply(X, 2, sd), n**-.4))
    
    link_ftn <- "probit"
    V_ftn <- "binomial"
    
    method = "NR"
    init = NULL
    tolerance = 1e-05
    max_iter = 1000
  }
  
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(link_ftn, V_ftn)
  if (d_lag > 0) estimates <- tail(estimates, -d_lag)
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  
  i_conti <- match(attr(df_lag, "conti_str"), names(df_lag)[-1])
  i_disc <- match(attr(df_lag, "disc_str"), names(df_lag)[-1])
  x_est <- point[i_conti]
  z_est <- point[i_disc]
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  if (nrow(df_lag) != length(estimates)) stop("check dimension of X and estimates")
  
  
  ### bootstrap ###
  X_boot <- subset(df_lag, select =  !(names(df_lag) %in% c("resp__var__", "resp__lagged__" %++% 1:d_lag)))
  
  # prepare to save bootstrap samples
  if (save_sample) boot_sample <- matrix(NA_integer_, nrow = length(estimates), ncol = n_boot)
  est_boot <- c()
  for (i_boot in 1:n_boot) {
    # generate bootstrap sample  
    Y_boot <- c()
    for (i_sample in 1:length(estimates)) {
      Y_boot[i_sample] <- rbinom(1, 1, estimates[i_sample])
    }
    if (save_sample) boot_sample[, i_boot] <- Y_boot
    
    # calculate bootstrap estimates
    df_lag_boot <- cmn_preproc(X_boot, Y_boot, d_lag)
    est_boot[i_boot] <- MNQLE_optimize(x_est, z_est, x_bw, z_bw, ftns, df_lag_boot, method, init, tolerance, max_iter)
  }
  
  if (save_sample) attr(est_boot, "samples") <- boot_sample
  return(est_boot)
}



MNQLE_optimize <- function(x_est, z_est, x_bw, z_bw, ftns, df_lag, method = "NR", init = NULL, tolerance = 1e-05, max_iter = 1000) {
  conti_str <- attr(df_lag, "conti_str")
  disc_str <- attr(df_lag, "disc_str")
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  
  ### start optimization ###
  if (method == "NR") {
    ## Newton Raphson algorithm ##
    f_new <- matrix(init, ncol = 1)
    f_old <- f_new + tolerance * 100
    i_iter <- 0
    while (norm(f_old - f_new, type = "F") > tolerance) {
      i_iter <- i_iter + 1
      if (i_iter > max_iter) break
      f_old <- f_new
      NR_cal <- MNQLE_NR_calculate(f_old, x_est, z_est, x_bw, z_bw, ftns, df_lag)
      
      if (abs(det(NR_cal$eq2)) <1e-05) {
        f_new <- f_old - solve(NR_cal$eq2 + diag(diag(NR_cal$eq2), 2) * .1) %*% NR_cal$eq1  
      } else {
        f_new <- f_old - solve(NR_cal$eq2) %*% NR_cal$eq1  
      }
    }
  }
  
  return(ginv(f_new[1]))
}

MNQLE_NR_calculate <- function(f_at, x_est, z_est, x_bw, z_bw, ftns, df_lag) {
  if (FALSE) {
    f_at = f_old
    x_est = point[1:length(conti_str)]
    z_est = point[-(1:length(conti_str))] 
    x_bw = bw[1:length(conti_str)]
    z_bw = bw[-(1:length(conti_str))]
    ftns = ftns
    df_lag = df_lag
  }
  conti_str <- attr(df_lag, "conti_str")
  disc_str <- attr(df_lag, "disc_str")
  
  wci_all <- x_kern(x_bw, x_est, df_lag[conti_str])
  cal_ind <- wci_all != 0
  wci <- wci_all[cal_ind]
  df_cal <- df_lag[cal_ind, ]
  
  wdi <- z_kern(z_bw, z_est, df_cal[disc_str])
  mi <- ftns$ginv(data.matrix(sweep(df_cal[conti_str], 2, x_est)) %*%
                    matrix(f_at[-1], ncol = 1) +
                    f_at[1]
  )[, 1]
  mi[mi >= 1 - 1e-05] <- 1 - 1e-05
  mi[mi <= 1e-05] <- 1e-05
  ### calculate eq1 val ###
  q1i <- (df_cal[, 1] - mi) / (ftns$Veval(mi) * ftns$gprime(mi))
  product1_i <- matrix(wci * wdi * q1i, nrow = 1)
  # dim(Wi) = n X d
  Wi <- cbind(intercept = 1, sweep(df_cal[conti_str], 2, x_est) %>% sweep(2, FUN = "/", x_bw)) %>%
    as.matrix()
  # eq1_val
  eq1_val <- (product1_i %*% Wi) / nrow(df_lag)
  
  ### calculate eq2 val ###
  q2i_tmp11 <- ftns$Vprime(mi) / (ftns$Veval(mi) * (ftns$gprime(mi)^2))
  q2i_tmp12 <- ftns$gprime2(mi) / (ftns$gprime(mi))^3
  q2i_tmp1 <- ftns$Veval(mi) / (q2i_tmp11 + q2i_tmp12)
  q2i_tmp2 <- ftns$Veval(mi) * ftns$gprime(mi)^2; q2i_tmp2 <- 1 / q2i_tmp2
  q2i <- -(df_cal[, 1] - mi) / q2i_tmp1 - q2i_tmp2
  
  product2_i <- matrix(wci * wdi * q2i, nrow = 1)
  
  eq2_val <- (sweep(t(Wi), 2, product2_i, FUN = "*") %*% Wi) / nrow(df_lag)
  
  return(list(eq1 = matrix(eq1_val, ncol = 1), eq2 = eq2_val))
}



# 2. generate data -----------------------------------------------------------
set.seed(111)
## 2-1. Example 1 ----------------------------------------------------------
n <- 501
y_init <- 0
x_bd <- c(-3, 3)
link_ftn <- "probit"
if (link_ftn == "probit") {
  ginv <- pnorm
  gprime <- function(x) 1 / dnorm(qnorm(x))
  gprime2 <- function(x) (x * dnorm(qnorm(x))) * (dnorm(x)) / (dnorm(qnorm(x)))^2
}

V_ftn <- "binomial"

y_vec <- c("0" = y_init)
mean_vec <- double()
x_vec = runif(n, min = x_bd[1], max = x_bd[2])
for (i in 1:n) {
  mean_val <- f1(x_vec[i], 
                 y_vec[as.character(i-1)]) %>%
    ginv
  mean_vec <- c(mean_vec, mean_val)
  y_vec <- c(y_vec, rbinom(1, 1, mean_val))
  names(y_vec)[length(y_vec)] <- as.character(i)
}

df_simul <- data.frame(x_1 = x_vec, 
                       z_1 = y_vec[-length(y_vec)],
                       y = y_vec[as.character(1:n)], 
                       mean = mean_vec)
df_analysis <- df_simul %>%
  select(x_1, y)

# check
plot(df_simul$x_1,
     df_simul$mean, 
     col = df_simul$z_1 + 1)



# 3. Estimation -----------------------------------------------------------
## inpt : define now
conti_str <- c("x_1")
disc_str <- character(0)
d_lag <- 1
x_bw_init <- n^{-1/5} * 1.06 * apply(df_analysis[conti_str], 2, sd)
z_bw_init <- n**-.4

df_simul$est <- NA
for (ii in 1:nrow(df_simul)) {
  set.seed(100)
  df_simul[ii, "est"] <- 
    MNQLE_est(X = df_analysis[conti_str],
              Y = df_analysis$y,
              d_lag = 1,
              point = c(df_simul[ii, "x_1"], df_simul[ii, "z_1"]),
              bw = c(x_bw_init, z_bw_init))
}

par(mfrow = c(1, 3))
plot(df_simul$x_1,
     df_simul$mean, 
     col = df_simul$z_1 + 1,
     ylim = 0:1)
plot(df_simul$x_1,
     df_simul$est, 
     col = df_simul$z_1 + 1,
     ylim = 0:1)
plot(df_simul$x_1,
     df_simul$est1, 
     col = df_simul$z_1 + 1,
     ylim = 0:1)
plot(df_simul$x_1,
     df_simul$est2, 
     col = df_simul$z_1 + 1,
     ylim = 0:1)

par(mfrow = c(1, 1))
plot(df_simul$mean - df_simul$est)
summary(df_simul$mean - df_simul$est)



# 4. bootstrap -----------------------------------------------------------
### check : coverage of pointwise confidence interval 
# calculate coverage for each points
n_boot <- 500
CI_alpha <- .05
n_coverage <- 200
df_CI <- df_simul[-1, ]

df_CI$count_coverage <- 0
cat("starts at", as.character(Sys.time()))
for (i_CI in 1:nrow(df_CI)) {
  if (i_CI %% (nrow(df_CI) %/% 10) == 0) cat("=")
  for (i_coverage in 1:n_coverage) {
    est_boot <- MNQLE_bootstrap1(n_boot = n_boot,
                                 X = df_analysis[conti_str], Y = df_analysis$y, 
                                 estimates = df_simul$est, d_lag = 1, 
                                 point = c(df_simul[i_CI, "x_1"], df_simul[i_CI, "z_1"]), 
                                 bw = c(x_bw_init, z_bw_init), 
                                 link_ftn = "probit", V_ftn = "binomial", 
                                 method = "NR", init = NULL, tolerance = 1e-05, max_iter = 1000, save_sample = FALSE)
    est_boot <- qnorm(est_boot)
    CI_boot <- quantile(est_boot, c(CI_alpha / 2, 1 - CI_alpha / 2))
    
    mean_true <- qnorm(df_CI[i_CI, "mean"])
    if (mean_true >= CI_boot[1] & mean_true <= CI_boot[2]) df_CI[i_CI, "count_coverage"] <- df_CI[i_CI, "count_coverage"] + 1
  }
}
cat("ends at", as.character(Sys.time()))