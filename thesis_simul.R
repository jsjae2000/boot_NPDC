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

x_kern <- function(bw, x, X, type = "epane"){
  base <- (x - X) / bw
  switch(type,
         epane = 3/4 * (1 - base)**2) / bw * (abs(base) <= 1)
}

z_kern <- function(bw, z, Z) {
  ifelse(z == Z, 1, bw)
}

QLNE_GD <- function(x_est, z_est, x_bw, z_bw, link_ftn, V_ftn, conti_str, disc_str, d_lag, y_str, df_analysis) {
  d <- length(conti_str)
  f_init <- rnorm(1 + d)
  if (link_ftn == "probit") {
    ginv <- pnorm
    gprime <- function(x) 1 / dnorm(qnorm(x))
    gprime2 <- function(x) (x * dnorm(qnorm(x))) * (dnorm(x)) / (dnorm(qnorm(x)))^2
  }
  
  
  if (length(x_est) != d) stop("dimension!")
  if (length(z_est) != length(disc_str) + d_lag) stop("dimension!")
  
  for (i_lag in 1:d_lag) {
    df_analysis[["y_" %++% i_lag]] <- dplyr::lag(df_analysis[[y_str]], i_lag)  
  }
  df_analysis <- na.omit(df_analysis)
  n <- nrow(df_analysis)
  disc_str2 <- union(disc_str, "y_" %++% 1:d_lag)
  
  f_grid <- expand.grid(f_0 = seq(-2, 2, length.out = 0), f_1 = seq(-2, 2, length.out = 50))
  for (ii in 1:nrow(f_grid)) {
    f_eval <- f_grid[ii, 1:(1+d)] %>% as.numeric
    mi <- ginv(data.matrix(sweep(df_analysis[conti_str], 2, x_est)) %*% 
                 matrix(f_eval[-1], ncol = 1) + 
                 f_eval[1]
    )[, 1]
    #if (ii %% 87 == 0) cat("grid search at", ii, mean(mi), '\n')
    obj_tmp <- x_kern(x_bw, x_est, df_analysis[[conti_str]]) * 
      #z_kern(z_bw, z_est, df_analysis[[disc_str2]]) * 
      (df_analysis[[y_str]] - mi) / (V_ftn(mi) * gprime(mi))
    f_grid$obj[[ii]] <- cbind(1, sweep(df_analysis[conti_str], 2, x_est) / x_bw ) %>% 
      sweep(., 1, obj_tmp, "*") %>% 
      apply(2, mean) %>%
      as.matrix() %>%
      norm(type = "F")
  }
  
  opt_ii <- which.min(abs(f_grid$obj))
  return(as.numeric(f_grid[opt_ii, 1:(1+d)]))
}

cmn_preproc <- function(X, Y, d_lag) {
  ### input
  # X : data frame of covariates which allowed to contain continuous and discrete variables
  # Y : response vector(0/1)
  
  
  Y <- data.frame(Y = Y)
  ### check dimension of X and Y ###
  if (nrow(X) != nrow(Y)) stop("check dimension!")
  
  X_isfactor <- unlist(lapply(X, is.factor))
  conti_str <- names(which(!X_isfactor))
  disc_str <- names(which(X_isfactor))
  if (is.null(init)) init <- rnorm(1 + length(conti_str))

  factor2double <- function(x) as.double(as.character(x))
  X <- mutate_if(X, is.factor, factor2double)
  df_lag <- cbind(Y, X)
  
  if (d_lag > 0) {
    for (i_lag in 1:d_lag) {
      df_lag[["Y_" %++% i_lag]] <- dplyr::lag(df_lag[["Y"]], i_lag)  
    }
    df_lag <- df_lag[-(1:d_lag), ]
    disc_str <- union(disc_str, "Y_" %++% 1:d_lag)
  }
  
  attr(df_lag, "conti_str") <- conti_str
  attr(df_lag, "disc_str") <- disc_str
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

MNQLE_NR_calculate <- function(f_at, x_est, z_est, x_bw, z_bw, link_ftn, V_ftn, conti_str, disc_str, y_str, df_lag) {
  wci_all <- x_kern(x_bw, x_est, df_lag[[conti_str]])
  cal_ind <- wci_all != 0
  wci <- wci_all[cal_ind]
  df_cal <- df_lag[cal_ind, ]
  
  wdi <- z_kern(z_bw, z_est, df_cal[[disc_str]])
  mi <- ginv(data.matrix(sweep(df_cal[conti_str], 2, x_est)) %*% 
               matrix(f_at[-1], ncol = 1) + 
               f_at[1]
  )[, 1]
  mi[mi >= 1 - 1e-02] <- 1 - 1e-02
  mi[mi <= 1e-02] <- 1e-02
  ### calculate eq1 val ###
  q1i <- (df_cal[[y_str]] - mi) / (Veval(mi) * gprime(mi))  
  product1_i <- matrix(wci * wdi * q1i, nrow = 1)
  # dim(Wi) = n X d
  Wi <- cbind(intercept = 1, sweep(df_cal[conti_str], 2, x_est) %>% sweep(2, FUN = "/", x_bw)) %>%
    as.matrix()
  # eq1_val
  eq1_val <- (product1_i %*% Wi) / nrow(df_lag)
  
  ### calculate eq2 val ###
  q2i_tmp11 <- Vprime(mi) / (Veval(mi) * (gprime(mi)^2))
  q2i_tmp12 <- gprime2(mi) / (gprime(mi))^3
  q2i_tmp1 <- Veval(mi) / (q2i_tmp11 + q2i_tmp12)
  q2i_tmp2 <- Veval(mi) * gprime(mi)^2; q2i_tmp2 <- 1 / q2i_tmp2
  q2i <- -(df_cal[[y_str]] - mi) / q2i_tmp1 - q2i_tmp2
  
  product2_i <- matrix(wci * wdi * q2i, nrow = 1)
  
  eq2_val <- (sweep(t(Wi), 2, product2_i, FUN = "*") %*% Wi) / nrow(df_lag)
  
  return(list(eq1 = matrix(eq1_val, ncol = 1), eq2 = eq2_val))
}

MNQLE_optimize <- function(X, Y, init = NULL, d_lag, point_est, bw, link_ftn, V_ftn, NR_tolerance = 1e-05, NR_max_iter = 1000, method = "NR") {
  if (FALSE) {
    X = df_analysis[c("x_1")]
    Y = df_analysis[, "y"]
    n <- nrow(df_analysis)
    point_est <- c(0.5, 1)
    bw = unname(c(n^{-1/5} * 1.06 * apply(X, 2, sd), n**-.4))
    
    link_ftn <- "probit"
    V_ftn <- "binomial"
    
    d_lag <- 1
    NR_tolerance = 1e-05
    NR_max_iter = 1000
    method = "BFGS"
    init = NULL
  }
  
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(link_ftn, V_ftn)
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  if (FALSE) {
    # Y <- data.frame(Y = Y)
    # 
    # ### check ncol(X) and length(bw) ###
    # if (ncol(X) + d_lag != length(bw)) stop("check dimension!")
    # if (nrow(X) != nrow(Y)) stop("check dimension!")
    # 
    # X_isfactor <- unlist(lapply(X, is.factor))
    # conti_str <- names(which(!X_isfactor))
    # disc_str <- names(which(X_isfactor))
    # if (is.null(init)) init <- rnorm(1 + length(conti_str))
    # if (link_ftn == "probit") {
    #   ginv <- pnorm
    #   gprime <- function(x) 1 / dnorm(qnorm(x))
    #   gprime2 <- function(x) (x * dnorm(qnorm(x))) * (dnorm(x)) / (dnorm(qnorm(x)))^2
    # }
    # if (V_ftn == "binomial") {
    #   Veval <- function(x) x * (1-x)
    #   Vprime <- function(x) 1 - 2*x
    # }
    # 
    # factor2double <- function(x) as.double(as.character(x))
    # X <- mutate_if(X, is.factor, factor2double)
    # df_lag <- cbind(X, Y)
    # 
    # for (i_lag in 1:d_lag) {
    #   df_lag[["Y_" %++% i_lag]] <- dplyr::lag(df_lag[["Y"]], i_lag)  
    # }
    # df_lag <- na.omit(df_lag)
    # disc_str2 <- union(disc_str, "Y_" %++% 1:d_lag)
  }

  
  ### start optimization ###
  if (method == "NR") {
    ## Newton Raphson algorithm ##
    f_new <- matrix(init, ncol = 1)
    f_old <- f_new + NR_tolerance * 100
    i_iter <- 0
    while (norm(f_old - f_new, type = "F") > NR_tolerance) {
      i_iter <- i_iter + 1
      if (i_iter > NR_max_iter) break
      f_old <- f_new
      NR_cal <- MNQLE_NR_calculate(f_old, 
                                   x_est = point_est[1:length(conti_str)], 
                                   z_est = point_est[-(1:length(conti_str))], 
                                   x_bw = bw[1:length(conti_str)],
                                   z_bw = bw[-(1:length(conti_str))], 
                                   link_ftn, V_ftn, 
                                   conti_str, disc_str2, 
                                   "Y", 
                                   df_lag)
      if (abs(det(NR_cal$eq2)) <1e-05) {
        f_new <- f_old - solve(NR_cal$eq2 + diag(diag(NR_cal$eq2), 2) * .1) %*% NR_cal$eq1  
      } else {
        f_new <- f_old - solve(NR_cal$eq2) %*% NR_cal$eq1  
      }
    }
  } 
  return(ginv(f_new))
}

MNQLE_bootstrap <- function(X, Y, init = NULL, d_lag, point_boot, bw, link_ftn, V_ftn, NR_tolerance = 1e-05, NR_max_iter = 1000, method = "NR") {
  
}

MNQLE_genboot <- function(X, Y, init = NULL, d_lag, point_boot, bw, link_ftn, V_ftn, NR_tolerance = 1e-05, NR_max_iter = 1000, method = "NR") {
  if (FALSE) {
    X = df_analysis[c("x_1")]
    Y = df_analysis[, "y"]
    n <- nrow(df_analysis)
    point_est <- c(0.5, 1)
    bw = unname(c(n^{-1/5} * 1.06 * apply(X, 2, sd), n**-.4))
    
    link_ftn <- "probit"
    V_ftn <- "binomial"
    
    d_lag <- 1
    NR_tolerance = 1e-05
    NR_max_iter = 1000
    method = "GS"
    init = NULL
  }
  
  cmn_preproc(X, Y, 1, environment())
  cmn_defftns(link_ftn, V_ftn, environment())
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  
  
  MNQLE_optimize(X = df_analysis[conti_str], 
                 Y = df_analysis$y, 
                 init = NULL, d_lag, 
                 point_est = c(df_simul[ii, "x_1"], df_simul[ii, "z_1"]), 
                 bw = c(x_bw_init, z_bw_init), 
                 link_ftn = "probit", V_ftn = "binomial", method = "NR")[1]
  
}

QLNE_maximize <- function(x_est, z_est, x_bw, z_bw, link_ftn, V_ftn, conti_str, disc_str, d_lag, y_str, df_analysis) {
  d <- length(conti_str)
  f_init <- rnorm(1 + d)
  if (link_ftn == "probit") {
    ginv <- pnorm
    gprime <- function(x) dnorm(qnorm(x))
  }
  
  if (length(x_est) != d) stop("dimension!")
  if (length(z_est) != length(disc_str) + d_lag) stop("dimension!")
  
  for (i_lag in 1:d_lag) {
    df_analysis[["y_" %++% i_lag]] <- dplyr::lag(df_analysis[[y_str]], i_lag)  
  }
  df_analysis <- na.omit(df_analysis)
  n <- nrow(df_analysis)
  disc_str2 <- union(disc_str, "y_" %++% 1:d_lag)
  
  f_grid <- expand.grid(f_0 = seq(-2.5, 1, length.out = 100), f_1 = seq(-2.5, 1, length.out = 100))
  for (ii in 1:nrow(f_grid)) {
    f_eval <- f_grid[ii, 1:(1+d)] %>% as.numeric
    
    mi <- ginv(data.matrix(sweep(df_analysis[conti_str], 2, x_est)) %*% 
                 matrix(f_eval[-1], ncol = 1) + 
                 f_eval[1])[, 1]
    #if (ii %% 87 == 0) cat("grid search at", ii, mean(mi), '\n')
    f_grid$obj[[ii]] <- mean(
      x_kern(x_bw, x_est, df_analysis[[conti_str]]) * 
        #z_kern(z_bw, z_est, df_analysis[[disc_str2]]) * 
        (df_analysis[[y_str]] * log(mi / (1-mi)) + log(1-mi))
    )
  }
  
  opt_ii <- which.max(f_grid$obj)
  return(as.numeric(f_grid[opt_ii, 1:(1+d)]))
}

# 2. generate data -----------------------------------------------------------
set.seed(111)
## 2-1. Example 1 ----------------------------------------------------------
n <- 501
y_init <- 0
x_bd <- c(-2, 2)
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
### old version
if (FALSE) {
  ## input : determined from data
  conti_str <- c("x_1")
  disc_str <- character(0)
  d_lag <- 1
  y_str <- c("y")
  
  ## inpt : define now
  x_est <- 0.5
  z_est <- 1
  x_bw_init <- n^{-1/5} * 1.06 * apply(df_analysis[conti_str], 2, sd)
  z_bw_init <- n**-.4
  link_ftn <- "probit"
  V_ftn <- function(x) x*(1-x)
  df_simul2 <- as.data.table(df_simul)
  
  ## running
  list_result <- list()
  df_bw <- expand.grid(x_bw = x_bw_init+c(-.2, -.1, 0, .1, .2), z_bw = c(z_bw_init, z_bw_init + .2))
  for (ii in 1:nrow(df_bw)) {
    x_bw <- df_bw[ii, 1]
    z_bw <- df_bw[ii, 2]
    cat('bandwidth (x, z)', x_bw, z_bw, ": starts at", as.character(Sys.time()), "\n")
    df_tmp <- df_simul
    df_tmp$QLNE <- df_tmp$QLNE2 <- NA
    for (i in 1:nrow(df_simul)) {
      if (i %% 100 == 0) cat("iterating", i, ": starts at", as.character(Sys.time()), "\n")
      df_tmp$QLNE[i] <- QLNE(df_simul[i,1], df_simul[i,2], x_bw, z_bw, link_ftn, V_ftn, conti_str, disc_str, d_lag, y_str, df_analysis)[1] %>% pnorm
      df_tmp$QLNE2[i] <- QLNE_maximize(df_simul[i,1], df_simul[i,2], x_bw, z_bw, link_ftn, V_ftn, conti_str, disc_str, d_lag, y_str, df_analysis)[1] %>% pnorm
    }
    list_result[[ii]] <- df_tmp
  }
}

### up-to-date version ###
## input : determined from data
conti_str <- c("x_1")
disc_str <- character(0)
d_lag <- 1
y_str <- c("y")

## inpt : define now
x_est <- 0.5
z_est <- 1
x_bw_init <- n^{-1/5} * 1.06 * apply(df_analysis[conti_str], 2, sd)
z_bw_init <- n**-.4
link_ftn <- "probit"
V_ftn = "binomial"
df_simul$est <- NA
for (ii in 2:nrow(df_simul)) {
  df_simul[ii, "est"] <- 
    MNQLE_optimize(X = df_analysis[conti_str], 
                   Y = df_analysis$y, 
                   init = NULL, d_lag, 
                   point_est = c(df_simul[ii, "x_1"], df_simul[ii, "z_1"]), 
                   bw = c(x_bw_init, z_bw_init), 
                   link_ftn = "probit", V_ftn = "binomial", method = "NR")[1]
}

par(mfrow = c(1,2))
plot(df_simul$x_1,
     df_simul$mean, 
     col = df_simul$z_1 + 1)
plot(df_simul$x_1,
     df_simul$est, 
     col = df_simul$z_1 + 1)
par(mfrow = c(1, 1))
plot(df_simul$mean - df_simul$est)
summary(df_simul$mean - df_simul$est)



# 4. bootstrap -----------------------------------------------------------

