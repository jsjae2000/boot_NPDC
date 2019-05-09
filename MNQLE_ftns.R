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

f2 <- function(x, y_1) {
  beta0 <- -.2
  beta1 <- -1.75
  beta2 <- 2
  
  out <- beta0 + sin(beta1*x + beta2*y_1)
  out <- unname(out)
  return(out)
}

f3<- function(x, y_1) {
  beta0 <- .3
  beta1 <- -.5
  beta2 <- .2
  gamma <- -.2
  
  out <- beta0 + beta1*x + beta2*y_1 + gamma*x^2
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
  # Z : data frame s.t. ncol(Z) == d with nrow(Z) observations
  
  # ### toy example
  # Z <- data.frame(c(1,1,1,2,2,2,3,3,3), c(1,1,1,2,2,2,3,3,3))
  # z <- c(3,1)
  # bw <- c(.5, .3)
  
  # admit the case of common bandwidths
  if (length(bw) == 1) bw <- rep(bw, ncol(Z))
  
  # modify for discrete kernel type
  count_categ <- apply(Z, 2, function(j_col) (length(unique(j_col)))) %>% unname
  eq_weight <- 1 - bw
  neq_weight <- bw / (count_categ - 1)
  
  return(apply(Z, 1, function(Z_vec) prod(c(eq_weight[z == Z_vec], neq_weight[z != Z_vec]))))
  
  # previous discrete kernel
  if (FALSE) {
    apply(Z, 1, function(Z_vec) prod(bw[z != Z_vec]))
  }
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
    df_lag_start <- df_lag[1:d_lag, ]
    df_lag <- df_lag[-(1:d_lag), ]
    disc_str <- union(disc_str, "resp__lagged__" %++% 1:d_lag)
  }
  
  attr(df_lag, "conti_str") <- conti_str
  attr(df_lag, "disc_str") <- disc_str
  attr(df_lag, "head") <- df_lag_start
  
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

MNQLE_est <- function(X, Y, d_lag, point, bw, link_ftn = "probit", V_ftn = "binomial", method = "NR", init = NULL, tolerance = 1e-05, max_iter = 200) {
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
  x_pt <- point[i_conti]
  z_pt <- point[i_disc]
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  return(MNQLE_optimize(x_pt, z_pt, x_bw, z_bw, ftns, df_lag, method, init, tolerance, max_iter))
}

MNQLE_train <- function(X, Y, d_lag, bw, link_ftn = "probit", V_ftn = "binomial", method = "NR", init = NULL, tolerance = 1e-05, max_iter = 200, verbose = TRUE) {
  if (FALSE) {
    X = df_analysis[c("x_1")]
    Y = df_analysis[, "y"]
    d_lag <- 1
    
    point <- c(0.5, 1)
    bw = unname(c(n_sample^{-1/5} * 1.06 * apply(X, 2, sd), n_sample**-.4))
    
    link_ftn <- "probit"
    V_ftn <- "binomial"
    
    method = "NR"
    init = NULL
    tolerance = 1e-05
    max_iter = 1000
    verbose = TRUE
  }
  ### TO execute autoregressive bootstrap,
  ### generate df_train object which consits of columns all continuous 
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(link_ftn, V_ftn)
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  
  i_conti <- match(attr(df_lag, "conti_str"), names(df_lag)[-1])
  i_disc <- match(attr(df_lag, "disc_str"), names(df_lag)[-1])
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  if (d_lag > 0) {
    lagresp <- lapply(df_lag["resp__lagged__" %++% 1:d_lag], unique) %>%
      expand.grid
    
    
    if (verbose) {
      i_verbose <- 0
      envhere <- environment()
      pb <- txtProgressBar(min = 0, max = nrow(df_lag) * nrow(lagresp), style = 3, width = 50)
    }
    
    train_model <- lapply(1:nrow(lagresp), 
                          function(i_lagresp) {
                            lagresp_value <- unlist(lagresp[i_lagresp, ])
                            apply(df_lag, 
                                  1,
                                  function(x_value) {
                                    est_value <- MNQLE_optimize(x_value[attr(df_lag, "conti_str")],
                                                                lagresp_value, 
                                                                x_bw, z_bw, ftns, df_lag, method, init, tolerance, max_iter)
                                    if (verbose) {
                                      assign("i_verbose", i_verbose + 1, envhere)
                                      setTxtProgressBar(pb, i_verbose)
                                    }
                                    return(c(mean = est_value$mean, est_value$reg))
                                  }) %>%
                              t() %>%
                              as.data.frame()
                          })
    
    attr(train_model, "df_lagresp") <- lagresp
  } else {
    stop("under construction")
  }
  
  attr(train_model, "model") <- list(X = X, Y = Y, d_lag = d_lag, bw = bw, link_ftn = link_ftn, V_ftn = V_ftn,
                                     method = method, init = init, tolerance = tolerance, max_iter = max_iter)
  
  return(train_model)
}

MNQLE_ABoot <- function(n_boot = 500, train_model, point, bw, link_ftn = "probit", V_ftn = "binomial", 
                        method = "NR", init = NULL, tolerance = 1e-05, max_iter = 1000, save_sample = TRUE, verbose = TRUE) {
  if (FALSE) {
    n_boot <- 500
    train_model <- MNQLE_tr
    
    d_lag <- attr(train_model, "model")$d_lag
    
    point <- c(0.5, 1)
    bw = unname(c(n^{-1/5} * 1.06 * apply(X, 2, sd), n**-.4))
    
    link_ftn <- "probit"
    V_ftn <- "binomial"
    
    method = "NR"
    init = NULL
    tolerance = 1e-05
    max_iter = 1000
  }
  X = attr(train_model, "model")$X
  Y = attr(train_model, "model")$Y
  d_lag <- attr(train_model, "model")$d_lag
  if (d_lag == 0) stop("use function for regression bootstrap")
  
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(link_ftn, V_ftn)
  if (is.null(init)) init <- rnorm(1 + length(attr(df_lag, "conti_str")))
  
  i_conti <- match(attr(df_lag, "conti_str"), names(df_lag)[-1])
  i_disc <- match(attr(df_lag, "disc_str"), names(df_lag)[-1])
  x_pt <- point[i_conti]
  z_pt <- point[i_disc]
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  ### bootstrap ###
  # calculate necessary index to calculate bootstrap samples to estimate at point
  df_conti <- df_lag[attr(df_lag, "conti_str")]
  i_calboot <- 1
  for (i in 1:ncol(df_conti)) {
    i_calboot <- max(i_calboot, max(which((df_conti[, i] <= x_est[i] + x_bw[i]) & (df_conti[, i] >= x_est[i] - x_bw[i]))))
  }
  
  # prepare to save bootstrap samples
  if (save_sample) boot_sample <- matrix(NA_integer_, nrow = i_calboot + d_lag, ncol = n_boot)
  if (verbose) pb <- txtProgressBar(min = 0, max = n_boot, style = 3, width = 50)
  
  est_boot <- c()
  for (i_boot in 1:n_boot) {
    Y_boot <- attr(df_lag, "head")[["resp__var__"]]
    for (i_df in 1:i_calboot) {
      # generate bootstrap sample
      lagresp_value <- rev(tail(Y_boot, d_lag))
      i_lagresp <- which(apply(attr(train_model, "df_lagresp"), 1, function(x) x == lagresp_value))
      Y_boot <- c(Y_boot, rbinom(1, 1, train_model[i_df, i_lagresp]))
    }
    
    if (save_sample) boot_sample[, i_boot] <- Y_boot
    
    # calculate bootstrap estimates
    df_lag_boot <- cmn_preproc(slice(X, 1:(d_lag + i_calboot)), Y_boot, d_lag)
    est_boot[i_boot] <- MNQLE_optimize(x_pt, z_pt, x_bw, z_bw, ftns, df_lag_boot, method, init, tolerance, max_iter)
    
    # print if verbose == TURE
    if (verbose) setTxtProgressBar(pb, i_boot)
  }
  
  if (save_sample) attr(est_boot, "samples") <- boot_sample
  return(est_boot)
}

##### function to generate Autoregression Bootstrap samples and estimators using one-step approximation
MNQLE_ABoot_step <- function(n_boot = 500, train_model, point, save_sample = TRUE, verbose = TRUE) {
  if (FALSE) {
    point <- point_est
    n_boot <- 500
    train_model <- MNQLE_tr
    save_sample = TRUE
    verbose = TRUE
  }
  
  ### defining and preprocessing ###
  X <- attr(train_model, "model")$X
  Y <- attr(train_model, "model")$Y
  d_lag <- attr(train_model, "model")$d_lag
  if (d_lag == 0) stop("use function for regression bootstrap")
  
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(attr(train_model, "model")$link_ftn, attr(train_model, "model")$V_ftn)
  
  i_conti <- match(attr(df_lag, "conti_str"), names(df_lag)[-1])
  i_disc <- match(attr(df_lag, "disc_str"), names(df_lag)[-1])
  x_pt <- point[i_conti]
  z_pt <- point[i_disc]
  
  bw <- attr(train_model, "model")$bw
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  ### estimation at point ###
  model_ <- attr(train_model, "model")
  model_$point <- point
  f_hat <- do.call(MNQLE_est, model_)$reg
  
  ### estimateion of D ###
  list_Dhat_cal <- MNQLE_onestep_cal1(f_hat, x_pt, z_pt, x_bw, z_bw, ftns, df_lag)
  
  ### bootstrap ###
  # calculate necessary index to calculate bootstrap samples to estimate at point
  df_conti <- df_lag[attr(df_lag, "conti_str")]
  i_calboot <- 1
  for (i in 1:ncol(df_conti)) {
    i_calboot <- max(i_calboot, max(which((df_conti[, i] <= x_pt[i] + x_bw[i]) & (df_conti[, i] >= x_pt[i] - x_bw[i]))))
  }
  
  # prepare to save bootstrap samples
  if (save_sample) {
    boot_sample <- matrix(NA_integer_, nrow = i_calboot + d_lag, ncol = n_boot)
    mat_boot_mean <- matrix(NA_integer_, nrow = i_calboot, ncol = n_boot) # mean used to generate boot_sample
  }
  if (verbose) pb <- txtProgressBar(min = 0, max = n_boot, style = 3, width = 50)
  
  est_boot <- matrix(NA_real_, nrow = n_boot, ncol = 1 + d_lag)
  for (i_boot in 1:n_boot) {
    Y_boot <- attr(df_lag, "head")[["resp__var__"]] # length  = n + d_lag
    Y_boot_mean <- c() # vector to save m hat (Xi, ZI*); length = n
    lagresp_value <- Y_boot
    for (i_df in 1:i_calboot) {
      # generate bootstrap sample
      #lagresp_value <- rev(tail(Y_boot, d_lag))
      i_lagresp <- which(apply(attr(train_model, "df_lagresp"), 1, function(x) x == lagresp_value))
      mhat <- train_model[[i_lagresp]][i_df, 1]
      Y_boot_1 <- rbinom(1, 1, mhat)
      lagresp_value <- c(lagresp_value[-1], Y_boot_1)
      Y_boot <- c(Y_boot, Y_boot_1)
      Y_boot_mean <- c(Y_boot_mean, mhat)
    }
    
    if (save_sample) {
      boot_sample[, i_boot] <- Y_boot
      mat_boot_mean[, i_boot] <- Y_boot_mean
    }
    
    # calculate bootstrap estimates
    df_lag_boot <- cmn_preproc(X[1:(i_calboot + d_lag),,drop = FALSE], Y_boot, d_lag)
    F_hat_boot <- MNQLE_onestep_cal2(z_pt, z_bw, ftns, list_Dhat_cal, df_lag_boot, Y_boot_mean)
    est_boot[i_boot, ] <- solve(list_Dhat_cal$nD_hat) %*% F_hat_boot
    
    # print if verbose == TURE
    if (verbose) setTxtProgressBar(pb, i_boot)
  }
  
  if (save_sample) {
    attr(est_boot, "samples") <- boot_sample
    attr(est_boot, "mat_mean") <- mat_boot_mean
  }
  attr(est_boot, "f_hat") <- f_hat
  return(est_boot)
}

MNQLE_RBoot_step <- function(n_boot = 500, train_model, point, save_sample = TRUE, verbose = TRUE) {
  if (FALSE) {
    point <- c(-1.5, 0)
    n_boot <- 500
    train_model <- MNQLE_tr
    save_sample = TRUE
    verbose = TRUE
  }
  
  ### defining and preprocessing ###
  X <- attr(train_model, "model")$X
  Y <- attr(train_model, "model")$Y
  d_lag <- attr(train_model, "model")$d_lag
  if (d_lag == 0) stop("use function for regression bootstrap")
  
  df_lag <- cmn_preproc(X, Y, d_lag)
  ftns <- cmn_defftns(attr(train_model, "model")$link_ftn, attr(train_model, "model")$V_ftn)
  
  i_conti <- match(attr(df_lag, "conti_str"), names(df_lag)[-1])
  i_disc <- match(attr(df_lag, "disc_str"), names(df_lag)[-1])
  x_pt <- point[i_conti]
  z_pt <- point[i_disc]
  
  bw <- attr(train_model, "model")$bw
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  ### estimation at point ###
  model_ <- attr(train_model, "model")
  model_$point <- point
  f_hat <- do.call(MNQLE_est, model_)$reg
  
  ### estimateion of D ###
  list_Dhat_cal <- MNQLE_onestep_cal1(f_hat, x_pt, z_pt, x_bw, z_bw, ftns, df_lag)
  
  ### modify train_model from AB to adapt to RB
  if (length(train_model) > 1) {
    train_RB <- matrix(nrow = 0, ncol = ncol(X) + d_lag)
    for (i_lagresp in 1:nrow(attr(train_model, "df_lagresp"))) {
      train_RB <- cbind(tail(attr(train_model, "model")$X, -d_lag), 
                        attr(train_model, "df_lagresp")[i_lagresp, , drop = FALSE], 
                        meanhat = train_model[[i_lagresp]]$mean) %>%
        rbind(train_RB, .)
    }
    
    train_RB <- left_join(df_lag, train_RB, by = colnames(train_RB)[-ncol(train_RB)])
    train_RB <- train_RB[, ncol(train_RB)]
  }
  
  ### bootstrap ###
  # calculate necessary index to calculate bootstrap samples to estimate at point
  df_conti <- df_lag[attr(df_lag, "conti_str")]
  i_calboot <- 1
  for (i in 1:ncol(df_conti)) {
    i_calboot <- max(i_calboot, max(which((df_conti[, i] <= x_pt[i] + x_bw[i]) & (df_conti[, i] >= x_pt[i] - x_bw[i]))))
  }
  
  # prepare to save bootstrap samples
  if (save_sample) {
    boot_sample <- matrix(NA_integer_, nrow = i_calboot + d_lag, ncol = n_boot)
    mat_boot_mean <- matrix(NA_integer_, nrow = i_calboot, ncol = n_boot) # mean used to generate boot_sample
  }
  if (verbose) pb <- txtProgressBar(min = 0, max = n_boot, style = 3, width = 50)
  
  est_boot <- matrix(NA_real_, nrow = n_boot, ncol = 1 + d_lag)
  for (i_boot in 1:n_boot) {
    Y_boot <- attr(df_lag, "head")[["resp__var__"]] # length  = n + d_lag
    Y_boot_mean <- c()
    for (i_df in 1:i_calboot) {
      # generate bootstrap sample
      #lagresp_value <- rev(tail(Y_boot, d_lag))
      mhat <- train_RB[i_df]
      Y_boot_1 <- rbinom(1, 1, mhat)
      Y_boot <- c(Y_boot, Y_boot_1)
      Y_boot_mean <- c(Y_boot_mean, mhat)
    }
    
    if (save_sample) {
      boot_sample[, i_boot] <- Y_boot
      mat_boot_mean[, i_boot] <- Y_boot_mean
    }
    
    # calculate bootstrap estimates
    df_lag_boot <- cmn_preproc(X[1:(i_calboot + d_lag),,drop = FALSE], Y_boot, d_lag)
    F_hat_boot <- MNQLE_onestep_cal2(z_pt, z_bw, ftns, list_Dhat_cal, df_lag_boot, Y_boot_mean)
    est_boot[i_boot, ] <- solve(list_Dhat_cal$nD_hat) %*% F_hat_boot
    
    # print if verbose == TURE
    if (verbose) setTxtProgressBar(pb, i_boot)
  }
  
  if (save_sample) {
    attr(est_boot, "samples") <- boot_sample
    attr(est_boot, "mat_mean") <- mat_boot_mean
  }
  attr(est_boot, "f_hat") <- f_hat
  return(est_boot)
}

MNQLE_RBoot <- function(n_boot = 500, X, Y, estimates, d_lag, point, bw, link_ftn = "probit", V_ftn = "binomial", 
                        method = "NR", init = NULL, tolerance = 5e-04, max_iter = 1000, save_sample = TRUE, verbose = TRUE) {
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
  x_pt <- point[i_conti]
  z_pt <- point[i_disc]
  x_bw <- bw[i_conti]
  z_bw <- bw[i_disc]
  
  if (nrow(df_lag) != length(estimates)) stop("check dimension of X and estimates")
  
  
  ### bootstrap ###
  X_boot <- subset(df_lag, select =  !(names(df_lag) %in% c("resp__var__", "resp__lagged__" %++% 1:d_lag)))
  
  # prepare to save bootstrap samples
  if (save_sample) boot_sample <- matrix(NA_integer_, nrow = length(estimates), ncol = n_boot)
  if (verbose) pb <- txtProgressBar(min = 0, max = n_boot, style = 3, width = 50)
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
    est_boot[i_boot] <- MNQLE_optimize(x_pt, z_pt, x_bw, z_bw, ftns, df_lag_boot, method, init, tolerance, max_iter)
    
    # print if verbose == TURE
    if (verbose) setTxtProgressBar(pb, i_boot)
  }
  
  if (save_sample) attr(est_boot, "samples") <- boot_sample
  return(est_boot)
}

MNQLE_optimize <- function(x_pt, z_pt, x_bw, z_bw, ftns, df_lag, method = "NR", init = NULL, tolerance = 1e-05, max_iter = 200) {
  conti_str <- attr(df_lag, "conti_str")
  disc_str <- attr(df_lag, "disc_str")
  if (is.null(init)) init <- rnorm(1 + length(conti_str))
  
  ### start optimization ###
  ## Newton Raphson algorithm ##
  f_new <- matrix(init, ncol = 1)
  f_old <- f_new + tolerance * 100
  i_iter <- 0
  
  while (norm(f_old - f_new, type = "F") > tolerance) {
    i_iter <- i_iter + 1
    if (i_iter > max_iter) {
      # if f_new does not converge, call MNQLE_optimize_search
      f_new <- MNQLE_optimize_search(x_pt, z_pt, x_bw, z_bw, ftns, df_lag, tolerance * 1e02)
      f_new <- matrix(f_new, ncol = 1)
      rownames(f_new) <- c("intercept", conti_str)
      break
    }
    f_old <- f_new
    NR_cal <- MNQLE_loglik_calcul(f_old, x_pt, z_pt, x_bw, z_bw, ftns, df_lag)
    
    cal_inv <- tryCatch({
      f_new <- f_old - solve(NR_cal$eq2) %*% NR_cal$eq1 # unique difference with NR
    }, error = function(e) e)
    
    if (inherits(cal_inv, "error") | any(is.na(f_new))) {
      # if NR_cal$eq2 is not invertible, call MNQLE_optimize_search
      f_new <- MNQLE_optimize_search(x_pt, z_pt, x_bw, z_bw, ftns, df_lag, tolerance * 1e02)
      f_new <- matrix(f_new, ncol = 1)
      rownames(f_new) <- c("intercept", conti_str)
      break
    }
    
    # old version
    # modified to NR and NR_fisher (20190422)
    if (FALSE) {
      if (abs(det(NR_cal$eq2)) <1e-10) {
        print("ridge")
        f_new <- f_old - solve(NR_cal$eq2 + diag(diag(NR_cal$eq2), 2) * .1) %*% NR_cal$eq1  
      } else {
        f_new <- f_old - solve(NR_cal$eq2) %*% NR_cal$eq1  
      }
    }
  } 
  
  if (method == "FALSE") {
    while (norm(f_old - f_new, type = "F") > tolerance) {
      i_iter <- i_iter + 1
      if (i_iter > max_iter) {
        print("maxiter with " %++% tolerance)
        #MNQLE_optimize(x_pt, z_pt, x_bw, z_bw, ftns, df_lag, method, init = NULL, tolerance = tolerance * 5, max_iter)
        break
      }
      f_old <- f_new
      NR_cal <- MNQLE_loglik_calcul(f_old, x_pt, z_pt, x_bw, z_bw, ftns, df_lag)
      
      f_new <- f_old - solve(NR_cal$eq2_infor) %*% NR_cal$eq1 # unique difference with NR
    } 
  }
  
  return(list(mean = ftns$ginv(f_new[1]), reg = f_new[, 1]))
}

MNQLE_optimize_search <- function(x_pt, z_pt, x_bw, z_bw, ftns, df_lag, tolerance, n_search) {
  # create date : 20190503
  # goal : find MNQLE using coordinatewise grid search(?)
  # note : if it appears error in MNQLE_optimzie, then call this function.
  # note : 
  # - "n_search" is the number of candidates for search
  # - "monitor_conv" is the 3 consequtive difference between f_old and f_new
  # - "int_search" regulates the length of intervals for grid in search. 
  #   It is determined by mean(monitor_conv)
  n_search <- 20
  max_search <- 100
  
  search_fun <- function(f_old, i_search, val_search) {
    f_old[i_search] <- val_search
    MNQLE_score_calcul(f_old, x_pt, z_pt, x_bw, z_bw, ftns, df_lag)
  }
  
  
  f_new <- rnorm(length(attr(df_lag, "conti_str")) * 2)
  f_dim <- length(f_new)
  monitor_conv <- rep(tolerance * 10, 3)
  i_iter <- 1
  i_search <- 1
  tryCatch({
    while (mean(monitor_conv) > tolerance) {
      if (i_iter > max_search) break
      f_old <- f_new
      
      int_search <- min(mean(monitor_conv) * 50, 5)
      for (i_search in 1:length(f_old)) {
        cand_search <- runif(n_search, f_old[i_search] - int_search, f_old[i_search] + int_search)
        cand_search <- sort(c(cand_search, f_old[i_search]))
        val_score <- sapply(cand_search, search_fun, f_old = f_old, i_search = i_search) %>% t()
        f_new[i_search] <- cand_search[which.min(apply(val_score, 1, function(x) sum(x**2)**.5))]
      }
      
      monitor_conv <- c(tail(monitor_conv, -1), norm(matrix(f_old - f_new), type = "F"))
      i_iter <- i_iter + 1
      i_search <- (i_search %% length(f_old)) + 1
    }
  }, error = function(e) {
    cat("error", x_pt, x_bw, "/", int_search, "/", cand_search, "/", val_score, "/", f_new, "/", f_old)
  })

  return(f_new)
}

MNQLE_score_calcul <- function(f_at, x_pt, z_pt, x_bw, z_bw, ftns, df_lag) {
  # create date : 20190503
  # goal : calculate score value at f_at
  # note : same with the first part of MNQLE_loglik_calcul
  
  conti_str <- attr(df_lag, "conti_str")
  disc_str <- attr(df_lag, "disc_str")
  
  wci_all <- x_kern(x_bw, x_pt, df_lag[conti_str])
  wdi_all <- z_kern(z_bw, z_pt, df_lag[disc_str])
  cal_ind <- wci_all != 0
  wci <- wci_all[cal_ind]
  wdi <- wdi_all[cal_ind]
  df_cal <- df_lag[cal_ind, ]
  if (nrow(df_cal) ==0) return(matrix(0, nrow = 1 + length(conti_str)))
  
  mi <- ftns$ginv(data.matrix(sweep(df_cal[conti_str], 2, x_pt)) %*%
                    matrix(f_at[-1], ncol = 1) +
                    f_at[1]
  )[, 1]
  
  if (TRUE) {
    mi[mi >= 1 - 1e-05] <- 1 - 1e-05
    mi[mi <= 1e-05] <- 1e-05
  }
  ### calculate eq1 val ###
  q1i <- (df_cal[, 1] - mi) / (ftns$Veval(mi) * ftns$gprime(mi))
  product1_i <- matrix(wci * wdi * q1i, nrow = 1)
  # dim(Wi) = n X d
  Wi <- cbind(intercept = 1, sweep(df_cal[conti_str], 2, x_pt) %>% sweep(2, FUN = "/", x_bw)) %>%
    as.matrix()
  # eq1_val
  eq1_val <- (product1_i %*% Wi) / nrow(df_lag)
  
  return(eq1_val)
}

MNQLE_loglik_calcul <- function(f_at, x_pt, z_pt, x_bw, z_bw, ftns, df_lag) {
  if (FALSE) {
    f_at = f_old
    x_pt = point[1:length(conti_str)]
    z_pt = point[-(1:length(conti_str))] 
    x_bw = bw[1:length(conti_str)]
    z_bw = bw[-(1:length(conti_str))]
    ftns = ftns
    df_lag = df_lag
  }
  
  conti_str <- attr(df_lag, "conti_str")
  disc_str <- attr(df_lag, "disc_str")
  
  wci_all <- x_kern(x_bw, x_pt, df_lag[conti_str])
  wdi_all <- z_kern(z_bw, z_pt, df_lag[disc_str])
  cal_ind <- wci_all != 0
  wci <- wci_all[cal_ind]
  wdi <- wdi_all[cal_ind]
  df_cal <- df_lag[cal_ind, ]
  if (nrow(df_cal) ==0) return(matrix(0, nrow = 1 + length(conti_str)))
  
  mi <- ftns$ginv(data.matrix(sweep(df_cal[conti_str], 2, x_pt)) %*%
                    matrix(f_at[-1], ncol = 1) +
                    f_at[1]
  )[, 1]
  
  if (FALSE) {
    mi[mi >= 1 - 1e-05] <- 1 - 1e-05
    mi[mi <= 1e-05] <- 1e-05
  }
  ### calculate eq1 val ###
  q1i <- (df_cal[, 1] - mi) / (ftns$Veval(mi) * ftns$gprime(mi))
  product1_i <- matrix(wci * wdi * q1i, nrow = 1)
  # dim(Wi) = n X d
  Wi <- cbind(intercept = 1, sweep(df_cal[conti_str], 2, x_pt) %>% sweep(2, FUN = "/", x_bw)) %>%
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
  
  q2i_infor <- -q2i_tmp2
  product2_i_infor <- matrix(wci * wdi * q2i_infor, nrow = 1)
  eq2_val_infor <- (sweep(t(Wi), 2, product2_i_infor, FUN = "*") %*% Wi) / nrow(df_lag)
  
  return(list(eq1 = matrix(eq1_val, ncol = 1), eq2 = eq2_val, eq2_infor = eq2_val_infor))
}

MNQLE_onestep_cal1 <- function(f_hat, x_pt, z_pt, x_bw, z_bw, ftns, df_lag) {
  ### functions to calculate valeus used in MNQLE_ABoot_step
  conti_str <- attr(df_lag, "conti_str")
  disc_str <- attr(df_lag, "disc_str")
  
  wci_all <- x_kern(x_bw, x_pt, df_lag[conti_str])
  wdi_all <- z_kern(z_bw, z_pt, df_lag[disc_str])
  cal_ind <- wci_all != 0
  wci <- wci_all[cal_ind]
  wdi <- wdi_all[cal_ind]
  df_cal <- df_lag[cal_ind, ]
  if (nrow(df_cal) ==0) return(matrix(0, nrow = 1 + length(conti_str)))

  # From here, all objects are reduced using cal_ind
  mi_hat <- ftns$ginv(data.matrix(sweep(df_cal[conti_str], 2, x_pt)) %*%
                        matrix(f_hat[-1], ncol = 1) +
                        f_hat[1]
  )[, 1]
  mi_hat[mi_hat >= 1 - 1e-05] <- 1 - 1e-05
  mi_hat[mi_hat <= 1e-05] <- 1e-05
  
  ### calculate D_hat ###
  Wi <- cbind(intercept = 1, sweep(df_cal[conti_str], 2, x_pt) %>% sweep(2, FUN = "/", x_bw)) %>%
    as.matrix()
  pre_Wi <- wci * wdi / (ftns$Veval(mi_hat) * ftns$gprime(mi_hat)^2) # values to calculate D_hat
  nD_hat <- (sweep(t(Wi), 2, pre_Wi, FUN = "*") %*% Wi)  # n * D_hat
  
  return(list(nD_hat = nD_hat, wci = wci, mi_hat = mi_hat, Wi = Wi, cal_ind = cal_ind))
}

MNQLE_onestep_cal2 <- function(z_pt, z_bw, ftns, out_cal1, df_lag_boot, boot_mean) {
  if (FALSE) {
    tmp <- list_Dhat_cal
    nD_hat <- tmp$nD_hat
    wci <- tmp$wci
    mi_hat <- tmp$mi_hat
    Wi <- tmp$Wi
    cal_ind <- tmp$cal_ind
    boot_mean <- Y_boot_mean
  }
  conti_str <- attr(df_lag_boot, "conti_str")
  disc_str <- attr(df_lag_boot, "disc_str")
  
  wci <- out_cal1$wci
  mi_hat <- out_cal1$mi_hat
  Wi <- out_cal1$Wi
  cal_ind <- out_cal1$cal_ind
  
  # From here, all objects are reduced using cal_ind
  wdi <- z_kern(z_bw, z_pt, df_lag_boot[disc_str])[cal_ind]
  df_cal <- df_lag_boot[cal_ind, ]
  
  ### calculate F_hat_boot ###
  preWi <- wci * wdi * (df_lag_boot$resp__var__ - boot_mean)[cal_ind] / 
    (ftns$Veval(mi_hat) * ftns$gprime(mi_hat)) # values to calculate F_hat_boot
  nF_hat_boot <- t(Wi) %*% matrix(preWi, ncol = 1)
  
  return(nF_hat_boot)
}

plot.progress <- function(...)	{
  vectOfBar <- c(...)*100
  numOfBar <- length(vectOfBar)
  plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
  for(i in 1:numOfBar) {
    rect(0, 0.1+i-1, vectOfBar[i], 0.9+i-1, col=rainbow(numOfBar)[i])
    text(0.5, 0.5+i-1, paste('Status ', i, ': ', round(vectOfBar[i],2), '%', sep=''), adj=0)
  }
  title('Progress...')
}

MC_gen <- function(reg_ftn, n, y_init, x_bd, list_ftns) {
  # start generating
  y_vec <- c("0" = y_init)
  mean_vec <- double()
  x_vec = runif(n, min = x_bd[1], max = x_bd[2])
  
  for (i in 1:n) {
    mean_val <- reg_ftn(x_vec[i], 
                        y_vec[as.character(i-1)]) %>%
      (list_ftns$ginv)
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
  
  return(df_analysis)
}