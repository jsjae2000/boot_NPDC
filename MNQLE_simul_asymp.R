### final update : 20190414
### Goal :
### - simulation of asympt after AB(with saved model_tr objects)
### considerations : 
### - add function MNQLE_ABoot_step and MNQLE_onestep_cal
### - no discrete covariates except lagged variables
### - 1 lagged variable
### update 20190321
### - save MNQLE_tr for RB after AB

# 1. settings ----------------------------------------------------------------
# anticipated settings:
# - same with AB
## library
library(dplyr)
source("Rsource/MNQLE_ftns.R")

cand_load_rdata <- list.files("outdata", pattern = ".Rdata")
cand_load_rdata <- list.files("outdata", pattern = "simul_simul2_501_*")
for (load_rdata in cand_load_rdata) {
  if (FALSE) load_rdata <- "simul_simul1_501_-1_0_.1_.05.Rdata"
  load("outdata/" %++% load_rdata)
  
  cat(match(load_rdata, cand_load_rdata), "th file\n")
  cat("simul_inputs =", simul_inputs, '\n')
  cat("starting time = ", as.character(Sys.time()), '\n')
  
  #2. bootstrap -----------------------------------------------------------
  if (reg_ftn_str == "simul1") {
    conti_str <- c("x_1")
    disc_str <- character(0)
    d_lag <- 1
    x_bdry <- c(-3, 3)
  }
  if (reg_ftn_str == "simul2") {
    conti_str <- c("x_1")
    disc_str <- character(0)
    d_lag <- 1
    x_bdry <- c(-2, 2)
  }
  
  # modification for asymp
  ftrue <- reg_ftn(point_est[1], point_est[2])
  MC_fhat <- save_rep[, 1] + ftrue
  
  save_rep <- matrix(0, nrow = n_rep, ncol = 3, dimnames = list(NULL, c("compare", "CI_L", "CI_U")))
  pb <- txtProgressBar(min = 0, max = n_rep, style = 3, width = 50) # verbose
  for (i_rep in 1:n_rep) {
    MNQLE_tr <- save_tr[[i_rep]]

    if (FALSE) {
      # check
      MNQLE_est(attr(MNQLE_tr, "model")$X, attr(MNQLE_tr, "model")$Y, attr(MNQLE_tr, "model")$d_lag, point_est, attr(MNQLE_tr, "model")$bw)
    }
    
    # necessary values
    n_CI <- n_sample - attr(MNQLE_tr, "model")$d_lag
    bw_x <- attr(MNQLE_tr, "model")$bw[1]
    bw_z <- attr(MNQLE_tr, "model")$bw[2]
    myftns <- cmn_defftns(attr(MNQLE_tr, "model")$link_ftn, attr(MNQLE_tr, "model")$V_ftn)
    
    # calcualte fhat, mhat
    fhat <- MC_fhat[i_rep]
    mhat <- myftns$ginv(fhat)
    
    # calculate phat
    df_lag <- cmn_preproc(attr(MNQLE_tr, "model")$X, attr(MNQLE_tr, "model")$Y, attr(MNQLE_tr, "model")$d_lag)
    wci_all <- x_kern(bw_x, point_est[1], df_lag[attr(df_lag, "conti_str")])
    wdi_all <- z_kern(bw_z, point_est[2], df_lag[attr(df_lag, "disc_str")])
    phat <- mean(wci_all * wdi_all)
    
    # calculate CI
    # integral of K^2 is .6 for Epanechnikov kernels
    asymp_sd <- (n_CI * bw_x)^.5 * (myftns$gprime(mhat)^2 * mhat * (1-mhat) / phat)^-.5 * (.6)^-.5
    CI_asymp <- qnorm(c(CI_alpha / 2, 1 - CI_alpha / 2)) / asymp_sd + fhat
    
    # compare CI and calculate coverage
    save_rep[i_rep, 1] <- ftrue
    save_rep[i_rep, 2:3] <- CI_asymp
    
    setTxtProgressBar(pb, i_rep) # verbose
  }
  close(pb)
  
  # save Rdata
  #save.image("outdata/simul_" %++% paste0(simul_inputs, collapse = "_") %++% ".Rdata")
  
  boot_coverage <- mean(save_rep[, "compare"] >= save_rep[, "CI_L"] & save_rep[, "compare"] <= save_rep[, "CI_U"])
  CI_avg_length <- mean(save_rep[, "CI_U"] - save_rep[, "CI_L"])
  
  # save csv
  df_out <- matrix(c("Asymp", simul_inputs, boot_coverage, CI_avg_length), nrow = 1) %>%
    data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("method", "setting", "n", "point_x", "point_z", "bw_log_exponent", "CI_alpha", "empirical_coverage", "CI_avg_length"))
  if (!file.exists("output/simulation_d1_k1.csv")) {
    write.table(df_out, "output/simulation_d1_k1.csv", sep = ",", row.names = FALSE)
  } else {
    write.table(df_out, "output/simulation_d1_k1.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}
