### final update : 20190411
### Goal :
### - simulation of RB after AB(with saved model_tr objects)
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

cand_load_rdata <- list.files("outdata", pattern = ".Rdata")[3]
cand_load_rdata <- list.files("outdata", pattern = "simul_simul2_501_*")[10:18]
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
  
  save_rep <- matrix(0, nrow = n_rep, ncol = 3, dimnames = list(NULL, c("compare", "CI_L", "CI_U")))
  pb <- txtProgressBar(min = 0, max = n_rep, style = 3, width = 50) # verbose
  for (i_rep in 1:n_rep) {
    # load MNQLE_tr from AB
    MNQLE_tr <- save_tr[[i_rep]]
    boot_RB <- MNQLE_RBoot_step(n_boot = n_boot, train_model = MNQLE_tr, point = point_est, verbose = FALSE)
    
    # compare CI and calculate coverage
    compare_val <- attr(boot_RB, "f_hat")[1] - reg_ftn(point_est[1], point_est[2])
    CI_boot <- quantile(boot_RB[, 1], c(CI_alpha / 2, 1 - CI_alpha / 2))
    
    save_rep[i_rep, 1] <- compare_val
    save_rep[i_rep, 2:3] <- CI_boot
    
    setTxtProgressBar(pb, i_rep) # verbose
  }
  close(pb)
  
  # save Rdata
  #save.image("outdata/simul_" %++% paste0(simul_inputs, collapse = "_") %++% ".Rdata")
  
  boot_coverage <- mean(save_rep[, "compare"] >= save_rep[, "CI_L"] & save_rep[, "compare"] <= save_rep[, "CI_U"])
  CI_avg_length <- mean(save_rep[, "CI_U"] - save_rep[, "CI_L"])
  
  # save csv
  df_out <- matrix(c("RB", simul_inputs, boot_coverage, CI_avg_length), nrow = 1) %>%
    data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("method", "setting", "n", "point_x", "point_z", "bw_log_exponent", "CI_alpha", "empirical_coverage", "CI_avg_length"))
  if (!file.exists("output/simulation_d1_k1.csv")) {
    write.table(df_out, "output/simulation_d1_k1.csv", sep = ",", row.names = FALSE)
  } else {
    write.table(df_out, "output/simulation_d1_k1.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}
