### final update : 20190321
### considerations : 
### - add function MNQLE_ABoot_step and MNQLE_onestep_cal
### - no discrete covariates except lagged variables
### - 1 lagged variable
### update 20190321
### - save MNQLE_tr for RB after AB

# 1. settings ----------------------------------------------------------------
# anticipated settings:
# 2 function
# 1 alpha : 95%
# 4 points : -1,1 -1,0 1.5,1 1.5,0 for simul1
# 1~2 sample size : n = 500,
# 3 methods : asymp / RB / AB
# 5 bandwidthd : undersmoothing with log(n)^-c(.1, .2, .3, .4, .5)
## library
library(dplyr)

source("Rsource/MNQLE_ftns.R")

# Rscript Rsource/MNQLE_simul.R simul1 501 -1 0 .1 .05
simul_inputs <- commandArgs(TRUE)
# simul_inputs <- c('simul1', '101', '-1', '1', '.1', '.05')
reg_ftn_str <- simul_inputs[1]
if (reg_ftn_str == "simul1") reg_ftn <- f2
if (reg_ftn_str == "simul2") reg_ftn <- f3
n_sample <- as.integer(simul_inputs[2])
point_est <- as.numeric(simul_inputs[3:4])
bw_factor <- log(n_sample)^(-as.numeric(simul_inputs[5]))
CI_alpha <- as.numeric(simul_inputs[6])

cat("simul_inputs =", simul_inputs, '\n')
cat("starting time = ", as.character(Sys.time()), '\n')

#2. bootstrap -----------------------------------------------------------
n_rep <- 500
n_boot <- 500

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

save_tr <- list()
save_rep <- matrix(0, nrow = n_rep, ncol = 3, dimnames = list(NULL, c("compare", "CI_L", "CI_U")))
pb <- txtProgressBar(min = 0, max = n_rep, style = 3, width = 50) # verbose
for (i_rep in 1:n_rep) {
  # MC data
  df_analysis <- MC_gen(reg_ftn = reg_ftn, n = n_sample * 10, y_init = 0, x_bd = x_bdry, list_ftns = cmn_defftns("probit", "binomial")) %>%
    tail(n_sample)
  # set bandwidths
  x_bw <- n_sample^{-1/5} * 1.06 * apply(df_analysis[conti_str], 2, sd) * bw_factor
  z_bw <- n_sample**-.4 * bw_factor
  
  # generate bootstrap estimators
  MNQLE_tr <- MNQLE_train(X = df_analysis[conti_str], Y = df_analysis$y, 
                          d_lag = d_lag, bw = c(x_bw, z_bw), verbose = FALSE)
  boot_AR <- MNQLE_ABoot_step(n_boot = n_boot, train_model = MNQLE_tr, point = point_est, verbose = FALSE)
  
  # compare CI and calculate coverage
  compare_val <- attr(boot_AR, "f_hat")[1] - reg_ftn(point_est[1], point_est[2])
  CI_boot <- quantile(boot_AR[, 1], c(CI_alpha / 2, 1 - CI_alpha / 2))

  save_tr[[i_rep]] <- MNQLE_tr
  save_rep[i_rep, 1] <- compare_val
  save_rep[i_rep, 2:3] <- CI_boot
  
  setTxtProgressBar(pb, i_rep) # verbose
}

# save Rdata
save.image("outdata/simul_" %++% paste0(simul_inputs, collapse = "_") %++% ".Rdata")

boot_coverage <- mean(save_rep[, "compare"] >= save_rep[, "CI_L"] & save_rep[, "compare"] <= save_rep[, "CI_U"])
CI_avg_length <- mean(save_rep[, "CI_U"] - save_rep[, "CI_L"])

# save csv
df_out <- matrix(c("AB", simul_inputs, boot_coverage, CI_avg_length), nrow = 1) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("method", "setting", "n", "point_x", "point_z", "bw_log_exponent", "CI_alpha", "empirical_coverage", "CI_avg_length"))
if (!file.exists("output/simulation_d1_k1.csv")) {
    write.table(df_out, "output/simulation_d1_k1.csv", sep = ",", row.names = FALSE)
} else {
    write.table(df_out, "output/simulation_d1_k1.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
}
