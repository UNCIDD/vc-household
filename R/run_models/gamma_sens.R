library(rstan)
library(tidyverse)
source("R/helper_scripts/cholera_utilities.R")

gammas <- c(0.5, 0.8, 1.2, 2)
gamma_list <- list()

for(i in 1:length(gammas)) {
  
  gamma <- gammas[i]
  dat <- readRDS(paste0("data/dat_11.26.24_nocov_320.rds"))
  
  dat$gamma_mult <- 1/gamma
  
  # res_gamma <- stan(file = "stan/HMM_cholera_nocov.stan",
  #                 data = dat,
  #                 cores = cores,
  #                 pars = c("beta_asym", "beta_ih", "beta_eh",
  #                          "eh_prob", "ih_prob_sym", "ih_prob_asym",
  #                          "gamma", "sigma", "p_symp", "init_probs_hh"),
  #                 iter = 2000,
  #                 init = rep(list(list(beta_asym = 0,
  #                                      beta_ih = logit(0.02),
  #                                      beta_eh = logit(0.02),
  #                                      logit_p_sypm = logit(0.5),
  #                                      logit_gamma = logit(1/2),
  #                                      logit_sigma = logit(1/1.5),
  #                                      init_probs_hh = c(0.2, 0.2, 0.2, 0.2, 0.2))), 4))
  
  res_gamma <- readRDS(paste0("../vc-household-working/model_results/final/res_nocov_gammamult_", gamma, ".rds"))

  ch <- rstan::extract(res_gamma)
  ch$gamma <- ch$gamma_sym
  gamma_list[[i]] <- process_runs(chains = ch, cov = FALSE) %>% mutate(asym_recovery_mult = gamma)
}

gamma_out <- bind_rows(gamma_list) %>%
  mutate(param = ifelse(param == "Recovery period", "Symptomatic\nrecovery period", param))
write.csv(gamma_out, "model_results/new_runs/res_gammasens.csv", row.names = FALSE)




