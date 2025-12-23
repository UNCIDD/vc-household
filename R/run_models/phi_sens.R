library(rstan)
library(tidyverse)
source("R/helper_scripts/cholera_utilities.R")

phis <- c(0.05, 0.10, 0.15, 0.20)
phi_list <- list()

for(i in 1:length(phis)) {
  
  phi <- phis[i]
  dat <- readRDS(paste0("data/dat_11.26.24_nocov_320.rds"))
  
  epsilon <- 1e-10
  dat$obs_prob_s <- matrix(0, nrow = 2, ncol = 5)
  dat$obs_prob_s[1,] <- c(1-phi, 1-phi, epsilon, 1-phi, 1-phi) # outcome = no infection
  dat$obs_prob_s[2,] <- c(phi, phi, 1-epsilon, phi, phi) # outcome = infection
  
  res_phi <- stan(file = "stan/HMM_cholera_nocov.stan",
                  data = dat,
                  cores = cores,
                  pars = c("beta_asym", "beta_ih", "beta_eh",
                           "eh_prob", "ih_prob_sym", "ih_prob_asym",
                           "gamma", "sigma", "p_symp", "init_probs_hh"),
                  iter = 2000,
                  init = rep(list(list(beta_asym = 0,
                                       beta_ih = logit(0.02),
                                       beta_eh = logit(0.02),
                                       logit_p_sypm = logit(0.5),
                                       logit_gamma = logit(1/2),
                                       logit_sigma = logit(1/1.5),
                                       init_probs_hh = c(0.2, 0.2, 0.2, 0.2, 0.2))), 4))

  ch <- rstan::extract(res_phi)
  phi_list[[i]] <- process_runs(chains = ch, cov = FALSE) %>% mutate(phi = phi)
}

phi_out <- bind_rows(phi_list)
write.csv(phi_out, "model_results/new_runs/res_phisens.csv", row.names = FALSE)




