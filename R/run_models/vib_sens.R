library(rstan)
library(tidyverse)
source("R/helper_scripts/cholera_utilities.R")

vib_cutoff <- c(160, 640, 1280)
vib_list <- list()

for(i in 1:length(vib_cutoff)) {
  
  vib <- vib_cutoff[i]
  dat <- readRDS(paste0("data/dat_11.26.24_nocov_",vib,".rds"))
  
  res_vib <- stan(file = "stan/HMM_cholera_nocov.stan",
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

  ch <- rstan::extract(res_vib)
  vib_list[[i]] <- process_runs(chains = ch, cov = FALSE) %>% mutate(vib = vib)
}

vib_out <- bind_rows(vib_list)
write.csv(vib_out, "model_results/new_runs/res_vibsens.csv", row.names = FALSE)




