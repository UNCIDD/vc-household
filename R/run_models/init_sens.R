library(rstan)
library(tidyverse)
source("R/helper_scripts/cholera_utilities.R")

init_s <- c(0.8, 0.85, 0.95)
init_list <- list()

for(i in 1:length(init_s)) {
  
  init <- init_s[i]
  dat <- readRDS(paste0("data/dat_11.26.24_nocov_320.rds"))
  dat$init_probs_index <- c(((1-init)-1.5*1e-10)*0.3, 1e-10, init-1.5*1e-10, 1e-10, ((1-init)-1.5*1e-10)*0.7)
  
  res_init <- stan(file = "stan/HMM_cholera_nocov.stan",
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
  
  ch <- rstan::extract(res_init)
  init_list[[i]] <- process_runs(chains = ch, cov = FALSE) %>% mutate(init = init)
}

init_out <- bind_rows(init_list)
write.csv(init_out, "model_results/new_runs/res_initsens.csv", row.names = FALSE)




