library(rstan)
source("R/helper_scripts/cholera_utilities.R")

dat <- readRDS("data/dat_11.26.24_nocov_320.rds")

res_nocov <- stan(file = "code/stan_models/HMM_cholera_nocov_fitgammaprior_fitinit_exposed2.stan",
                  data = dat,
                  cores = cores,
                  pars = c("beta_asym", "beta_ih", "beta_eh",
                           "eh_prob", "ih_prob_sym", "ih_prob_asym",
                           "gamma", "sigma", "p_symp", "init_probs_hh",
                           "alpha"),
                  iter = 2000,
                  init = rep(list(list(beta_asym = 0,
                                       beta_ih = logit(0.02),
                                       beta_eh = logit(0.02),
                                       logit_p_sypm = logit(0.5),
                                       logit_gamma = logit(1/2),
                                       logit_sigma = logit(1/1.5),
                                       init_probs_hh = c(0.2, 0.2, 0.2, 0.2, 0.2))), 4))

saveRDS(res_nocov, paste0("model_results/new_runs/res_nocov_320_stanfit.rds"))

ch <- rstan::extract(res_nocov)
nocov_320 <- process_runs(chains = ch, cov = FALSE)
write.csv(out_nocov, "model_results/new_runs/res_nocov_320.csv", row.names = FALSE)

