
library(rstan)
library(tidyverse)
source("R/helper_scripts/cholera_utilities.R")

# Simulated data, no covariates
# for(i in 1:10) {
#   dat <- sim_seir(eh_prob = 0.01, ih_prob_asym = 0.025, ih_prob_sym = 0.05,
#                       p_sym = 0.4, phi = 0.17, prop_under1 = 0.3,
#                       n_hh = 500, hh_size = 2:8, days = 30,
#                       gamma = 1/2, sigma = 1/1.5, 
#                       covs_eh = c(0, 0), covs_ih = c(0, 0))
#   
#   saveRDS(dat$complete_obs, paste0("data/simulated/s", i, "_nocov_complete.rds"))
#   saveRDS(dat$dat, paste0("data/simulated/s", i, "_nocov_standat.rds"))
# }

param_ests <- list()
state_probs <- list()
true_sum <- list()

for(snum in 1:10) {
  
  print(snum)
  
  # Code to simulate data is above, uncomment to run your own simulations
  # For replicability, default to the original simulations
  dat <- readRDS(paste0("data/simulated/s", snum, "_nocov_standat.rds"))
  dat$init_probs_hh <- NULL # fititng this
  dat$gamma <- NULL # fitting this
  dat$prop_under1 <- 0.3 # Set proportion of incubation periods under 1 day
  
  res_sim <- stan(file = "code/stan_models/HMM_cholera_nocov_fitgammaprior_fitinit_exposed2.stan",
                    data = dat,
                    cores = 4,
                    pars = c("beta_asym", "beta_ih", "beta_eh",
                             "eh_prob", "ih_prob_sym", "ih_prob_asym",
                             "gamma", "sigma", "p_symp", "init_probs_hh", "alpha"),
                    init = rep(list(list(beta_asym = 0,
                                         beta_ih = logit(0.02),
                                         beta_eh = logit(0.02),
                                         logit_p_sypm = logit(0.5),
                                         logit_gamma = logit(1/2),
                                         logit_sigma = logit(1/1.5),
                                         init_probs_hh = c(0.2, 0.2, 0.2, 0.2, 0.2))), 4))

  ch <- rstan::extract(res_sim)
  
  # extract parameter estimates
  param_ests[[snum]] <- process_runs(chains = ch, cov = FALSE) %>% mutate(snum = snum)
  
  # Computate state probabilities over time
  true <- readRDS(paste0("data/simulated/s", snum, "_nocov_complete.rds")) %>%
    mutate(comp = ifelse(state == 1, "S",
                         ifelse(state == 2, "I_a",
                                ifelse(state == 3, "I_s",
                                       ifelse(state == 4, "R", "E")))))
  
  true_sum[[snum]] <- true %>%
    group_by(comp, day) %>%
    summarize(tot = n()) %>%
    mutate(sim_num = snum)
  
  probs <- data.frame(day = numeric(length = 3*30*sum(dat$hh_size)),
                      pS = numeric(length = 3*30*sum(dat$hh_size)),
                      pI_a = numeric(length = 3*30*sum(dat$hh_size)),
                      pI_s = numeric(length = 3*30*sum(dat$hh_size)),
                      pR = numeric(length = 3*30*sum(dat$hh_size)),
                      pE = numeric(length = 3*30*sum(dat$hh_size)),
                      type = character(length = 3*30*sum(dat$hh_size)))
  
  ct <- 1
  for(p in 1:sum(dat$hh_size)) {
    for(d in 1:30) {
      temp <- ch$alpha[,((p-1)*5+1):((p-1)*5+5),d]
      
      probs$day[ct] <- d
      probs$pS[ct] <- median(temp[,1])
      probs$pI_a[ct] <- median(temp[,2])
      probs$pI_s[ct] <- median(temp[,3])
      probs$pR[ct] <- median(temp[,4])
      probs$pE[ct] <- median(temp[,5])
      probs$type[ct] <- "med"
      ct <- ct+1
      
      probs$day[ct] <- d
      probs$pS[ct] <- quantile(temp[,1], 0.025)
      probs$pI_a[ct] <- quantile(temp[,2], 0.025)
      probs$pI_s[ct] <- quantile(temp[,3], 0.025)
      probs$pR[ct] <- quantile(temp[,4], 0.025)
      probs$pE[ct] <- quantile(temp[,5], 0.025)
      probs$type[ct] <- "low"
      ct <- ct+1
      
      probs$day[ct] <- d
      probs$pS[ct] <- quantile(temp[,1], 0.975)
      probs$pI_a[ct] <- quantile(temp[,2], 0.975)
      probs$pI_s[ct] <- quantile(temp[,3], 0.975)
      probs$pR[ct] <- quantile(temp[,4], 0.975)
      probs$pE[ct] <- quantile(temp[,5], 0.975)
      probs$type[ct] <- "high"
      ct <- ct+1
    }
  }
  
  state_probs[[snum]] <- probs %>%
    pivot_longer(c(-day, -type)) %>%
    group_by(day, type, name) %>%
    summarize(prob = sum(value)) %>%
    mutate(comp = substr(name, 2, nchar(name)),
           sim_num = snum)
}


res_sims <- bind_rows(param_ests) %>%
  mutate(param_type = ifelse(param %in% c("Extra-household", "Intra-household\nasymptomatic", "Intra-household\nsymptomatic"),
                             "Infection probability",
                             ifelse(param %in% c("Pr_S", "Pr_Ia", "Pr_Is", "Pr_R"), "Initial probability", param)))

state_probs <- bind_rows(state_probs)
true_sum <- bind_rows(true_sum)

write.csv(res_sims, "model_results/new_runs/sim_param_ests.csv", row.names = F)
write.csv(state_probs, "model_results/new_runs/sim_stateprobs.rds", row.names = F)
write.csv(true_sum, "model_results/new_runs/sim_truth.rds", row.names = F)


