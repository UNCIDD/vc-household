library(rstan)
source("R/helper_scripts/cholera_utilities.R")

dat <- readRDS("data/dat_11.26.24_cov_320_new.rds")

# Drop categorical age covariate
dat$x_ih <- dat$x_ih %>% dplyr::select(-age_0.5, -age_5.17, -age_unk)
dat$x_eh <- dat$x_eh %>% dplyr::select(-age_0.5, -age_5.17, -age_unk)
dat$k_ih <- ncol(dat$x_ih)
dat$k_eh <- ncol(dat$x_eh)

# Fully adjusted run
res_cov <- stan(file = "stan/HMM_cholera_cov.stan",
                data = dat,
                cores = 4,
                pars = c("beta_asym", "beta_ih", "beta_eh",
                         "eh_prob", "ih_prob_sym", "ih_prob_asym",
                         "gamma", "sigma", "p_symp", "init_probs_hh"),
                init = rep(list(list(beta_asym = 0,
                                     beta_ih = logit(c(0.02, rep(0.5, dat$k_ih - 1))),
                                     beta_eh = logit(c(0.02, rep(0.5, dat$k_eh - 1))),
                                     logit_p_sypm = logit(0.5),
                                     logit_gamma = logit(1/2),
                                     logit_sigma = logit(1/1.5),
                                     init_probs_hh = c(0.2, 0.2, 0.2, 0.2, 0.2))), 4))
saveRDS(res_cov, paste0("model_results/new_runs/res_cov.rds"))

ch <- rstan::extract(res_cov)
res_multivar <- process_runs(chains = list(ch$beta_eh,
                                           ch$beta_ih),
                             params = c("Extra-household", "Intra-household"),
                             coef = c("Intercept", "Female", "Unknown",
                                      "Private tap", "Tubewell", "Public tap\nor other", "No soap",
                                      "8,500-11,999", "12,000-17,999", "18,000+",
                                      "At home", "Under 18", "Unknown",
                                      "Child", "Sibling", "Spouse", "Unknown or self"),
                             coef_group = c("Intercept", "Sex", "Sex",
                                            rep("Water\nsource", 3), "Soap\nin home",
                                            rep("Monthly\nincome", 3), rep("Occupation\ngroup", 3),
                                            rep("Relation\nto index", 4))) %>%
  mutate(model_type = "Adjusted")

# Unadjusted runs
cov_list <- c("age", "sex", "water", "soap", "income", "occ", "rti")
coef <- list(c("0-4", "5-17", "Unknown"),
             c("Female", "Unknown"),
             c("Private tap", "Tubewell", "Public tap\nor other"),
             c("No soap"),
             c("8,500-11,999", "12,000-17,999", "18,000+"),
             c("At home", "Under 18", "Unknown"),
             c("Child", "Sibling", "Spouse", "Unknown or self"))
coef_group <- c("Age", "Sex", "Water\nsource", "Soap\nin home", "Monthly\nincome", "Occupation\ngroup", "Relation\nto index")
res_univar <- list()
for(i in 1:length(cov_list)) {
  cov_name <- cov_list[i]
  
  dat <- readRDS("data/dat_11.26.24_cov_320_new.rds")
  
  if(cov_name == "age") {
    # age covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, age_0.5, age_5.17, age_unk)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, age_0.5, age_5.17, age_unk)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  if(cov_name == "sex") {
    # sex covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, female, sex_unk)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, female, sex_unk)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  if(cov_name == "water") {
    # water covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, water_private, water_tubewell, water_other)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, water_private, water_tubewell, water_other)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  if(cov_name == "soap") {
    # soap covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, no_soap)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, no_soap)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  if(cov_name == "income") {
    # income covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, inc_8500.11999, inc_12000.17999, inc_18000.)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, inc_8500.11999, inc_12000.17999, inc_18000.)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  if(cov_name == "occ") {
    # occupation covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, occ_home, occ_under18, occ_unk)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, occ_home, occ_under18, occ_unk)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  if(cov_name == "rti") {
    # relation to index covariates
    dat$x_ih <- dat$x_ih %>% dplyr::select(int, rel_child, rel_sibling, rel_spouse, rel_unk_self)
    dat$x_eh <- dat$x_eh %>% dplyr::select(int, rel_child, rel_sibling, rel_spouse, rel_unk_self)
    dat$k_ih <- ncol(dat$x_ih)
    dat$k_eh <- ncol(dat$x_eh)
  }
  
  # res_cov_uni <- stan(file = "code/stan_models/HMM_cholera_cov_fitgammaprior_fitinit_exposed2.stan",
  #                     data = dat,
  #                     cores = 4,
  #                     pars = c("beta_asym", "beta_ih", "beta_eh",
  #                              "eh_prob", "ih_prob_sym", "ih_prob_asym",
  #                              "gamma", "sigma", "p_symp", "init_probs_hh"),
  #                     init = rep(list(list(beta_asym = 0,
  #                                          beta_ih = logit(c(0.02, rep(0.5, dat$k_ih - 1))),
  #                                          beta_eh = logit(c(0.02, rep(0.5, dat$k_eh - 1))),
  #                                          logit_p_sypm = logit(0.5),
  #                                          logit_gamma = logit(1/2),
  #                                          logit_sigma = logit(1/1.5),
  #                                          init_probs_hh = c(0.2, 0.2, 0.2, 0.2, 0.2))), 4))
  # saveRDS(res_cov_uni, paste0("model_results/new_runs/res_cov_", cov_name, ".rds"))
  
  res_cov_uni <- readRDS(paste0("../vc-household-working/model_results/final/res_cov_", cov_name, ".rds"))
  ch <- rstan::extract(res_cov_uni)
  res_univar[[i]] <- process_runs(chains = list(ch$beta_eh,
                                                ch$beta_ih),
                                  params = c("Extra-household", "Intra-household"),
                                  coef = c("Intercept", coef[[i]]),
                                  coef_group = c("Intercept", rep(coef_group[i], length(coef[[i]]))))
  
}

res_univar <- bind_rows(res_univar) %>%
  mutate(model_type = "Unadjusted")

cov_res <- bind_rows(res_multivar,
                     res_univar)

write.csv(cov_res, "model_results/new_runs/res_cov.csv", row.names = FALSE)
