
require(tidyverse)
library(cowplot)

sim_cholera <- function(eh_prob = 0.01, ih_prob_asym = 0.15, ih_prob_sym = 0.25,
                        p_sym = 0.2, n_hh = 100, hh_size = 2:8, days = 30, phi = 0.01,
                        gamma = 1/5, full_obs = TRUE) {
  
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH
  
  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    enroll_per_hh[i] <- sample(2:hh_size[i], 1)
  }
  
  epsilon <- 1e-10
  enroll_ids <- list()
  
  # Infection state for all household members on all days
  complete_obs <- data.frame(day = numeric(),
                             part_id = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())
  
  for(i in 1:length(hh_size)) {
    
    # Move HH members through SIR states
    for(d in 1:days) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        state = c(3, rep(1, hh_size[i]-1)),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i]),
                                        sample_c = c(1, sample(c(0,1), enroll_per_hh[i]-1, replace = T, prob = c(0.7, 0.3)), rep(0, hh_size[i] - enroll_per_hh[i])),
                                        sample_s = c(1, sample(c(0,1), enroll_per_hh[i]-1, replace = T, prob = c(0.7, 0.3)), rep(0, hh_size[i] - enroll_per_hh[i]))))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_asym <- sum(prior == 2)
        prior_sym <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            no_inf_prob <- (1-eh_prob)*(1-ih_prob_asym)^prior_asym*(1-ih_prob_sym)^prior_sym 
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(no_inf_prob - epsilon/3,
                                               (1-no_inf_prob)*(1-p_sym) - epsilon/3,
                                               (1-no_inf_prob)*p_sym - epsilon/3,
                                               epsilon))
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               1-gamma - epsilon,
                                               epsilon,
                                               gamma - epsilon))
          } else if(prior[part] == 3) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               1-gamma-epsilon,
                                               gamma-epsilon))
          } else {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               epsilon,
                                               1-3*epsilon))
          }
        }
        if(d <= 10) {
          samp_c <- c(sample(c(0,1), enroll_per_hh[i], replace = T, prob = c(0.7, 0.3)), rep(0, hh_size[i] - enroll_per_hh[i]))
          samp_s <- c(sample(c(0,1), enroll_per_hh[i], replace = T, prob = c(0.7, 0.3)), rep(0, hh_size[i] - enroll_per_hh[i]))
        } else if(d == 30) {
          samp_c <- c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i] - enroll_per_hh[i]))
          samp_c <- c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i] - enroll_per_hh[i]))
        }
        
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i]),
                                        sample_c = samp_c,
                                        sample_s = samp_s))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_states
      }
    }
  }
  
  complete_obs <- complete_obs %>%
    arrange(hh_id, day, part_id) %>%
    mutate(row_id = 1:nrow(complete_obs),
           non_cholera_symp = sample(c(0,1), nrow(complete_obs), replace = T, prob = c(1-phi, phi)),
           c = ifelse(state %in% c(2,3), 2, 1),
           s = ifelse(state == 3 | non_cholera_symp == 1, 2, 1))
  
  obs_prob_c <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_c[1,] <- c(1-epsilon, epsilon, epsilon,  1-epsilon) # outcome = no infection
  obs_prob_c[2,] <- c(epsilon, 1-epsilon, 1-epsilon, epsilon) # outcome = infection
  
  obs_prob_s <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_s[1,] <- c(1-phi, 1-phi, epsilon, 1-phi) # outcome = no infection
  obs_prob_s[2,] <- c(phi, phi, 1-epsilon, phi) # outcome = infection
  
  
  if(full_obs) {
    obs_sum <- complete_obs %>%
      group_by(hh_id) %>%
      summarize(obs_per_hh = n(),
                start_ind = min(row_id),
                end_ind = max(row_id))
    
    dat_sim <- list(n_hh =  n_hh,
                    hh_size = hh_size,
                    n_enroll = hh_size,
                    n_obs_c = nrow(complete_obs),
                    n_obs_s = nrow(complete_obs),
                    n_unique_obs_c = 2,
                    n_unique_obs_s = 2,
                    y_c = complete_obs$c,
                    y_s = complete_obs$s,
                    part_id_c = complete_obs$part_id,
                    part_id_s = complete_obs$part_id,
                    t_day_c = complete_obs$day,
                    t_day_s = complete_obs$day,
                    obs_per_hh_c = obs_sum$obs_per_hh,
                    obs_per_hh_s = obs_sum$obs_per_hh,
                    hh_start_ind_c = obs_sum$start_ind,
                    hh_end_ind_c = obs_sum$end_ind,
                    hh_start_ind_s = obs_sum$start_ind,
                    hh_end_ind_s = obs_sum$end_ind,
                    obs_prob_c = obs_prob_c,
                    obs_prob_s = obs_prob_s,
                    init_probs_index = c(epsilon, epsilon, 1-3*epsilon, epsilon),
                    init_probs_hh = c(1-3*epsilon, epsilon, epsilon, epsilon),
                    p_symp = p_sym,
                    gamma = gamma,
                    epsilon = epsilon)
  } else {
    obs <- complete_obs %>% filter(sample_c == 1 | sample_s == 1)
    
    hh_sum <- obs %>%
      group_by(hh_id) %>%
      summarize(n_enroll = length(unique(part_id)))
    
    obs_c <- complete_obs %>% filter(sample_c == 1)
    obs_c$row_id <- 1:nrow(obs_c)
    obs_s <- complete_obs %>% filter(sample_s == 1)
    obs_s$row_id <- 1:nrow(obs_s)
    
    obs_sum_c <- obs_c %>%
      group_by(hh_id) %>%
      summarize(obs_per_hh = n(),
                start_ind = min(row_id),
                end_ind = max(row_id))
    
    obs_sum_s <- obs_s %>%
      group_by(hh_id) %>%
      summarize(obs_per_hh = n(),
                start_ind = min(row_id),
                end_ind = max(row_id))
    
    dat_sim <- list(n_hh =  n_hh,
                    hh_size = hh_size,
                    n_enroll = hh_sum$n_enroll,
                    n_obs_c = nrow(obs_c),
                    n_obs_s = nrow(obs_s),
                    n_unique_obs_c = 2,
                    n_unique_obs_s = 2,
                    y_c = obs_c$c,
                    y_s = obs_s$s,
                    part_id_c = obs_c$part_id,
                    part_id_s = obs_s$part_id,
                    t_day_c = obs_c$day,
                    t_day_s = obs_s$day,
                    obs_per_hh_c = obs_sum_c$obs_per_hh,
                    obs_per_hh_s = obs_sum_s$obs_per_hh,
                    hh_start_ind_c = obs_sum_c$start_ind,
                    hh_end_ind_c = obs_sum_c$end_ind,
                    hh_start_ind_s = obs_sum_s$start_ind,
                    hh_end_ind_s = obs_sum_s$end_ind,
                    obs_prob_c = obs_prob_c,
                    obs_prob_s = obs_prob_s,
                    init_probs_index = c(epsilon, epsilon, 1-3*epsilon, epsilon),
                    init_probs_hh = c(1-3*epsilon, epsilon, epsilon, epsilon),
                    p_symp = p_sym,
                    gamma = gamma,
                    epsilon = epsilon)
  }
  
  return(dat_sim)
  
}

process_runs <- function(chains, params = NULL, coef = NULL, coef_group = NULL, cov = TRUE) {
  
  if(cov) {
    beta_name <- coef
    beta_group <- coef_group
    param_name <- params
    
    plt_dat <- data.frame(Est = numeric(0),
                          Est_med = numeric(0),
                          CI_high = numeric(0),
                          CI_low = numeric(0),
                          cov = character(0),
                          cov_group = character(0),
                          param = character(0))
    
    for(p in 1:length(params)) {
      for(b in 1:length(coef)) {
        plt_dat <- plt_dat %>%
          bind_rows(data.frame(Est = exp(mean(chains[[p]][,b])),
                               Est_med = exp(median(chains[[p]][,b])),
                               CI_high = exp(quantile(chains[[p]][,b], .975)),
                               CI_low = exp(quantile(chains[[p]][,b], .025)),
                               cov = coef[b],
                               cov_group = coef_group[b],
                               param = params[p]))
      }
    }
    
    refs <- unique(coef_group[coef_group != "Intercept"])
    ref_names <- c("Age" = "18+",
                   "Sex" = "Male",
                   "Water\nsource" = "Boiled\nwater",
                   "Soap\nin home" = "Soap",
                   "Relation\nto index" = "Parent",
                   "Monthly\nincome" = "0-8,499",
                   "Occupation\ngroup" = "Outside home")
    ref_labs <- ref_names[refs]
    
    
    plt_dat <- plt_dat %>%
      bind_rows(data.frame(Est = rep(rep(1, length(refs)), length(params)),
                           Est_med = rep(rep(NA, length(refs)), length(params)),
                           CI_high = rep(rep(NA, length(refs)), length(params)),
                           CI_low = rep(rep(NA, length(refs)), length(params)),
                           cov = rep(ref_labs, length(params)),
                           cov_group = rep(refs, length(params)),
                           param = rep(param_name, each = length(refs)))) %>%
      mutate(cov_group = factor(cov_group, levels = c(unique(coef_group), "Symptom status")),
             cov = factor(cov, levels = c("Unknown", unname(ref_names), unique(coef[which(coef != "Unknown")]), "Asymptomatic", "Symptomatic")))  
  } else {
    plt_dat <- data.frame(param = c("Intra-household\nsymptomatic", "Intra-household\nasymptomatic", "OR - asmpt/symp",
                                    "Extra-household", "Recovery period", "Incubation period", "Probability of symptoms",
                                    "Pr_S", "Pr_Ia", "Pr_Is", "Pr_R", "Pr_E"),
                          est = c(inv_logit(mean(chains$beta_ih)),
                                  inv_logit(mean(chains$beta_ih+chains$beta_asym)),
                                  exp(mean(chains$beta_asym)),
                                  inv_logit(mean(chains$beta_eh)),
                                  mean(chains$gamma),
                                  mean(chains$sigma),
                                  mean(chains$p_symp),
                                  hh_s = mean(chains$init_probs_hh[,1]),
                                  hh_ia = mean(chains$init_probs_hh[,2]),
                                  hh_is = mean(chains$init_probs_hh[,3]),
                                  hh_r = mean(chains$init_probs_hh[,4]),
                                  hh_r = mean(chains$init_probs_hh[,5])),
                          est_med = c(inv_logit(median(chains$beta_ih)),
                                      inv_logit(median(chains$beta_ih+chains$beta_asym)),
                                      exp(median(chains$beta_asym)),
                                      inv_logit(median(chains$beta_eh)),
                                      median(chains$gamma),
                                      median(chains$sigma),
                                      median(chains$p_symp),
                                      hh_s = median(chains$init_probs_hh[,1]),
                                      hh_ia = median(chains$init_probs_hh[,2]),
                                      hh_is = median(chains$init_probs_hh[,3]),
                                      hh_r = median(chains$init_probs_hh[,4]),
                                      hh_r = median(chains$init_probs_hh[,5])),
                          ci_high = c(inv_logit(quantile(chains$beta_ih, 0.975)),
                                      inv_logit(quantile(chains$beta_ih+chains$beta_asym, 0.975)),
                                      exp(quantile(chains$beta_asym, 0.975)),
                                      inv_logit(quantile(chains$beta_eh, 0.975)),
                                      quantile(chains$gamma, 0.975),
                                      quantile(chains$sigma, 0.975),
                                      quantile(chains$p_symp, 0.975),
                                      hh_s = quantile(chains$init_probs_hh[,1], 0.025),
                                      hh_ia = quantile(chains$init_probs_hh[,2], 0.025),
                                      hh_is = quantile(chains$init_probs_hh[,3], 0.025),
                                      hh_r = quantile(chains$init_probs_hh[,4], 0.025),
                                      hh_r = quantile(chains$init_probs_hh[,5], 0.025)),
                          ci_low = c(inv_logit(quantile(chains$beta_ih, 0.025)),
                                     inv_logit(quantile(chains$beta_ih+chains$beta_asym, 0.025)),
                                     exp(quantile(chains$beta_asym, 0.025)),
                                     inv_logit(quantile(chains$beta_eh, 0.025)),
                                     quantile(chains$gamma, 0.025),
                                     quantile(chains$sigma, 0.025),
                                     quantile(chains$p_symp, 0.025),
                                     hh_s = quantile(chains$init_probs_hh[,1], 0.975),
                                     hh_ia = quantile(chains$init_probs_hh[,2], 0.975),
                                     hh_is = quantile(chains$init_probs_hh[,3], 0.975),
                                     hh_r = quantile(chains$init_probs_hh[,4], 0.975),
                                     hh_e = quantile(chains$init_probs_hh[,5], 0.975)))
    
  }
  
  
  
  return(plt_dat)
  
}

inv_logit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

logit <- function(p) {
  return(log(p/(1-p)))
}

sim_cholera_vib <- function(eh_prob = 0.01, ih_prob_asym = 0.15, ih_prob_sym = 0.25,
                            p_sym = 0.4, n_hh = 100, hh_size = 2:8, days = 30, phi = 0.01,
                            gamma = 1/2, full_obs = TRUE) {


  
  epsilon <- 1e-10
  sens <- 0.82
  spec <- 0.943
  obs_prob_c <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_c[1,] <- c(spec, 1-sens, 1-sens,  spec) # outcome = no infection
  obs_prob_c[2,] <- c(1-spec, sens, sens, 1-spec) # outcome = infection
  
  obs_prob_s <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_s[1,] <- c(1-phi, 1-phi, epsilon, 1-phi) # outcome = no infection
  obs_prob_s[2,] <- c(phi, phi, 1-epsilon, phi) # outcome = infection
  
  obs_prob_v <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_v[1,] <- c(1-epsilon, 0.9, 0.9, 0.3) # outcome = low titer
  obs_prob_v[2,] <- c(epsilon, 0.1, 0.1, 0.7) # outcome = high titer
  
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH
  
  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    enroll_per_hh[i] <- sample(2:hh_size[i], 1)
  }
  
  epsilon <- 1e-10
  enroll_ids <- list()
  
  # Infection state for all household members on all days
  complete_obs <- data.frame(day = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())
  
  for(i in 1:length(hh_size)) {
    
    # Move HH members through SIR states
    for(d in 1:days) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = c(3, sample(1:4, hh_size[i]-1, replace = T, prob = c(0.4, 0.1, 0.1, 0.4))),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_asym <- sum(prior == 2)
        prior_sym <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            no_inf_prob <- (1-eh_prob)*(1-ih_prob_asym)^prior_asym*(1-ih_prob_sym)^prior_sym 
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(no_inf_prob - epsilon/3,
                                               (1-no_inf_prob)*(1-p_sym) - epsilon/3,
                                               (1-no_inf_prob)*p_sym - epsilon/3,
                                               epsilon))

           
              
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               1-gamma - epsilon,
                                               epsilon,
                                               gamma - epsilon))
          } else if(prior[part] == 3) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               1-gamma-epsilon,
                                               gamma-epsilon))
          } else {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               epsilon,
                                               1-3*epsilon))
          }
        }
        
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_states
      }
    }
  }
  
  complete_obs <- complete_obs %>%
    arrange(hh_id, day, part_id)
  
  obs <- complete_obs %>% filter(enroll == 1)
  
  hh_sum <- obs %>%
    group_by(hh_id) %>%
    summarize(n_enroll = length(unique(part_id)))
    
  obs_c <- complete_obs %>% filter(day %in% 2:10 | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_c$row_id <- 1:nrow(obs_c)
  obs_v <- complete_obs %>% filter(day %in% c(2, 7, 30)) %>%
    mutate(obs = 1)
  obs_v$row_id <- 1:nrow(obs_v)
  obs_s <- complete_obs %>% filter(day %in% c(2:10, 30) | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_s$row_id <- 1:nrow(obs_s)
    
  for(i in 1:nrow(obs_c)) {
    obs_c$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_c[,obs_c$state[i]]))
  }
  for(i in 1:nrow(obs_v)) {
    obs_v$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_v[,obs_v$state[i]]))
  }
  for(i in 1:nrow(obs_s)) {
    obs_s$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_s[,obs_s$state[i]]))
  }
  
  obs_sum_c <- obs_c %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
    
  obs_sum_s <- obs_s %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_v <- obs_v %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
    
  dat_sim <- list(n_hh =  n_hh,
                  hh_size = hh_size,
                  n_enroll = hh_sum$n_enroll,
                  n_obs_c = nrow(obs_c),
                  n_obs_s = nrow(obs_s),
                  n_obs_v = nrow(obs_v),
                  n_unique_obs_c = 2,
                  n_unique_obs_s = 2,
                  n_unique_obs_v = 2,
                  y_c = obs_c$obs,
                  y_s = obs_s$obs,
                  y_v = obs_v$obs,
                  part_id_c = obs_c$part_id,
                  part_id_s = obs_s$part_id,
                  part_id_v = obs_v$part_id,
                  t_day_c = obs_c$day,
                  t_day_s = obs_s$day,
                  t_day_v = obs_v$day,
                  obs_per_hh_c = obs_sum_c$obs_per_hh,
                  obs_per_hh_s = obs_sum_s$obs_per_hh,
                  obs_per_hh_v = obs_sum_v$obs_per_hh,
                  hh_start_ind_c = obs_sum_c$start_ind,
                  hh_end_ind_c = obs_sum_c$end_ind,
                  hh_start_ind_s = obs_sum_s$start_ind,
                  hh_end_ind_s = obs_sum_s$end_ind,
                  hh_start_ind_v = obs_sum_v$start_ind,
                  hh_end_ind_v = obs_sum_v$end_ind,
                  obs_prob_c = obs_prob_c,
                  obs_prob_s = obs_prob_s,
                  obs_prob_v = obs_prob_v,
                  init_probs_index = c(epsilon, epsilon, 1-3*epsilon, epsilon),
                  init_probs_hh = c(0.4,0.1, 0.1, 0.4),
                  p_symp = p_sym,
                  gamma = gamma,
                  epsilon = epsilon)

 return(list(dat = dat_sim,
             obs = obs,
             complete_obs = complete_obs))
  
}

sim_cholera_vib_cov <- function(eh_prob = 0.01, ih_prob_asym = 0.025, ih_prob_sym = 0.05,
                            p_sym = 0.4, n_hh = 100, hh_size = 2:8, days = 30, phi = 0.01,
                            gamma = 1/2, covs_eh = c(0.7, 1.6), covs_ih = c(0.6, 1.4)) {

  ##Simple infection/exposure matrix. This stores the exposures and 
  ##the results of a particular infection event.
  infect_evt_df <- NULL
  
  
  
  epsilon <- 1e-10
  sens <- 0.82
  spec <- 0.943
  obs_prob_c <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_c[1,] <- c(spec, 1-sens, 1-sens,  spec) # outcome = no infection
  obs_prob_c[2,] <- c(1-spec, sens, sens, 1-spec) # outcome = infection
  
  obs_prob_s <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_s[1,] <- c(1-phi, 1-phi, epsilon, 1-phi) # outcome = no infection
  obs_prob_s[2,] <- c(phi, phi, 1-epsilon, phi) # outcome = infection
  
  obs_prob_v <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_v[1,] <- c(1-epsilon, 0.9, 0.9, 0.3) # outcome = low titer
  obs_prob_v[2,] <- c(epsilon, 0.1, 0.1, 0.7) # outcome = high titer
  
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH
  
  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))
  
  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4)
  }
  
  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    enroll_per_hh[i] <- sample(2:hh_size[i], 1)
  }
  
  epsilon <- 1e-10
  enroll_ids <- list()
  
  # Infection state for all household members on all days
  complete_obs <- data.frame(day = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())
  
  last_x <- 0
  
  for(i in 1:length(hh_size)) {
    
    # Move HH members through SIR states
    for(d in 1:days) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = c(3, sample(1:4, hh_size[i]-1, replace = T, prob = c(0.4, 0.1, 0.1, 0.4))),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_asym <- sum(prior == 2)
        prior_sym <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            eh_prob_x <- inv_logit(logit(eh_prob)+sum(x[last_x+part,]*covs_eh))
            ih_prob_asym_x <- inv_logit(logit(ih_prob_asym)+sum(x[last_x+part,]*covs_ih))
            ih_prob_sym_x <- inv_logit(logit(ih_prob_sym)+sum(x[last_x+part,]*covs_ih))
            no_inf_prob <- (1-eh_prob_x)*(1-ih_prob_asym_x)^prior_asym*(1-ih_prob_sym_x)^prior_sym 
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(no_inf_prob - epsilon/3,
                                               (1-no_inf_prob)*(1-p_sym) - epsilon/3,
                                               (1-no_inf_prob)*p_sym - epsilon/3,
                                               epsilon))

             infect_evt_df<- bind_rows(infect_evt_df, tibble (
                inf = new_states[part] == 2,
                asym =  prior_asym,
                sym = prior_sym
              )
              )
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               1-gamma - epsilon,
                                               epsilon,
                                               gamma - epsilon))
          } else if(prior[part] == 3) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               1-gamma-epsilon,
                                               gamma-epsilon))
          } else {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               epsilon,
                                               1-3*epsilon))
          }
        }
        
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }
  
  complete_obs <- complete_obs %>%
    arrange(hh_id, day, part_id)
  
  obs <- complete_obs %>% filter(enroll == 1)
  
  hh_sum <- obs %>%
    group_by(hh_id) %>%
    summarize(n_enroll = length(unique(part_id)))
  
  obs_c <- complete_obs %>% filter(day %in% 2:10 | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_c$row_id <- 1:nrow(obs_c)
  obs_v <- complete_obs %>% filter(day %in% c(2, 7, 30)) %>%
    mutate(obs = 1)
  obs_v$row_id <- 1:nrow(obs_v)
  obs_s <- complete_obs %>% filter(day %in% c(2:10, 30) | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_s$row_id <- 1:nrow(obs_s)
  
  for(i in 1:nrow(obs_c)) {
    obs_c$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_c[,obs_c$state[i]]))
  }
  for(i in 1:nrow(obs_v)) {
    obs_v$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_v[,obs_v$state[i]]))
  }
  for(i in 1:nrow(obs_s)) {
    obs_s$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_s[,obs_s$state[i]]))
  }
  
  obs_sum_c <- obs_c %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_s <- obs_s %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_v <- obs_v %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  dat_sim <- list(n_hh =  n_hh,
                  hh_size = hh_size,
                  n_enroll = hh_sum$n_enroll,
                  n_obs_c = nrow(obs_c),
                  n_obs_s = nrow(obs_s),
                  n_obs_v = nrow(obs_v),
                  n_unique_obs_c = 2,
                  n_unique_obs_s = 2,
                  n_unique_obs_v = 2,
                  y_c = obs_c$obs,
                  y_s = obs_s$obs,
                  y_v = obs_v$obs,
                  part_id_c = obs_c$part_id,
                  part_id_s = obs_s$part_id,
                  part_id_v = obs_v$part_id,
                  t_day_c = obs_c$day,
                  t_day_s = obs_s$day,
                  t_day_v = obs_v$day,
                  obs_per_hh_c = obs_sum_c$obs_per_hh,
                  obs_per_hh_s = obs_sum_s$obs_per_hh,
                  obs_per_hh_v = obs_sum_v$obs_per_hh,
                  hh_start_ind_c = obs_sum_c$start_ind,
                  hh_end_ind_c = obs_sum_c$end_ind,
                  hh_start_ind_s = obs_sum_s$start_ind,
                  hh_end_ind_s = obs_sum_s$end_ind,
                  hh_start_ind_v = obs_sum_v$start_ind,
                  hh_end_ind_v = obs_sum_v$end_ind,
                  obs_prob_c = obs_prob_c,
                  obs_prob_s = obs_prob_s,
                  obs_prob_v = obs_prob_v,
                  init_probs_index = c(epsilon, epsilon, 1-3*epsilon, epsilon),
                  init_probs_hh = c(0.4,0.1, 0.1, 0.4),
                  p_symp = p_sym,
                  gamma = gamma,
                  epsilon = epsilon,
                  x_eh = cbind(rep(1, sum(hh_size)), x),
                  x_ih = cbind(rep(1, sum(hh_size)), x),
                  k_eh = ncol(x)+1,
                  k_ih = ncol(x)+1,
                  n_days = 30)
  
  return(list(dat = dat_sim,
              obs = obs,
              complete_obs = complete_obs,
              infect_evt_df = infect_evt_df))
  
}

process_sims_nocov <- function(sims) {
  out <- data.frame(est = numeric(),
                    ci_high = numeric(),
                    ci_low = numeric(),
                    param = character(),
                    sim_num = numeric())
  
  for(i in 1:length(sims)) {
    ch <- extract(sims[[i]])
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_ih)),
                                ci_low = inv_logit(quantile(ch$beta_ih, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_ih, 0.975)),
                                param = "Intra-household - symptomatic",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_ih + ch$beta_asym)),
                                ci_low = inv_logit(quantile(ch$beta_ih + ch$beta_asym, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_ih + ch$beta_asym, 0.975)),
                                param = "Intra-household - asymptomatic",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_eh)),
                                ci_low = inv_logit(quantile(ch$beta_eh, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_eh, 0.975)),
                                param = "Extra-household",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = exp(mean(ch$beta_asym)),
                                ci_low = exp(quantile(ch$beta_asym, 0.025)),
                                ci_high = exp(quantile(ch$beta_asym, 0.975)),
                                param = "Asym vs. sym OR",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = mean(ch$p_symp),
                                ci_low = quantile(ch$p_symp, 0.025),
                                ci_high = quantile(ch$p_symp, 0.975),
                                param = "Symptomatic proportion",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = mean(1/ch$gamma),
                                ci_low = quantile(1/ch$gamma, 0.025),
                                ci_high = quantile(1/ch$gamma, 0.975),
                                param = "Duration of infectiousness",
                                sim_num = i))
  }
  return(out)
}

sim_cholera_vib_cov <- function(eh_prob = 0.01, ih_prob_asym = 0.025, ih_prob_sym = 0.05,
                            p_sym = 0.4, n_hh = 100, hh_size = 2:8, days = 30, phi = 0.01,
                            gamma = 1/2, covs_eh = c(0.7, 1.6), covs_ih = c(0.6, 1.4)) {

  ##Simple infection/exposure matrix. This stores the exposures and 
  ##the results of a particular infection event.
  infect_evt_df <- NULL
  
  
  
  epsilon <- 1e-10
  sens <- 0.82
  spec <- 0.943
  obs_prob_c <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_c[1,] <- c(spec, 1-sens, 1-sens,  spec) # outcome = no infection
  obs_prob_c[2,] <- c(1-spec, sens, sens, 1-spec) # outcome = infection
  
  obs_prob_s <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_s[1,] <- c(1-phi, 1-phi, epsilon, 1-phi) # outcome = no infection
  obs_prob_s[2,] <- c(phi, phi, 1-epsilon, phi) # outcome = infection
  
  obs_prob_v <- matrix(0, nrow = 2, ncol = 4)
  obs_prob_v[1,] <- c(1-epsilon, 0.9, 0.9, 0.3) # outcome = low titer
  obs_prob_v[2,] <- c(epsilon, 0.1, 0.1, 0.7) # outcome = high titer
  
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH
  
  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))
  
  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4)
  }
  
  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    enroll_per_hh[i] <- sample(2:hh_size[i], 1)
  }
  
  epsilon <- 1e-10
  enroll_ids <- list()
  
  # Infection state for all household members on all days
  complete_obs <- data.frame(day = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())
  
  last_x <- 0
  
  for(i in 1:length(hh_size)) {
    
    # Move HH members through SIR states
    for(d in 1:days) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = c(3, sample(1:4, hh_size[i]-1, replace = T, prob = c(0.4, 0.1, 0.1, 0.4))),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_asym <- sum(prior == 2)
        prior_sym <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            eh_prob_x <- inv_logit(logit(eh_prob)+sum(x[last_x+part,]*covs_eh))
            ih_prob_asym_x <- inv_logit(logit(ih_prob_asym)+sum(x[last_x+part,]*covs_ih))
            ih_prob_sym_x <- inv_logit(logit(ih_prob_sym)+sum(x[last_x+part,]*covs_ih))
            no_inf_prob <- (1-eh_prob_x)*(1-ih_prob_asym_x)^prior_asym*(1-ih_prob_sym_x)^prior_sym 
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(no_inf_prob - epsilon/3,
                                               (1-no_inf_prob)*(1-p_sym) - epsilon/3,
                                               (1-no_inf_prob)*p_sym - epsilon/3,
                                               epsilon))

             infect_evt_df<- bind_rows(infect_evt_df, tibble (
                inf = new_states[part] == 2,
                asym =  prior_asym,
                sym = prior_sym
              )
              )
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               1-gamma - epsilon,
                                               epsilon,
                                               gamma - epsilon))
          } else if(prior[part] == 3) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               1-gamma-epsilon,
                                               gamma-epsilon))
          } else {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               epsilon,
                                               1-3*epsilon))
          }
        }
        
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }
  
  complete_obs <- complete_obs %>%
    arrange(hh_id, day, part_id)
  
  obs <- complete_obs %>% filter(enroll == 1)
  
  hh_sum <- obs %>%
    group_by(hh_id) %>%
    summarize(n_enroll = length(unique(part_id)))
  
  obs_c <- complete_obs %>% filter(day %in% 2:10 | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_c$row_id <- 1:nrow(obs_c)
  obs_v <- complete_obs %>% filter(day %in% c(2, 7, 30)) %>%
    mutate(obs = 1)
  obs_v$row_id <- 1:nrow(obs_v)
  obs_s <- complete_obs %>% filter(day %in% c(2:10, 30) | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_s$row_id <- 1:nrow(obs_s)
  
  for(i in 1:nrow(obs_c)) {
    obs_c$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_c[,obs_c$state[i]]))
  }
  for(i in 1:nrow(obs_v)) {
    obs_v$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_v[,obs_v$state[i]]))
  }
  for(i in 1:nrow(obs_s)) {
    obs_s$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_s[,obs_s$state[i]]))
  }
  
  obs_sum_c <- obs_c %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_s <- obs_s %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_v <- obs_v %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  dat_sim <- list(n_hh =  n_hh,
                  hh_size = hh_size,
                  n_enroll = hh_sum$n_enroll,
                  n_obs_c = nrow(obs_c),
                  n_obs_s = nrow(obs_s),
                  n_obs_v = nrow(obs_v),
                  n_unique_obs_c = 2,
                  n_unique_obs_s = 2,
                  n_unique_obs_v = 2,
                  y_c = obs_c$obs,
                  y_s = obs_s$obs,
                  y_v = obs_v$obs,
                  part_id_c = obs_c$part_id,
                  part_id_s = obs_s$part_id,
                  part_id_v = obs_v$part_id,
                  t_day_c = obs_c$day,
                  t_day_s = obs_s$day,
                  t_day_v = obs_v$day,
                  obs_per_hh_c = obs_sum_c$obs_per_hh,
                  obs_per_hh_s = obs_sum_s$obs_per_hh,
                  obs_per_hh_v = obs_sum_v$obs_per_hh,
                  hh_start_ind_c = obs_sum_c$start_ind,
                  hh_end_ind_c = obs_sum_c$end_ind,
                  hh_start_ind_s = obs_sum_s$start_ind,
                  hh_end_ind_s = obs_sum_s$end_ind,
                  hh_start_ind_v = obs_sum_v$start_ind,
                  hh_end_ind_v = obs_sum_v$end_ind,
                  obs_prob_c = obs_prob_c,
                  obs_prob_s = obs_prob_s,
                  obs_prob_v = obs_prob_v,
                  init_probs_index = c(epsilon, epsilon, 1-3*epsilon, epsilon),
                  init_probs_hh = c(0.4,0.1, 0.1, 0.4),
                  p_symp = p_sym,
                  gamma = gamma,
                  epsilon = epsilon,
                  x_eh = cbind(rep(1, sum(hh_size)), x),
                  x_ih = cbind(rep(1, sum(hh_size)), x),
                  k_eh = ncol(x)+1,
                  k_ih = ncol(x)+1,
                  n_days = 30)
  
  return(list(dat = dat_sim,
              obs = obs,
              complete_obs = complete_obs,
              infect_evt_df = infect_evt_df))
  
}

process_sims_nocov <- function(sims) {
  out <- data.frame(est = numeric(),
                    ci_high = numeric(),
                    ci_low = numeric(),
                    param = character(),
                    sim_num = numeric())
  
  for(i in 1:length(sims)) {
    ch <- extract(sims[[i]])
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_ih)),
                                ci_low = inv_logit(quantile(ch$beta_ih, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_ih, 0.975)),
                                param = "Intra-household - symptomatic",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_ih + ch$beta_asym)),
                                ci_low = inv_logit(quantile(ch$beta_ih + ch$beta_asym, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_ih + ch$beta_asym, 0.975)),
                                param = "Intra-household - asymptomatic",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_eh)),
                                ci_low = inv_logit(quantile(ch$beta_eh, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_eh, 0.975)),
                                param = "Extra-household",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = exp(mean(ch$beta_asym)),
                                ci_low = exp(quantile(ch$beta_asym, 0.025)),
                                ci_high = exp(quantile(ch$beta_asym, 0.975)),
                                param = "Asym vs. sym OR",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = mean(ch$p_symp),
                                ci_low = quantile(ch$p_symp, 0.025),
                                ci_high = quantile(ch$p_symp, 0.975),
                                param = "Symptomatic proportion",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = mean(1/ch$gamma),
                                ci_low = quantile(1/ch$gamma, 0.025),
                                ci_high = quantile(1/ch$gamma, 0.975),
                                param = "Duration of infectiousness",
                                sim_num = i))
  }
  return(out)
}

sim_seir <- function(eh_prob = 0.01, ih_prob_asym = 0.025, ih_prob_sym = 0.05,
                     p_sym = 0.4, phi = 0.01, prop_under1 = 0.3,
                     n_hh = 100, hh_size = 2:8, days = 30,
                     gamma = 1/2, sigma = 1/1.5, 
                     covs_eh = c(0.7, 1.6), covs_ih = c(0.6, 1.4)) {
  
  epsilon <- 1e-10
  sens <- 0.82
  spec <- 0.943
  obs_prob_c <- matrix(0, nrow = 2, ncol = 5)
  obs_prob_c[1,] <- c(spec, 1-sens, 1-sens,  spec, spec) # outcome = no infection
  obs_prob_c[2,] <- c(1-spec, sens, sens, 1-spec, 1-spec) # outcome = infection
  
  obs_prob_s <- matrix(0, nrow = 2, ncol = 5)
  obs_prob_s[1,] <- c(1-phi, 1-phi, epsilon, 1-phi, 1-phi) # outcome = no infection
  obs_prob_s[2,] <- c(phi, phi, 1-epsilon, phi, phi) # outcome = infection
  
  obs_prob_v <- matrix(0, nrow = 2, ncol = 5)
  obs_prob_v[1,] <- c(1-epsilon, 0.9, 0.9, 0.3, 1-epsilon) # outcome = low titer
  obs_prob_v[2,] <- c(epsilon, 0.1, 0.1, 0.7, epsilon) # outcome = high titer
  
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH
  
  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))
  
  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4)
  }
  
  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    enroll_per_hh[i] <- sample(2:hh_size[i], 1)
  }
  
  epsilon <- 1e-10
  enroll_ids <- list()
  
  # Infection state for all household members on all days
  complete_obs <- data.frame(day = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())
  
  last_x <- 0

  for(i in 1:length(hh_size)) {
    
    # Move HH members through SIR states
    for(d in 1:days) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = c(3, sample(1:5, hh_size[i]-1, replace = T, prob = c(0.85/2, 0.05, 0.05, 0.85/2, 0.05))),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_asym <- sum(prior == 2)
        prior_sym <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) { # start in S
            eh_prob_x <- inv_logit(logit(eh_prob)+sum(x[last_x+part,]*covs_eh))
            ih_prob_asym_x <- inv_logit(logit(ih_prob_asym)+sum(x[last_x+part,]*covs_ih))
            ih_prob_sym_x <- inv_logit(logit(ih_prob_sym)+sum(x[last_x+part,]*covs_ih))
            no_inf_prob <- (1-eh_prob_x)*(1-ih_prob_asym_x)^prior_asym*(1-ih_prob_sym_x)^prior_sym 
            new_states[part] = sample(x = c(1, 2, 3, 4, 5),
                                      size = 1,
                                      prob = c(no_inf_prob - epsilon/4,
                                               (1-no_inf_prob)*(1-p_sym)*prop_under1 - epsilon/4,
                                               (1-no_inf_prob)*p_sym*prop_under1 - epsilon/4,
                                               epsilon,
                                               (1-prop_under1)*(1-no_inf_prob)-epsilon/4))
            
          } else if(prior[part] == 2) { # start in Ia
            new_states[part] = sample(x = c(1, 2, 3, 4, 5),
                                      size = 1,
                                      prob = c(epsilon,
                                               1-gamma - 3*epsilon/2,
                                               epsilon,
                                               gamma - 3*epsilon/2,
                                               epsilon))
          } else if(prior[part] == 3) { # start in Is
            new_states[part] = sample(x = c(1, 2, 3, 4, 5),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               1-gamma - 3*epsilon/2,
                                               gamma - 3*epsilon/2,
                                               epsilon))
          } else if(prior[part] == 4) { # start in R
            new_states[part] = sample(x = c(1, 2, 3, 4, 5),
                                      size = 1,
                                      prob = c(epsilon,
                                               epsilon,
                                               epsilon,
                                               1-4*epsilon,
                                               epsilon))
          } else { # start in E
            new_states[part] = sample(x = c(1, 2, 3, 4, 5),
                                      size = 1,
                                      prob = c(epsilon,
                                               (1-p_sym)*sigma - 2*epsilon/3,
                                               p_sym*sigma - 2*epsilon/3,
                                               epsilon,
                                               1-sigma - 2*epsilon/3))
          }
        }
        
        new_obs <- bind_rows(data.frame(day = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))
        
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }
  
  complete_obs <- complete_obs %>%
    arrange(hh_id, day, part_id)
  
  obs <- complete_obs %>% filter(enroll == 1)
  
  hh_sum <- obs %>%
    group_by(hh_id) %>%
    summarize(n_enroll = length(unique(part_id)))
  
  obs_c <- complete_obs %>% filter(day %in% 2:10 | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_c$row_id <- 1:nrow(obs_c)
  obs_v <- complete_obs %>% filter(day %in% c(2, 7, 30)) %>%
    mutate(obs = 1)
  obs_v$row_id <- 1:nrow(obs_v)
  obs_s <- complete_obs %>% filter(day %in% c(2:10, 30) | (day == 1 & part_id == 1)) %>%
    mutate(obs = 1)
  obs_s$row_id <- 1:nrow(obs_s)
  
  for(i in 1:nrow(obs_c)) {
    obs_c$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_c[,obs_c$state[i]]))
  }
  for(i in 1:nrow(obs_v)) {
    obs_v$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_v[,obs_v$state[i]]))
  }
  for(i in 1:nrow(obs_s)) {
    obs_s$obs[i] <- sample(1:2, 1, prob = as.vector(obs_prob_s[,obs_s$state[i]]))
  }
  
  obs_sum_c <- obs_c %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_s <- obs_s %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  obs_sum_v <- obs_v %>%
    group_by(hh_id) %>%
    summarize(obs_per_hh = n(),
              start_ind = min(row_id),
              end_ind = max(row_id))
  
  dat_sim <- list(n_hh =  n_hh,
                  hh_size = hh_size,
                  n_enroll = hh_sum$n_enroll,
                  n_obs_c = nrow(obs_c),
                  n_obs_s = nrow(obs_s),
                  n_obs_v = nrow(obs_v),
                  n_unique_obs_c = 2,
                  n_unique_obs_s = 2,
                  n_unique_obs_v = 2,
                  y_c = obs_c$obs,
                  y_s = obs_s$obs,
                  y_v = obs_v$obs,
                  part_id_c = obs_c$part_id,
                  part_id_s = obs_s$part_id,
                  part_id_v = obs_v$part_id,
                  t_day_c = obs_c$day,
                  t_day_s = obs_s$day,
                  t_day_v = obs_v$day,
                  obs_per_hh_c = obs_sum_c$obs_per_hh,
                  obs_per_hh_s = obs_sum_s$obs_per_hh,
                  obs_per_hh_v = obs_sum_v$obs_per_hh,
                  hh_start_ind_c = obs_sum_c$start_ind,
                  hh_end_ind_c = obs_sum_c$end_ind,
                  hh_start_ind_s = obs_sum_s$start_ind,
                  hh_end_ind_s = obs_sum_s$end_ind,
                  hh_start_ind_v = obs_sum_v$start_ind,
                  hh_end_ind_v = obs_sum_v$end_ind,
                  obs_prob_c = obs_prob_c,
                  obs_prob_s = obs_prob_s,
                  obs_prob_v = obs_prob_v,
                  init_probs_index = c(epsilon, epsilon, 1-4*epsilon, epsilon, epsilon),
                  init_probs_hh = c(0.85/2, 0.05, 0.05, 0.85/2, 0.05),
                  p_symp = p_sym,
                  gamma = gamma,
                  epsilon = epsilon,
                  x_eh = cbind(rep(1, sum(hh_size)), x),
                  x_ih = cbind(rep(1, sum(hh_size)), x),
                  k_eh = ncol(x)+1,
                  k_ih = ncol(x)+1,
                  n_days = 30)
  
  return(list(dat = dat_sim,
              obs = obs,
              complete_obs = complete_obs))
  
}

process_sims_nocov <- function(sims) {
  out <- data.frame(est = numeric(),
                    ci_high = numeric(),
                    ci_low = numeric(),
                    param = character(),
                    sim_num = numeric())
  
  for(i in 1:length(sims)) {
    ch <- extract(sims[[i]])
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_ih)),
                                ci_low = inv_logit(quantile(ch$beta_ih, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_ih, 0.975)),
                                param = "Intra-household - symptomatic",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_ih + ch$beta_asym)),
                                ci_low = inv_logit(quantile(ch$beta_ih + ch$beta_asym, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_ih + ch$beta_asym, 0.975)),
                                param = "Intra-household - asymptomatic",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = inv_logit(mean(ch$beta_eh)),
                                ci_low = inv_logit(quantile(ch$beta_eh, 0.025)),
                                ci_high = inv_logit(quantile(ch$beta_eh, 0.975)),
                                param = "Extra-household",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = exp(mean(ch$beta_asym)),
                                ci_low = exp(quantile(ch$beta_asym, 0.025)),
                                ci_high = exp(quantile(ch$beta_asym, 0.975)),
                                param = "Asym vs. sym OR",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = mean(ch$p_symp),
                                ci_low = quantile(ch$p_symp, 0.025),
                                ci_high = quantile(ch$p_symp, 0.975),
                                param = "Symptomatic proportion",
                                sim_num = i))
    
    out <- bind_rows(out,
                     data.frame(est = mean(1/ch$gamma),
                                ci_low = quantile(1/ch$gamma, 0.025),
                                ci_high = quantile(1/ch$gamma, 0.975),
                                param = "Duration of infectiousness",
                                sim_num = i))
  }
  return(out)
}



