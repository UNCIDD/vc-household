library(tidyverse)

#################### Generate stan data #############################################

dat <- bl %>%
  left_join(fu) %>%
  group_by(household_id) %>%
  mutate(part_id = individual_id - min(individual_id) + 1) %>%
  ungroup() %>%
  mutate(symp = ifelse(watery_diarrhea == 1, 1, 0)) %>%
  arrange(household_id, day, part_id)
  
dat$symp[dat$day == 1 & dat$part_id == 1] <- 1 # all index cases start symptomatic
  
dat_c <- dat %>%
  dplyr::select(household_id, part_id, day, culture_pos) %>%
  drop_na()
dat_c$row_id <- 1:nrow(dat_c)
  
dat_s <- dat %>%
  dplyr::select(household_id, part_id, day, symp) %>%
  drop_na()
  
# Add dummy row for household 350 to avoid indexing issues
dat_s <- dat_s %>%
  bind_rows(data.frame(household_id = 350,
                       part_id = 99,
                       day = 99,
                       symp = 99)) %>%
  arrange(household_id, day, part_id)
dat_s$row_id <- 1:nrow(dat_s)

dat_v <- dat %>%
  dplyr::select(household_id, part_id, day, vib_titer) %>%
  mutate(vib_titer1280 = ifelse(vib_titer < 1280, 1, 2),
         vib_titer640 = ifelse(vib_titer < 640, 1, 2),
         vib_titer320 = ifelse(vib_titer < 320, 1, 2),
         vib_titer160 = ifelse(vib_titer < 160, 1, 2)) %>%
  drop_na()
  
# Add dummy row for household without vib observations to avoid indexing issues
no_vib_hh <- unique(dat$household_id[which(!(dat$household_id %in% dat_v$household_id))])
dummy_v <- data.frame(household_id = no_vib_hh,
                      part_id = 99,
                      day = 99,
                      vib_titer160 = 99,
                      vib_titer320 = 99,
                      vib_titer640 = 99,
                      vib_titer1280 = 99)
dat_v <- dat_v %>%
  bind_rows(dummy_v) %>%
  arrange(household_id, day, part_id)
dat_v$row_id <- 1:nrow(dat_v)
  
  
hh_sum <- dat %>%
  group_by(household_id) %>%
  summarize(hh_size = household_size[1],
            n_enroll = max(part_id))
  
sum_c <- dat_c %>%
  group_by(household_id) %>%
  summarize(obs_per_hh = n(),
            start_ind = min(row_id),
            end_ind = max(row_id))
  
sum_s <- dat_s %>%
  group_by(household_id) %>%
  summarize(obs_per_hh = n(),
            start_ind = min(row_id),
            end_ind = max(row_id))
  
sum_v <- dat_v %>%
  group_by(household_id) %>%
  summarize(obs_per_hh = n(),
            start_ind = min(row_id),
            end_ind = max(row_id))
  
phi <- 0.13
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
  
dat_stan <- list(n_hh =  nrow(hh_sum),
                 hh_size = hh_sum$hh_size,
                 n_enroll = hh_sum$n_enroll,
                 n_obs_c = nrow(dat_c),
                 n_obs_s = nrow(dat_s),
                 n_obs_v = nrow(dat_v),
                 n_unique_obs_c = 2,
                 n_unique_obs_s = 2,
                 n_unique_obs_v = 2,
                 y_c = dat_c$culture_pos+1,
                 y_s = dat_s$symp+1,
                 y_v = dat_v$vib_titer320,
                 part_id_c = dat_c$part_id,
                 part_id_s = dat_s$part_id,
                 part_id_v = dat_v$part_id,
                 t_day_c = dat_c$day + 1,
                 t_day_s = dat_s$day + 1,
                 t_day_v = dat_v$day + 1,
                 obs_per_hh_c = sum_c$obs_per_hh,
                 obs_per_hh_s = sum_s$obs_per_hh,
                 obs_per_hh_v = sum_v$obs_per_hh,
                 hh_start_ind_c = sum_c$start_ind,
                 hh_end_ind_c = sum_c$end_ind,
                 hh_start_ind_s = sum_s$start_ind,
                 hh_end_ind_s = sum_s$end_ind,
                 hh_start_ind_v = sum_v$start_ind,
                 hh_end_ind_v = sum_v$end_ind,
                 obs_prob_c = obs_prob_c,
                 obs_prob_s = obs_prob_s,
                 obs_prob_v = obs_prob_v,
                 init_probs_index = c((0.1-1.5*1e-10)*0.3, 1e-10, 0.9-1.5*1e-10, 1e-10, (0.1-1.5*1e-10)*0.7),
                 p_symp = 0.1,
                 epsilon = 1e-10,
                 n_days = 31,
                 prop_under1 = 0.3)
  
saveRDS(dat_stan, "data/dat_11.26.24_nocov_320.rds")
  
dat_stan_640 <- dat_stan
dat_stan_640$y_v <- dat_v$vib_titer640
saveRDS(dat_stan_640, "data/dat_11.26.24_nocov_640.rds")
  
dat_stan_160 <- dat_stan
dat_stan_160$y_v <- dat_v$vib_titer160
saveRDS(dat_stan_160, "data/dat_11.26.24_nocov_160.rds")
  
dat_stan_1280 <- dat_stan
dat_stan_1280$y_v <- dat_v$vib_titer1280
saveRDS(dat_stan_1280, "data/dat_11.26.24_nocov_1280.rds")
  
  
# Data with expanded covariates
# Set up data where unenrolled is it's own category
bl_stan <- bl %>%
  arrange(household_id, individual_id) %>%
  dplyr::select(household_id, individual_id, household_size, age_group, sex, soap_in_home, drinking_water,
                monthly_income, relation_to_index, occupat_group) %>%
  group_by(household_id) %>%
  mutate(n_enroll = n(),
         part_id = individual_id - min(individual_id) + 1) %>%
  ungroup() %>%
  mutate(int = 1,
         index = ifelse(part_id == 1, 1, 0),
         age_0.5 = ifelse(age_group == "0-5", 1, 0),
         age_5.17 = ifelse(age_group == "5-17", 1, 0),
         female = ifelse(sex == "female", 1, 0),
         no_soap = ifelse(soap_in_home == "No", 1, 0),
         age_unk = 0,
         sex_unk = 0,
         water_boiled = ifelse(drinking_water == "Boiled water", 1, 0),
         water_private = ifelse(drinking_water == "Private tap", 1, 0),
         water_tubewell = ifelse(drinking_water == "Tubewell", 1, 0),
         water_other = ifelse(drinking_water %in% c("Other", "Public Tap"), 1, 0),
         inc_8500.11999 = ifelse(monthly_income >= 8500 & monthly_income < 12000, 1, 0),
         inc_12000.17999 = ifelse(monthly_income >= 12000 & monthly_income < 18000, 1, 0),
         inc_18000. = ifelse(monthly_income >= 18000, 1, 0),
         occ_home = ifelse(occupat_group == "AtHome", 1, 0),
         occ_under18 = ifelse(occupat_group == "Under18", 1, 0),
         occ_unk = 0,
         rel_child = ifelse(relation_to_index == "Child", 1, 0),
         rel_sibling = ifelse(relation_to_index == "Sibling", 1, 0),
         rel_spouse = ifelse(relation_to_index == "Spouse", 1, 0),
         rel_unk_self = ifelse(relation_to_index == "Unknown" | relation_to_index == "Self", 1, 0)) %>%
  dplyr::select(-age_group, -sex, -soap_in_home, -drinking_water, -relation_to_index, -occupat_group, -monthly_income)
  
# Fill in covariates for unenrolled household members
bl_stan_aug <- data.frame()
hh_ids <- unique(bl_stan$household_id)
for(i in 1:length(hh_ids)) {
  enrolled <- bl_stan %>% filter(household_id == hh_ids[i])
  if(enrolled$household_size[1] != enrolled$n_enroll[1]) {
    enrolled <- enrolled %>%
      bind_rows(data.frame(household_id = hh_ids[i],
                           part_id = (enrolled$n_enroll[1]+1):(enrolled$household_size[1]),
                           household_size = enrolled$household_size[1],
                           n_enroll = enrolled$n_enroll[1],
                           int = 1,
                           index = 0,
                           age_0.5 = 0,
                           age_5.17 = 0,
                           age_unk = 1,
                           female = 0,
                           sex_unk = 1,
                           no_soap = enrolled$no_soap[1],
                           water_boiled = enrolled$water_boiled[1],
                           water_private = enrolled$water_private[1],
                           water_tubewell = enrolled$water_tubewell[1],
                           water_other = enrolled$water_other[1],
                           inc_8500.11999 = enrolled$inc_8500.11999[1],
                           inc_12000.17999 = enrolled$inc_12000.17999[1],
                           inc_18000. = enrolled$inc_18000.[1],
                           occ_home = 0,
                           occ_under18 = 0,
                           occ_unk = 1,
                           rel_child = 0,
                           rel_sibling = 0,
                           rel_spouse = 0,
                           rel_unk_self = 1))
  }
  bl_stan_aug = bind_rows(bl_stan_aug, enrolled)
}
  
dat_stan <- readRDS("data/dat_11.26.24_nocov_320.rds")
dat_stan$x_ih <- bl_stan_aug %>% dplyr::select(int, age_0.5, age_5.17, age_unk, female, sex_unk, water_private, water_tubewell, water_other, no_soap,
                                               inc_8500.11999, inc_12000.17999, inc_18000., occ_home, occ_under18, occ_unk,
                                               rel_child, rel_sibling, rel_spouse, rel_unk_self)
dat_stan$x_eh <- bl_stan_aug %>% dplyr::select(int, age_0.5, age_5.17, age_unk, female, sex_unk, water_private, water_tubewell, water_other, no_soap,
                                               inc_8500.11999, inc_12000.17999, inc_18000., occ_home, occ_under18, occ_unk,
                                               rel_child, rel_sibling, rel_spouse, rel_unk_self)
dat_stan$k_ih <- ncol(dat_stan$x_ih)
dat_stan$k_eh <- ncol(dat_stan$x_eh)
saveRDS(dat_stan, "data/dat_11.26.24_cov_320_new.rds")
