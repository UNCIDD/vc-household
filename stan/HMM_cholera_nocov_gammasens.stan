

data {
  
  // Household information
  int n_hh; // number of households
  int hh_size[n_hh]; // household size
  int n_days;
    
  // Household-level data
  int n_obs_c; // number of culture observations
  int n_obs_s; // number of symptom observations
  int n_obs_v; // number of vibriocidal observations
  int n_unique_obs_c; // number of unique outcomes for culture
  int n_unique_obs_s; // number of unique outcomes for symptoms
  int n_unique_obs_v; // number of unique outcomes for vib
  int y_c[n_obs_c]; // culture outcome, ordered by 1) household, then 2) time, then 3) individual
  int y_s[n_obs_s]; // symptom outcome, ordered by 1) household, then 2) time, then 3) individual
  int y_v[n_obs_v]; // vib outcome, ordered by 1) household, then 2) time, then 3) individual
  int part_id_c[n_obs_c]; // particpants associated with the culture observation
  int part_id_s[n_obs_s]; // particpants associated with the symptom observation
  int part_id_v[n_obs_v]; // particpants associated with the vib observation
  int t_day_c[n_obs_c]; // observation times for culture (day)
  int t_day_s[n_obs_s]; // observation times for symptoms (day)
  int t_day_v[n_obs_v]; // observation times for vib (day)
  int obs_per_hh_c[n_hh]; // total number of culture observations per HH
  int obs_per_hh_s[n_hh]; // total number of symptom observations per HH
  int obs_per_hh_v[n_hh]; // total number of vib observations per HH
  int hh_start_ind_c[n_hh]; // starting index for the HH - culture observations
  int hh_end_ind_c[n_hh]; // ending index for the HH - culture observations
  int hh_start_ind_s[n_hh]; // starting index for the HH - symptom observations
  int hh_end_ind_s[n_hh]; // ending index for the HH - symptom observations
  int hh_start_ind_v[n_hh]; // starting index for the HH - vib observations
  int hh_end_ind_v[n_hh]; // ending index for the HH - vib observations
  
  // Observation probabilities  
  matrix[n_unique_obs_c, 5] obs_prob_c; // observation process for culture
  matrix[n_unique_obs_s, 5] obs_prob_s; // observation process symptoms
  matrix[n_unique_obs_v, 5] obs_prob_v; // observation process vib
  simplex[5] init_probs_index; // starting state probabilities forindex household member
  
  real epsilon; // very small number to avoid log(0) issues
  real<lower=0, upper=1> prop_under1; // Proportion of cases with incubation period less than one day
  real gamma_mult;
  
}

parameters {
  real beta_asym;
  real beta_ih;
  real beta_eh;
  real logit_p_symp; // proportion of infections that are symptomatic
  real logit_gamma_sym; // recovery rate
  real logit_sigma; // E > I transition rate
  simplex[5] init_probs_hh; // starting state probabilities for non-index household member
}

transformed parameters {
  matrix[sum(hh_size), n_days] llik; // lik contribution for enrolled per participant and time
  matrix[sum(hh_size)*5, n_days] logalpha; // log forward probability
  matrix[sum(hh_size)*5, n_days] alpha; // forward prob, normalized
  real<lower=0, upper = 1> eh_prob;
  real<lower=0, upper = 1> ih_prob_sym;
  real<lower=0, upper = 1> ih_prob_asym;
  real gamma_sym;
  real gamma_asym;
  real p_symp;
  real sigma;
  
  eh_prob = inv_logit(beta_eh);
  ih_prob_sym = inv_logit(beta_ih);
  ih_prob_asym = inv_logit(beta_ih + beta_asym);
  gamma_sym = inv_logit(logit_gamma_sym);
  gamma_asym = gamma_sym*gamma_mult;
  p_symp = inv_logit(logit_p_symp);
  sigma = inv_logit(logit_sigma);

  llik = rep_matrix(0, sum(hh_size), n_days);
  
  for(h in 1:n_hh) { // loop through household
    
    int y_c_hh[obs_per_hh_c[h]];
    int y_s_hh[obs_per_hh_s[h]];
    int y_v_hh[obs_per_hh_v[h]];
    int part_id_c_hh[obs_per_hh_c[h]];
    int part_id_s_hh[obs_per_hh_s[h]];
    int part_id_v_hh[obs_per_hh_v[h]];
    int t_day_c_hh[obs_per_hh_c[h]];
    int t_day_s_hh[obs_per_hh_s[h]];
    int t_day_v_hh[obs_per_hh_v[h]];
    
    int index_c; // index for next observation - cultuure
    int index_s; // index for next observation - symptoms
    int index_v; // index for next observation - vib
    int sym_rows[hh_size[h]]; // rows in alpha corresponding to I
    int asym_rows[hh_size[h]]; // rows in alpha corresponding to I
    int last_lik;
    int obs_switch_c; // indicatior if there's an observation for this time step and participant - culture
    int obs_switch_s; // indicatior if there's an observation for this time step and participant - symptoms
    int obs_switch_v; // indicatior if there's an observation for this time step and participant - vib
    real obs_c[5]; // observation probabilities - culture
    real obs_s[5]; // observation probabilities - culture
    real obs_v[5]; // observation probabilities - culture
    vector[hh_size[h]] no_hh_inf; // probability of avoiding infection from within the household
    real no_inf; // probability of avoiding all infection
    vector[5] init_probs;
    
    if(h == 1) {
      last_lik = 0;
    } else {
      last_lik = sum(hh_size[1:(h-1)]);
    }
    
    // subset to data only for the given HH
    y_c_hh = y_c[(hh_start_ind_c[h]):(hh_end_ind_c[h])];
    y_s_hh = y_s[(hh_start_ind_s[h]):(hh_end_ind_s[h])];
    y_v_hh = y_v[(hh_start_ind_v[h]):(hh_end_ind_v[h])];
    t_day_c_hh = t_day_c[(hh_start_ind_c[h]):(hh_end_ind_c[h])];
    t_day_s_hh = t_day_s[(hh_start_ind_s[h]):(hh_end_ind_s[h])];
    t_day_v_hh = t_day_v[(hh_start_ind_v[h]):(hh_end_ind_v[h])];
    part_id_c_hh = part_id_c[(hh_start_ind_c[h]):(hh_end_ind_c[h])];
    part_id_s_hh = part_id_s[(hh_start_ind_s[h]):(hh_end_ind_s[h])];
    part_id_v_hh = part_id_v[(hh_start_ind_v[h]):(hh_end_ind_v[h])];
    
    index_c = 1;
    index_s = 1;
    index_v = 1;
    
    { // START FORWARD ALGORITHM
    // fill first column of alpha using starting probabilities
    
    for(i in 1:hh_size[h]) {
        if(i == 1) {
          init_probs = init_probs_index;
        } else {
          init_probs = init_probs_hh;
        }
      
      obs_switch_c = 0;
      obs_switch_s = 0;
      obs_switch_v = 0;
      
      if(t_day_c_hh[index_c] == 1) {
        if(part_id_c_hh[index_c] == i) {
          obs_switch_c = 1;
        }
      }
      
      if(t_day_s_hh[index_s] == 1) {
        if(part_id_s_hh[index_s] == i) {
          obs_switch_s = 1;
        }
      }
      
      if(t_day_v_hh[index_v] == 1) {
        if(part_id_v_hh[index_v] == i) {
          obs_switch_v = 1;
        }
      }
      
      if(obs_switch_c == 1) {
        obs_c[1] = obs_prob_c[y_c_hh[index_c], 1]; // Pr(y_t | S)
        obs_c[2] = obs_prob_c[y_c_hh[index_c], 2]; // Pr(y_t | Ia)
        obs_c[3] = obs_prob_c[y_c_hh[index_c], 3]; // Pr(y_t | Is) 
        obs_c[4] = obs_prob_c[y_c_hh[index_c], 4]; // Pr(y_t | R)
        obs_c[5] = obs_prob_c[y_c_hh[index_c], 5]; // Pr(y_t | E)
      } else {
        obs_c = rep_array(1, 5);
      }
        
      if(obs_switch_s == 1) {
        obs_s[1] = obs_prob_s[y_s_hh[index_s], 1]; // Pr(y_t | S)
        obs_s[2] = obs_prob_s[y_s_hh[index_s], 2]; // Pr(y_t | Ia)
        obs_s[3] = obs_prob_s[y_s_hh[index_s], 3]; // Pr(y_t | Is) 
        obs_s[4] = obs_prob_s[y_s_hh[index_s], 4]; // Pr(y_t | R) 
        obs_s[5] = obs_prob_s[y_s_hh[index_s], 5]; // Pr(y_t | E) 
      } else {
        obs_s = rep_array(1, 5);
      }
        
      if(obs_switch_v == 1) {
        obs_v[1] = obs_prob_v[y_v_hh[index_v], 1]; // Pr(y_t | S)
        obs_v[2] = obs_prob_v[y_v_hh[index_v], 2]; // Pr(y_t | Ia)
        obs_v[3] = obs_prob_v[y_v_hh[index_v], 3]; // Pr(y_t | Is) 
        obs_v[4] = obs_prob_v[y_v_hh[index_v], 4]; // Pr(y_t | R) 
        obs_v[5] = obs_prob_v[y_v_hh[index_v], 5]; // Pr(y_t | E) 
       } else {
        obs_v = rep_array(1, 5);
       }   
      
      logalpha[5*last_lik+5*(i-1)+1, 1] = log(init_probs[1]) + log(obs_c[1]) + log(obs_s[1]) + log(obs_v[1]);
      logalpha[5*last_lik+5*(i-1)+2, 1] = log(init_probs[2]) + log(obs_c[2]) + log(obs_s[2]) + log(obs_v[2]);
      logalpha[5*last_lik+5*(i-1)+3, 1] = log(init_probs[3]) + log(obs_c[3]) + log(obs_s[3]) + log(obs_v[3]);
      logalpha[5*last_lik+5*(i-1)+4, 1] = log(init_probs[4]) + log(obs_c[4]) + log(obs_s[4]) + log(obs_v[4]);
      logalpha[5*last_lik+5*(i-1)+5, 1] = log(init_probs[5]) + log(obs_c[5]) + log(obs_s[5]) + log(obs_v[5]);
      asym_rows[i] = 5*last_lik+5*(i-1)+2;
      sym_rows[i] = 5*last_lik+5*(i-1)+3;
      
      llik[last_lik + i, 1] = log_sum_exp(logalpha[(5*last_lik+5*(i-1)+1):(5*last_lik+5*(i-1)+5),1]);
      
      // normalize and convert to the probability scale
      alpha[(5*last_lik+5*(i-1)+1):(5*last_lik+5*(i-1)+5), 1] = softmax(logalpha[(5*last_lik+5*(i-1)+1):(5*last_lik+5*(i-1)+5),1]);
    
      if(obs_switch_c == 1) {
        index_c = min(index_c + 1, hh_end_ind_c[h]-hh_start_ind_c[h]+1);   
      }
      if(obs_switch_s == 1) {
        index_s = min(index_s + 1, hh_end_ind_s[h]-hh_start_ind_s[h]+1);   
      }
      if(obs_switch_v == 1) {
        index_v = min(index_v + 1, hh_end_ind_v[h]-hh_start_ind_v[h]+1);   
      }
    }
        
    for (tt in 2:n_days) {
      for(p in 1:hh_size[h]) {
        real no_inf_prob; // probability of avoiding all infections
        vector[hh_size[h]] no_hh_inf_prob; // probability of avoiding infection from each HH member
        
        obs_switch_c = 0;
        obs_switch_s = 0;
        obs_switch_v = 0;
        
        if(t_day_c_hh[index_c] == tt) {
          if(part_id_c_hh[index_c] == p) {
            obs_switch_c = 1;
          }
        }
        
        if(t_day_s_hh[index_s] == tt) {
          if(part_id_s_hh[index_s] == p) {
            obs_switch_s = 1;
          }
        }
        
        if(t_day_v_hh[index_v] == tt) {
          if(part_id_v_hh[index_v] == p) {
            obs_switch_v = 1;
          }
        }
        
        // Observation probabilities  
      if(obs_switch_c == 1) {
        obs_c[1] = obs_prob_c[y_c_hh[index_c], 1]; // Pr(y_t | S)
        obs_c[2] = obs_prob_c[y_c_hh[index_c], 2]; // Pr(y_t | Ia)
        obs_c[3] = obs_prob_c[y_c_hh[index_c], 3]; // Pr(y_t | Is) 
        obs_c[4] = obs_prob_c[y_c_hh[index_c], 4]; // Pr(y_t | R)
        obs_c[5] = obs_prob_c[y_c_hh[index_c], 5]; // Pr(y_t | E)
      } else {
        obs_c = rep_array(1, 5);
      }
        
      if(obs_switch_s == 1) {
        obs_s[1] = obs_prob_s[y_s_hh[index_s], 1]; // Pr(y_t | S)
        obs_s[2] = obs_prob_s[y_s_hh[index_s], 2]; // Pr(y_t | Ia)
        obs_s[3] = obs_prob_s[y_s_hh[index_s], 3]; // Pr(y_t | Is) 
        obs_s[4] = obs_prob_s[y_s_hh[index_s], 4]; // Pr(y_t | R)
        obs_s[5] = obs_prob_s[y_s_hh[index_s], 5]; // Pr(y_t | E)
      } else {
        obs_s = rep_array(1, 5);
      }
      
      if(obs_switch_v == 1) {
        obs_v[1] = obs_prob_v[y_v_hh[index_v], 1]; // Pr(y_t | S)
        obs_v[2] = obs_prob_v[y_v_hh[index_v], 2]; // Pr(y_t | Ia)
        obs_v[3] = obs_prob_v[y_v_hh[index_v], 3]; // Pr(y_t | Is) 
        obs_v[4] = obs_prob_v[y_v_hh[index_v], 4]; // Pr(y_t | R)
        obs_v[5] = obs_prob_v[y_v_hh[index_v], 5]; // Pr(y_t | E)
      } else {
        obs_v = rep_array(1, 5);
      }
        
      no_hh_inf = to_vector(alpha[asym_rows, tt-1]).*rep_vector(1-ih_prob_asym, hh_size[h]) +
                  to_vector(alpha[sym_rows, tt-1]).*rep_vector(1-ih_prob_sym, hh_size[h]) + 
                  (rep_vector(1, hh_size[h]) - to_vector(alpha[asym_rows, tt-1]) - to_vector(alpha[sym_rows, tt-1]));  
      no_hh_inf[p] = 1; // Particpant can't infect themselves
      no_inf = (1-eh_prob)*prod(no_hh_inf); // Pr avoiding all infections
        
        // Compute the probability of each state given prior probabilities and observation  
        logalpha[5*last_lik+5*(p-1)+1, tt] = log(exp(logalpha[5*last_lik+5*(p-1)+1, tt-1])*(no_inf-epsilon/4) + //Pr(S>>S)
                                      exp(logalpha[5*last_lik+5*(p-1)+2, tt-1])*epsilon + //Pr(Ia>>S)
                                      exp(logalpha[5*last_lik+5*(p-1)+3, tt-1])*epsilon + //Pr(Is>>S)
                                      exp(logalpha[5*last_lik+5*(p-1)+4, tt-1])*epsilon + //Pr(R>>S)
                                      exp(logalpha[5*last_lik+5*(p-1)+5, tt-1])*epsilon) + //Pr(E>>S)
                                      log(obs_c[1]) + log(obs_s[1]) + log(obs_v[1]);
        logalpha[5*last_lik+5*(p-1)+2, tt] = log(exp(logalpha[5*last_lik+5*(p-1)+1, tt-1])*(prop_under1*(1-p_symp)*(1-no_inf)-epsilon/4) + // Pr (S>>Ia)
                                      exp(logalpha[5*last_lik+5*(p-1)+2, tt-1])*(1-gamma_asym-3*epsilon/2) + // Pr (Ia>>Ia)
                                      exp(logalpha[5*last_lik+5*(p-1)+3, tt-1])*epsilon + // Pr (Is>>Ia)
                                      exp(logalpha[5*last_lik+5*(p-1)+4, tt-1])*epsilon + // Pr (R>>Ia)
                                      exp(logalpha[5*last_lik+5*(p-1)+5, tt-1])*(sigma*(1-p_symp)-2*epsilon/3)) + // Pr (E>>Ia)
                                      log(obs_c[2]) + log(obs_s[2]) + log(obs_v[2]);
        logalpha[5*last_lik+5*(p-1)+3, tt] = log(exp(logalpha[5*last_lik+5*(p-1)+1, tt-1])*(prop_under1*p_symp*(1-no_inf)-epsilon/4) + // Pr (S>>Is)
                                      exp(logalpha[5*last_lik+5*(p-1)+2, tt-1])*epsilon + // Pr (Ia>>Is)
                                      exp(logalpha[5*last_lik+5*(p-1)+3, tt-1])*(1-gamma_sym-3*epsilon/2) + // Pr (Is>>Is)
                                      exp(logalpha[5*last_lik+5*(p-1)+4, tt-1])*epsilon + // Pr (R>>Is)
                                      exp(logalpha[5*last_lik+5*(p-1)+5, tt-1])*(sigma*p_symp-2*epsilon/3)) + // Pr (E>>Is)
                                      log(obs_c[3]) + log(obs_s[3]) + log(obs_v[3]);
        logalpha[5*last_lik+5*(p-1)+4, tt] = log(exp(logalpha[5*last_lik+5*(p-1)+1, tt-1])*epsilon + // Pr (S>>R)
                                      exp(logalpha[5*last_lik+5*(p-1)+2, tt-1])*(gamma_asym-3*epsilon/2) + // Pr (Ia>>R)
                                      exp(logalpha[5*last_lik+5*(p-1)+3, tt-1])*(gamma_sym-3*epsilon/2) + // Pr (Is>>R)
                                      exp(logalpha[5*last_lik+5*(p-1)+4, tt-1])*(1-4*epsilon) + // Pr (R>>R)
                                      exp(logalpha[5*last_lik+5*(p-1)+5, tt-1])*epsilon) + // Pr (E>>R)
                                      log(obs_c[4]) + log(obs_s[4]) + log(obs_v[4]);
        logalpha[5*last_lik+5*(p-1)+5, tt] = log(exp(logalpha[5*last_lik+5*(p-1)+1, tt-1])*((1-no_inf)*(1-prop_under1)-epsilon/4) + // Pr (S>>E)
                                      exp(logalpha[5*last_lik+5*(p-1)+2, tt-1])*epsilon + // Pr (Ia>>E)
                                      exp(logalpha[5*last_lik+5*(p-1)+3, tt-1])*epsilon + // Pr (Is>>E)
                                      exp(logalpha[5*last_lik+5*(p-1)+4, tt-1])*epsilon + // Pr (R>>E)
                                      exp(logalpha[5*last_lik+5*(p-1)+5, tt-1])*(1-sigma-2*epsilon/3)) + // Pr (E>>E)
                                      log(obs_c[5]) + log(obs_s[5]) + log(obs_v[5]);
        
        // normalize and convert to probability scale
        alpha[(5*last_lik+5*(p-1)+1):(5*last_lik+5*(p-1)+5), tt] = softmax(logalpha[(5*last_lik+5*(p-1)+1):(5*last_lik+5*(p-1)+5),tt]);
        
        llik[last_lik + p, tt] = log_sum_exp(logalpha[(5*last_lik+5*(p-1)+1):(5*last_lik+5*(p-1)+5),tt]);
        
        if(obs_switch_c == 1) {
          index_c = min(index_c + 1, hh_end_ind_c[h]-hh_start_ind_c[h]+1);  
        }
        if(obs_switch_s == 1) {
          index_s = min(index_s + 1, hh_end_ind_s[h]-hh_start_ind_s[h]+1);  
        }
        if(obs_switch_v == 1) {
          index_v = min(index_v + 1, hh_end_ind_v[h]-hh_start_ind_v[h]+1);  
        }
        
      } // end participant loop
        
    } // end time loop
    } // END FORWARD ALGORITHM
    
  } // end household loop
  
}


model {
  
  // beta coefficients
  beta_eh ~ normal(-3,3);
  beta_ih ~ normal(-3,3);
  beta_asym ~ normal(0,1);
  init_probs_hh ~  dirichlet(rep_vector(0.5, 5));
  logit_sigma ~ normal(0.9, 0.5);
  logit_gamma_sym ~ normal(0, 0.5);
  
  // Only increment by final alpha
  target += sum(llik[,n_days]);

}

