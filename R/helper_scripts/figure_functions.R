
sample_timing <- function() {
  return(bind_rows(data.frame(day = c(1:10,30),
                              part_type = "Index",
                              obs_type = "Symptoms"),
                   data.frame(day = c(2:10,30),
                              part_type = "Household contact",
                              obs_type = "Symptoms"),
                   data.frame(day = c(2, 7,30),
                              part_type = "Index",
                              obs_type = "Vibriocidal titer"),
                   data.frame(day = c(2, 7, 30),
                              part_type = "Household contact",
                              obs_type = "Vibriocidal titer"),
                   data.frame(day = 1,
                              part_type = "Index",
                              obs_type = "Culture"),
                   data.frame(day = c(2:10),
                              part_type = "Household contact",
                              obs_type = "Culture")))
}

make_fig1 <- function() {
  
  # Data frame with observation timing
  plot <- sample_timing()
  
  ggplot(plot) +
    geom_point(aes(x = factor(day, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "", "30")),
                   y = obs_type, color = part_type, shape = part_type), size = 2)  +
    scale_color_brewer(palette = "Set1") +
    scale_x_discrete(drop = FALSE,
                     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "", "30")) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    labs(x = "Day of follow-up", y = "")
  
}

make_fig3 <- function(dat) {
  
  ggplot(dat %>% filter(str_detect(param, "household")),
         aes(x = param, y = est_med, ymin = ci_low, ymax = ci_high)) +
    geom_point() +
    # geom_point(aes(y = est), shape = 8) +
    geom_errorbar(width = 0.5) +
    theme_bw() +
    labs(x = "", y = "Per day probability of infection") +
    theme(axis.text.x = element_text(size = 11),)
  
}

make_fig4 <- function(dat) {
  
  dat <- dat %>%
    mutate(cov = factor(cov, levels = c("Intercept", "Unknown",
                                        "0-4", "5-17", "18+",
                                        "Male", "Female",
                                        "0-8,499", "8,500-11,999", "12,000-17,999", "18,000+",
                                        "Under 18", "At home", "Outside home",
                                        "Unknown or self", "Spouse", "Sibling", "Child", "Parent",
                                        "Soap", "No soap")),
           cov_group = factor(cov_group, levels = c("Sex", "Water\nsource", "Soap\nin home", "Monthly\nincome",
                                                    "Occupation\ngroup", "Relation\nto index", "Age"))) 
  
  ggplot(dat %>% filter(cov != "Intercept"),
         aes(color = model_type)) +
    geom_point(aes(y = cov, x = Est_med), position=position_dodge(width=0.5)) +
    geom_errorbar(aes(y = cov, xmin = CI_low, xmax = CI_high),
                  position=position_dodge(width=0.5), width = 0.25) +
    facet_grid(cols = vars(param), rows = vars(cov_group), scales = "free_y") +
    theme_bw() +
    scale_x_log10() +
    labs(x = "Odds ratio",
         y = "") +
    geom_vline(xintercept = 1, lty = "dashed") +
    theme(legend.title = element_blank())
  
}

make_table1 <- function(bl, fu) {
  
  hh_info <- bl %>%
    left_join(fu) %>%
    group_by(household_id) %>%
    summarize(household_size = household_size[1],
              monthly_income = monthly_income[1],
              drinking_water = drinking_water[1],
              soap_in_home = soap_in_home[1]) %>%
    mutate(hh_size_cat = ifelse(household_size <= 4, "2-4",
                                ifelse(household_size <= 7, "5-7",
                                       ifelse(household_size <= 10, "8-10", "11-14"))),
           income_cat = ifelse(monthly_income <= 4999, "0-4,999",
                               ifelse(monthly_income <= 9999, "5,000-9,999",
                                      ifelse(monthly_income <= 14999, "10,000-14,999",
                                             ifelse(monthly_income <= 19999, "15,000-19,999", "20,000+")))))
  
  out <- bind_rows(as.data.frame(table(hh_info$hh_size_cat)) %>% mutate(variable = "Household size"),
                   as.data.frame(table(hh_info$income_cat)) %>% mutate(variable = "Monthly income (BDT)"),
                   as.data.frame(table(hh_info$drinking_water)) %>% mutate(variable = "Drinking water source"),
                   as.data.frame(table(hh_info$soap_in_home)) %>% mutate(variable = "Use of soap in home")) %>%
    rename(level = Var1,
           count = Freq) %>%
    relocate(variable)
  
  return(out)

}

make_table2 <- function(bl, fu) {
  
  dat <- bl %>% 
    left_join(fu) %>%
    mutate(part_type = ifelse(household_contact == 0, "Index", "Household contact")) %>%
    group_by(individual_id, part_type) %>%
    summarize(wd = ifelse(sum(is.na(watery_diarrhea)) == n(), NA, max(watery_diarrhea[day <= 10], na.rm = T)),
              vom = ifelse(sum(is.na(vomiting)) == n(), NA, max(vomiting[day <= 10], na.rm = T)),
              culture_res = ifelse(sum(is.na(culture_pos)) == n(), NA, ifelse(sum(culture_pos == 1, na.rm = T) >= 1, 1, 0)),
              age = age[1],
              sex = sex[1],
              occ_cat = occupat_group[1],
              relation_to_index = relation_to_index[1],
              income_cat = ifelse(monthly_income[1] <= 4999, "0-4,999",
                                  ifelse(monthly_income[1] <= 9999, "5,000-9,999",
                                         ifelse(monthly_income[1] <= 14999, "10,000-14,999",
                                                ifelse(monthly_income[1] <= 19999, "15,000-19,999", "20,000+")))),
              water_source = drinking_water[1],
              soap_in_home = soap_in_home[1]) %>%
    mutate(symp_type = ifelse(is.na(vom) & is.na(wd), "Missing",
                              ifelse(vom == 1 & wd == 1, "Watery diarrhea & vominiting",
                                     ifelse(vom == 1, "Vomiting only",
                                            ifelse(wd == 1, "Watery_diarrhea only", "None")))))
  
  sum_table <- function(x, colname) {
    out  <- bind_rows(data.frame(variable = "Age",
                                 Var1 = c("Mean", "IQR_25", "IQR_75"),
                                 Freq = c(mean(x$age, na.rm = T), quantile(x$age, 0.25), quantile(x$age, 0.75))),
                      as.data.frame(table(x$sex)) %>% mutate(variable = "Sex"),
                      as.data.frame(table(x$occ_cat)) %>% mutate(variable = "Occupation group"),
                      as.data.frame(table(x$relation_to_index)) %>% mutate(variable = "Relation to index"),
                      as.data.frame(table(x$income_cat)) %>% mutate(variable = "Monthly household income"),
                      as.data.frame(table(x$water_source)) %>% mutate(variable = "Drinking water source"),
                      as.data.frame(table(x$soap_in_home)) %>% mutate(variable = "Soap in home"),
                      as.data.frame(table(x$symp_type)) %>% mutate(variable = "Symptoms in first 10 days")) %>%
      mutate(Freq = round(Freq, 1)) %>%
      rename(level = Var1,
             !!sym(colname) := Freq)
    
    return(out)
  }
  
  out <- sum_table(dat, "All_participants") %>%
    left_join(sum_table(dat %>% filter(part_type == "Index"), "Index_cases")) %>%
    left_join(sum_table(dat %>% filter(part_type == "Household contact"), "Index_cases")) %>%
    left_join(sum_table(dat %>% filter(part_type == "Household contact", culture_res == 1), "HH_contact_culture_pos")) %>%
    left_join(sum_table(dat %>% filter(part_type == "Household contact", culture_res == 0), "HH_contact_culture_neg")) %>%
    left_join(sum_table(dat %>% filter(part_type == "Household contact", is.na(culture_res)), "HH_contact_culture_missing"))
  
  return(out)
  
}

make_figS1 <- function(bl, fu) {
  
  plot <- sample_timing()
  
  counts <- bl %>%
    left_join(fu) %>%
    pivot_longer(c(vib_titer, culture_pos, watery_diarrhea)) %>%
    mutate(obs_type = ifelse(name == "vib_titer", "Vibriocidal titer",
                             ifelse(name == "culture_pos", "Culture", "Symptoms")),
           part_type = ifelse(household_contact == 0, "Index", "Household contact")) %>%
    group_by(household_contact, day, obs_type) %>%
    summarize(n_obs = sum(!is.na(value))) %>%
    mutate(part_type = ifelse(household_contact == 1, "Household contact", "Index"))
  
  plot <- left_join(plot, counts)
  
  ggplot(plot) +
    geom_bar(aes(x = factor(day), y = n_obs, fill = part_type), stat = "identity") +
    facet_grid(rows = vars(obs_type)) +
    scale_fill_manual(values = c("Household contact" = "green4", "Index" = "purple4")) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    labs(x = "Day", y = "Number of observations")
  
}

make_figS2 <- function(bl, fu) {
  
  plt <- fu %>%
    left_join(bl) %>%
    mutate(index = ifelse(household_contact == 0, "Index case", "Household contact"),
           high_vib = ifelse(vib_titer >= 320, 1, 0),
           symp_cult = ifelse(culture_pos == 1 & watery_diarrhea == 1, 1, 0)) %>%
    pivot_longer(cols = c(high_vib, culture_pos, watery_diarrhea, symp_cult)) %>%
    mutate(res_type = ifelse(name == "high_vib", "High vibriocidal titer",
                             ifelse(name == "culture_pos", "Positive culture",
                                    ifelse(name == "watery_diarrhea", "Symptoms", "Positive culture + symptoms")))) %>%
    group_by(day, res_type, index) %>%
    summarize(prop = mean(value, na.rm = T)) %>%
    filter(!(res_type == "Positive culture" & day == 1))
  
  ggplot() +
    geom_point(data = plt %>% filter(!is.na(prop), res_type != "Positive culture + symptoms"),
               aes(x = factor(day, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "", "30")),
                   y = prop, color=index, shape = index, group = index ),
               size = 2.5) +
    facet_grid(rows = vars(res_type), scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title=element_blank()) +
    labs(x = "Day of follow-up",
         y = "Proportion") +
    scale_color_brewer(palette = "Dark2") +
    scale_x_discrete(drop = FALSE,
                     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "", "30"))
}

make_figS3 <-  function(ih_chains) {
  
  diff <- ih_chains$ih_sym - ih_chains$ih_asym
  
  xx <- seq(0, .075, .001)
  p_diff <- data.frame(xx = xx,
                       pr = numeric(length(xx)))
  for(i in 1:length(xx)) {
    p_diff$pr[i] <- sum(diff > xx[i])/length(diff)
  }
  
  ggplot(p_diff) +
    geom_line(aes(x = xx, y = pr, group = 1)) +
    theme_bw() +
    labs(x = "Difference threshold (symptomatic - asymptomatic)",
         y = "Probability")
  
}
