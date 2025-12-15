make_fig1 <- function() {
  
  # Data frame with observation timing
  plot <- bind_rows(data.frame(day = c(1:10,30),
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
                               obs_type = "Culture"))
  
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
