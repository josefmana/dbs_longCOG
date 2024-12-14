#
# In this script I fit and compare models used to provide some preliminary answers to oponents' questions
# in the defense of the thesis.
#


library(here)
library(tidyverse)
library(brms)
library(performance)
library(tidybayes)

# read data
data <- readRDS( here("_data","cogPRED_data.rds") )
df <- data$d0

# scale variables of interest
df$drs <- with( data$scl, (df$drs_tot - M$drs) / SD$drs )
df$time <- df$time_y + data$scl$Md$time
df$year <- as.numeric( scale(df$stim_year, center = T, scale = T) )
df$cohort <- factor(df$stim_year, ordered = T)

# prepare linear models
lm <- list(
  
  fixed_slope = bf( drs ~ 1 + time + (1 | id) ),
  mixed_slope = bf( drs ~ 1 + time + (1 + time | id) ),
  scaleloc1 = bf( drs ~ 1 + time + (1 | id), sigma ~ 1 + (1 | id) ),
  scaleloc2 = bf( drs ~ 1 + time + (1 | id), sigma ~ 1 + time + (1 | id) ),
  cohort_slope = bf( drs ~ 1 + time * year + (1 + time | id) ),
  cohort_nested = bf( drs ~ 1 + time + (1 + time | cohort / id ) )
  
)

# put regularizing priors on slopes
p <- prior("normal(0,0.5)", class = "b")

# fit them
fit <- lapply(
  
  set_names( names(lm) ),
  function(x) brm(
    
    formula = lm[[x]],
    data = df,
    family = gaussian(),
    prior = p,
    cores = parallel::detectCores(),
    save_pars = save_pars(all = T),
    seed = 87542
    
  )
  
)

# print ICC estimates for fixed scale models
variance_decomposition(fit$fixed_slope)
variance_decomposition(fit$mixed_slope)

# check predictive performance of the 'cohort effect' compared to 'vanilla' LMM
loo(fit$mixed_slope, fit$cohort_slope, fit$cohort_nested, moment_match = T)

# plot 'cohort effects' from the nested model
fit$cohort_nested %>%
  
  spread_draws(r_cohort[year, parameter]) %>%
  mutate(
    deviation = r_cohort * data$scl$SD$drs,
    parameter = if_else(parameter == "time", "Slope", parameter)
  ) %>%
  median_qi( deviation, .width = c(.95, .66) ) %>%
  
  ggplot() +
  aes(y = deviation, x = year, ymin = .lower, ymax = .upper) +
  geom_pointinterval() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  scale_x_continuous(
    breaks = sort( unique(df$stim_year) ),
    labels = sort( unique(df$stim_year) ),
    limits = c( min(df$stim_year), max(df$stim_year) )
  ) +
  labs(y = "Deviation from mean", x = "Year of surgery") +
  facet_wrap( ~ parameter, nrow = 2, scales = "free_y") +
  theme_bw(base_size = 14)

# save the figure
ggsave(
  plot = last_plot(),
  filename = here("_figs","_cohort_effect.jpg"),
  dpi = 300,
  width = 10.4,
  height = 4.9
)

