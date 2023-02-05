# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "dplyr", "tidyverse", # for data wrangling
           "brms", "tidybayes", # for Bayesian analyses 
           "ggplot2", "patchwork" # for plotting
           )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set some values for later
s = 87542 # seed for reproducibility

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2500 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .99 # adapt_delta parameter

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_minimal(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# read the data set
d0 <- read.csv( "data/20220508_dbs_longCOG_data.csv" , sep = "," )
df <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data


# ---- pre-processing  ----

# before merging compute scaling values for DRS-2, BDI-II, LEDD, age and time
scl <- list( M = list( drs = mean( df$drs_tot , na.rm = T ), # 136.90
                       bdi = mean( df$bdi , na.rm = T ), # 10.95
                       led = mean( df$ledd_mg , na.rm = T ), # 1197.21
                       age = mean( df$age_ass_y, na.rm = T ) # 59.63
                       ),
             SD = list( drs = sd( df$drs_tot , na.rm = T ), # 7.72
                        bdi = sd( df$bdi , na.rm = T ), # 7.19
                        led = sd( df$ledd_mg , na.rm = T ), # 687.20
                        age = sd( df$age_ass_y, na.rm = T ) # 8.27
                        ),
             # add median time of pre-surgery assessment for GLMMs intercepts s
             Md = list(
               time = -median( df[df$ass_type == "pre", ]$time_y , na.rm = T ) # 0.30
               )
             )

# scale the data
df <- df %>%
  filter( complete.cases(drs_tot) ) %>% # get rid of three dummy rows due to more than one stimulation parameter (no DRS-2)
  # scale all DRS-2, BDI-II, LEDD and time already
  mutate( time = time_y + scl$Md$time,
          drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
          cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ) # right censoring for DRS == 144
          ) %>%
  # keep only variables of interest
  select( id, time, drs, cens_drs )


# ---- models set-up ---- 

# set-up the linear model
f <- list( m0_linear = bf( drs | cens(cens_drs) ~ time + (1 + time | id) ),
           m0_spline = bf( drs | cens(cens_drs) ~ t2(time) + (1 + time | id) )
           )

# set-up priors (using brms default non-informative to allow for information from data to prevail)
p <- NULL


# ---- models fitting ----

# prepare a list for single models
m <- list()

# conduct model fitting
for ( i in names(f) ) m[[i]] <- brm( formula = f[[i]], family = student(), prior = p,
                                     data = df, sample_prior = T, seed = s, chains = ch,
                                     iter = it, warmup = wu, control = list( adapt_delta = ad ),
                                     file = paste0( "models/",i,".rds" ), save_model = paste0("models/",i,".stan")
                                     )


# ---- soft model checking ----

# add PSIS-LOO to both model for influential variables check and model comparisons
for ( i in names(m) ) m[[i]] <- add_criterion( m[[i]] , criterion = c("loo","waic") )

# check the highest Rhat for chains convergence and Pareto-k for influential outliers
cbind.data.frame(
  R_hat = sapply( names(m) , function(i) max( rhat(m[[i]]), na.rm = T ) ) %>% round(3),
  Pareto_k = sapply( names(m) , function(i) max( loo(m[[i]])$diagnostics$pareto_k, na.rm = T ) ) %>% round(3)
)

# check Pareto-k visually as well
par( mfrow = c(1,2) )
for ( i in names(m) ) plot( loo(m[[i]]), main = i )
par( mfrow = c(1,1) )


# ---- fig 2 linear vs non-linear fit ----

# prepare data to be predicted
d_seq <- data.frame( time_y = seq(-2,12,length.out = 50), id = NA ) %>% mutate( time = time_y + scl$Md$time )

# add predictions of expectation (epreds, based on fixed-effects only) from both linear and non-linear models
preds <- lapply( names(m) , function(i)
  d_seq %>%
    add_epred_draws( m[[i]], re_formula = NA ) %>%
    mutate( .epred = .epred * scl$SD$drs + scl$M$drs ) %>%
    median_hdi( .width = .95 ) %>%
    add_column( Model = factor( ifelse(i == "m0_linear", "Linear", "Smooths"), levels = c("Smooths","Linear"), ordered = T ) )
)

# collapse the linear and spline predictions for plotting purposes to a single file
preds <- do.call( rbind.data.frame, preds )

# plot linear and spline fits over each other
preds %>%
  ggplot( aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, color = Model, fill = Model) ) +
  geom_ribbon( alpha = .2 , color = NA ) +
  geom_line( size = 2.5 , alpha = .75 ) +
  scale_y_continuous(name = "DRS-2", limits = c(119,145), breaks = seq(120,150,10), labels = seq(120,150,10) ) +
  scale_x_continuous(name = "Time from surgery (years)", limits = c(-2,12), breaks = seq(-2,12,2), labels = seq(-2,12,2) ) +
  scale_color_manual( values = c("black",cbPal[8]) ) +
  scale_fill_manual( values = cbPal[c(1,8)] ) +
  theme( legend.position = c(0.15,0.21), legend.key.width = unit(2.6,"cm"), legend.key.height = unit(1.5,"cm") )

# save as Fig 2
ggsave( "figures/Fig 2 linear vs non-linear fit.tiff", dpi = 300, width = 9.64, height = 6.54 )
ggsave( "figures/Fig 2 linear vs non-linear fit.png", dpi = 600, width = 9.64, height = 6.54 )

# ---- stats for in-text reporting ----

# compare the models via PSIS-LOO
loo_compare( m$m0_linear, m$m0_spline )[ 2, paste0( c("elpd","se"), "_diff" ) ] %>% # extract ELPD_dif and its SE
  # re-format the object for data wrangling
  t() %>% as.data.frame() %>%
  # flip the sign of elpd_diff such that positive means m0_linear had better predictive performance
  mutate( elpd_diff = ifelse( rownames( loo_compare(m$m0_linear,m$m0_spline) )[1] == "m$m0_linear" , -elpd_diff, elpd_diff ) ) %>%
  # add 95% CI
  mutate( ELPD_dif = sprintf( "%.2f", round( elpd_diff, 2) ),
          SE_dif = sprintf( "%.2f", round( se_diff, 2 ) ),
          `95% CI` = paste0( "[", sprintf( "%.2f", qnorm( .025, elpd_diff, se_diff ) %>% round(2) ),
                             ", ", sprintf( "%.2f", qnorm( .975, elpd_diff, se_diff ) %>% round(2) ), "]"
                             ) ) %>%
  # print the results (note, the ne)
  select( ELPD_dif, SE_dif, `95% CI` ) %>% t()

# expected cognitive decline in DRS-2 points/year as approximated by the linear model (i.e., the estimand #1)
spread_draws( m$m0_linear, `b_.*` , regex = T ) %>% # extract parameter estimates
  # calculate model group-level intercept and slope summaries in the correct scale
  median_hdi( b_Intecept = (b_Intercept * scl$SD$drs ) + scl$M$drs,
              b_slope = b_time * scl$SD$drs,
              .width = .95 ) %>%
  # tidy-up the table
  select( starts_with("b_") ) %>% mutate_all( function(x) x = sprintf( "%.2f", round(x, 2) ) ) %>%
  matrix( nrow = 2, ncol = 3, byrow = T, dimnames = list( c("intercept","slope"), c("b","PPI_low","PPI_upp") ) ) %>%
  as.data.frame() %>% mutate( `b [95% PPI]` = paste0(b," [",PPI_low,", ",PPI_upp,"]") ) %>% select( `b [95% PPI]`)


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/desc_trajectories.txt" )
