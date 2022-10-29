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
imp = 100 # number of multiple imputations to account for missing pre-surgery data
s = 87542 # seed for reproducibility
nf = 7 # number of factors extracted via factor analysis

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
theme_set( theme_classic(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# read the files needed (raw data, imputed data sets, scaling values, test and domain names)
for( i in names( readRDS( "data/longitudinal_df.rds" ) ) ) assign( i , readRDS( "data/longitudinal_df.rds" )[[i]] )

# read a file containing mapping of variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")


# ---- models set-up ---- 

# set-up the linear models for drs
f.drs <- list(
  m3_doms_cov = ( paste0( "drs | cens(cens_drs) ~ 1 + age + mi(bdi) + mi(led) + ", paste("time", doms, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf() ) + student(),
  m4_tests_cov = ( paste0( "drs | cens(cens_drs) ~ 1 + age + mi(bdi) + mi(led) + ", paste("time", tests, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf() ) + student()
)

# set-up linear models for covariates with missing values
f.bdi <- bf( bdi | mi() ~ 1 + time + sex + age + mi(led) + (1 + time | id) ) + gaussian()
f.led <- bf( led | mi() ~ t2(time) + (1 | id) ) + gaussian()

# set-up priors (using the same priors for both models)
p <- c(
  # DRS-2
  prior( normal(0.3, .1), class = Intercept, resp = drs),
  prior( lasso(1), class = b, resp = drs ),
  prior( normal(0, .1), class = sd, coef = Intercept, group = id, resp = drs ),
  prior( normal(0, .1), class = sd, coef = time, group = id, resp = drs ),
  prior( exponential(1), class = sigma, resp = drs ),
  prior( gamma(2, 0.1), class = nu, resp = drs ),
  # BDI-II
  prior( normal(.6, .5), class = Intercept, resp = bdi ),
  prior( normal(0, .5), class = b, coef = time, resp = bdi ),
  prior( normal(0, .5), class = b, coef = sexmale, resp = bdi ),
  prior( normal(0, .5), class = b, coef = miled, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  # LEDD
  prior( normal(0, 100), class = Intercept, resp = led ),
  prior( normal(0, 100), class = b, coef = t2time_1, resp = led ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = led ),
  prior( exponential(1), class = sigma, resp = led ),
  # covariance structure
  prior( lkj(2), class = cor)
)


# ---- models fitting ----

# prepare a list for single models
m <- list()

# conduct model fitting
for ( i in names(f.drs) ) m[[i]] <- brm_multiple( formula = f.drs[[i]] + f.bdi + f.led, prior = p,
                                                  data = df, sample_prior = T, seed = s, chains = ch,
                                                  iter = it, warmup = wu, control = list( adapt_delta = ad ),
                                                  file = paste0( "models/",i,".rds" ),
                                                  save_model = paste0("models/",i,".stan")
                                                  )


# ----------- soft model checking -----------

# print the highest Rhat to get an idea whether chains converged
sapply( names(m), function(i)
        sapply( paste0("_",c("drs","bdi","led"),"_" ),  function(j) max( m[[i]]$rhats %>% select(contains(j)), na.rm = T ) )
        )
