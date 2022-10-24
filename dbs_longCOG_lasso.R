# Ran in R version 4.2.0 (2022-04-22), on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.6.

# I used the following versions of packages employed: dplyr_1.0.9, tidyverse_1.3.1, brms_2.17.0,
# tidybayes_3.0.2, ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", "tidyverse", # for data wrangling
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
ad = .95 # adapt_delta parameter

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_classic(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# read the data set
d <- read.csv( "data/20220508_dbs_longCOG_data.csv", sep = "," ) %>% filter( included == 1 ) # raw data
d.imp <- lapply( 1:imp, function(i) read.csv( paste0("data/imputed/imputed_df_",i,".csv"), sep = "," ) ) # imputed pre-surgery data sets
efa <- readRDS( "models/factanal.rds" ) # EFA models including regression-based factor scores

# read a file containing mapping of variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")


# ----------- pre-processing  -----------

# before merging compute scaling values for DRS-2, BDI-II, LEDD, age and time
scl <- list( M = list( drs = mean( d$drs_tot , na.rm = T ), # 136.90
                       bdi = mean( d$bdi , na.rm = T ), # 10.95
                       led = mean( d$ledd_mg , na.rm = T ), # 1197.21
                       age = mean( d$age_ass_y, na.rm = T ) # 59.63
                       ),
             SD = list( drs = sd( d$drs_tot , na.rm = T ), # 7.72
                        bdi = sd( d$bdi , na.rm = T ), # 7.19
                        led = sd( d$ledd_mg , na.rm = T ), # 687.20
                        age = sd( d$age_ass_y, na.rm = T ) # 8.27
                        ),
             # add median time of pre-surgery assessment for GLMMs intercepts s
             Md = list(
               time = -median( d[d$ass_type == "pre", ]$time_y , na.rm = T ) # 0.30
               )
             )

# merge longitudinal d with baseline factor scores in efa and test scores in d.imp (joining by id)
df <- lapply( 1:imp, function(i) d %>%
                # first prepare a pre-surgery df with id and factor scores for each patient
                left_join( cbind.data.frame( id = d.imp[[i]]$id, efa[[i]][[nf-2]]$scores ) , by = "id" ) %>%
                filter( complete.cases(drs_tot) ) %>% # get rid of three dummy rows due to more than one stimulation parameter (no DRS-2)
                # scale all DRS-2, BDI-II, LEDD and time already
                mutate( time = time_y + scl$Md$time,
                        drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
                        bdi = ( bdi - scl$M$bdi ) / scl$SD$bdi,
                        led = ( ledd_mg - scl$M$led ) / scl$SD$led,
                        age = ( age_ass_y - scl$M$age ) / scl$SD$age,
                        sex = as.factor( sex ), # for better estimation of BDI in the second (covariate) model
                        cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ) # right censoring for DRS == 144
                        ) %>%
                # keep only variables of interest
                select( id, time, drs, cens_drs, bdi, led, age, sex, # outcomes, demographics, clinics
                        exec_fun, epis_mem, verb_wm, visp_mem, set_shift, anxiety, visp_wm # pre-surgery cognition
                        ) %>%
                # add raw test scores for model comparisons
                left_join( cbind.data.frame( d.imp[[i]] ), by = "id" )
              )

# extract domain and test names for further analyses
doms <- efa[[1]][[nf-2]]$loadings %>% colnames()
tests <- d.imp[[1]] %>% select(-id) %>% colnames()

# loop through all imputations to get means and SDs of the pre-surgery cognitive domains and tests,
# then transform the pre-surgery cognition in each data set to (pre-surgery) zero mean, unit SD variables
for ( i in 1:imp ) {
  
  # start with pre-processing the cognitive domains
  for ( j in doms ) {
    
    # calculate scaling values
    scl$M[[j]][[i]] <- efa[[i]][[nf-2]]$scores[,j] %>% mean()
    scl$SD[[j]][[i]] <- efa[[i]][[nf-2]]$scores[,j] %>% sd()
    
    # scale in the jth imputed data set
    df[[i]][[j]] <- case_when(
      # all but anxiety measures will be inverse such that parameters
      # can be interpreted as effect of deficit in said measure
      j == "anxiety" ~ ( df[[i]][[j]] - scl$M[[j]][[i]] ) / scl$SD[[j]][[i]],
      j != "anxiety" ~ ( scl$M[[j]][[i]] - df[[i]][[j]] ) / scl$SD[[j]][[i]]
    )
  }
  
  # next pre-process single cognitive tests
  for ( j in tests ) {
    
    # calculate scaling values
    scl$M[[j]][[i]] <- d.imp[[i]][,j] %>% mean()
    scl$SD[[j]][[i]] <- d.imp[[i]][,j] %>% sd()
    
    # scale in the jth imputed data set
    if ( j %in% c( paste0("sc_tmt_",c("a","b")), paste0("sc_pst_",c("d","w","c")), paste0("sc_staix",1:2) ) ) {
      # all but reaction speed and anxiety measures will be inversed such that parameters
      # can be interpreted as effect of deficit in said measure
      df[[i]][[j]] <- ( df[[i]][[j]] - scl$M[[j]][[i]] ) / scl$SD[[j]][[i]]
    } else df[[i]][[j]] <- ( scl$M[[j]][[i]] - df[[i]][[j]] ) / scl$SD[[j]][[i]]
  }
  
}

# save the data and scaling values
saveRDS( list(d = d, df = df, scl = scl, tests = tests, doms = doms), "data/longitudinal_df.rds" )


# ---- models set-up ---- 

# set-up the linear models
f <- list(
  m1_lasso_doms = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", doms, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf(),
  m2_lasso_tests = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", tests, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf()
)

# set-up priors
p <- list()

# fill-in with identical prior templates for all three models
for ( i in names(f) ) p[[i]] <- c(
  # fixed effects
  prior( normal(0.3, .1), class = Intercept ),
  prior( lasso(1), class = b ),
  # random effects
  prior( normal(0, .1), class = sd, coef = Intercept, group = id ),
  prior( normal(0, .1), class = sd, coef = time, group = id ),
  prior( lkj(2), class = cor ),
  # other distributional parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu )
)


# ---- models fitting ----

# prepare a list for single models
m <- list()

# conduct model fitting
for ( i in names(f) ) m[[i]] <- brm_multiple( formula = f[[i]], family = student(), prior = p[[i]],
                                              data = df, sample_prior = T, seed = s, chains = ch,
                                              iter = it, warmup = wu, control = list( adapt_delta = ad ),
                                              file = paste0( "models/",i,".rds" ),
                                              save_model = paste0("models/",i,".stan")
                                              )


# ----------- soft model checking -----------

# print the highest Rhat to get an idea whether chains converged
sapply( names(m) , function(i) max(m[[i]]$rhats, na.rm = T ) )


# ---- model comparisons ----

# clean the environment
rm( list = ls()[ !( ls() %in% c( "df", "imp", "m" ) ) ] )
gc()

# compute PSIS-LOO for each imputation in the primary models
# first read loo if there is already something computed
if ( file.exists("models/lasso_psis_loo.rds") ) l <- readRDS( "models/lasso_psis_loo.rds" )

# if no PSIS-LOO was already computed, start from a scratch
if (!exists("l") ) {
  
  l <- list()
  for ( i in names(m) ) {
    for ( j in 1:imp ) {
      l[[i]][[j]] <- loo( m[[i]] , newdata = df[[j]] )
      saveRDS( l, "models/lasso_psis_loo.rds" ) # save after each iteration
      print( paste0(i,", dataset #",j) ) # print the number of the last data set with PSIS-LOO
    }
  }

  # otherwise continue from the next data set after the last one with already computed PSIS-LOO
} else {
  for ( i in names(m) ) {
    for ( j in 1:(imp-1) ) {
      j = length(l[[i]])
      if ( j < imp ) {
        l[[i]][[j+1]] <- loo( m[[i]] , newdata = df[[j+1]] )
        saveRDS( l, "models/lasso_psis_loo.rds" )
        print(paste0(i,", dataset #",j+1) ) # print the number of the last data set with PSIS-LOO
      }
    }
  }
}
