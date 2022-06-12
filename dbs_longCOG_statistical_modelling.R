# All analyses reported in the article (Mana et al., in review) ran in R version 4.2.0 (2022-04-22),
# on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.4.

# I used the following versions of packages employed: dplyr_1.0.9, missMDA_1.18, psych_2.2.5, brms_2.17.0

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "missMDA", # for imputation
  "psych", # for EFA
  "brms" # for Bayesian model fitting / interface with Stan
)

# load required packages
# prints NULL if a package is already installed
sapply(
  pkgs, # packages to be loaded/installed
  function(i)
    if ( !require( i , character.only = T ) ){
      # it's important to have the 'character.only = T' command here
      install.packages(i)
      library(i)
    }
)

# set number of multiple imputations to account for missing pre-surgery data
imp = 100
s = 87542 # seed for reproducibility

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# read the data set and prepare subsets for individual analyses
d0 <- read.csv( "data/20220508_dbs_longCOG_data.csv" , sep = "," )
d1 <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data
d2 <- d1[ d1$ass_type == "pre" , ] # only pre-surgery assessments of included patients


# ----------- pre-surgery cognitive profile  -----------

# for EFA keep only id and cognitive tests in d2
d2 <- d2[ , c( 2, which(names(d2) == "tmt_a"):which(names(d2) == "fp_dr"), which(names(d2) %in% paste0("staix",1:2)) ) ]

# change staix1 and staix2 names such that they get correct label in post-processing
d2 <- d2 %>% rename( "stai_x1" = "staix1" , "stai_x2" = "staix2" ) 

# log-transform reaction times before the analysis
for ( i in c( paste0("tmt_", c("a","b")), paste0("pst_", c("d","w","c")) ) ) d2[[i]] <- log( d2[[i]] )

# normalize (center and scale) all test scores before analysis
for ( i in colnames(d2)[-1] ) d2[,i] <- as.vector( scale(d2[,i], center = T, scale = T) )

# find out the optimal number of components for multiple imputation
nb <- estim_ncpPCA( d2[,-1] , ncp.min = 0, ncp.max = 10 , nbsim = imp )

# impute via PCA-based multiple imputation (n = 100 imputations)
set.seed(s) # set seed for reproducibility
d2.imp <- MIPCA( d2[ ,- 1] , ncp = nb$ncp , nboot = imp )

# calculate parallel test for each imputed data set 
p.test <- lapply( 1:imp , function(i) fa.parallel( d2.imp$res.MI[[i]] , plot = F ) )

# look at the results of parallel tests
# note, does not replicate even with a seed set
table( sapply( 1:imp , function(i) p.test[[i]]$nfact ) )

# fit EFA to each imputed data set with 3:8 factors
efa <- lapply(
  # loop through all 100 data sets
  1:imp, function(i)
    lapply(
      # loop through three to eight latent factors for each imputation
      3:8, function(j)
        fa(
          d2.imp$res.MI[[i]], # one-by-one use each imputed data set
          nfactors = j, # fit 3-8 factor solutions
          rotate = "varimax", # rotate varimax to enforce orthogonality and for interpretation purposes
          scores = "regression" # compute regression scores for each patient
        )
    )
)

# one-by-one inspect all six-factor and seven-factor solutions' loading matrices
# create a convenience function so that I don't go crazy immediately
print_load <- function( i, c = .4, f = 7 ) print( efa[[i]][[f-2]]$loadings, cutoff = c, sort = T )

# choosing 7-factor solution due to good performance indexes,
# and theoretically sound loading patterns across imputed data sets
nf = 7

# prepare an array for labels of the seven-factor solution factors
# create an empty 2 (names/signs) x 7 (factors) x 100 (imputations) array
doms_sum <- array( data = NA, dim = c(2, nf, imp), dimnames = list( c("nms","sgn"), paste0("F", 1:nf), 1:imp ) )

# read the table with seven-factor labels
doms_sum["nms", , ] <- t( read.csv( "data/dbs_longCOG_efa_labels.csv" , sep = "," , row.names = 1, header = T) )

# fill-in signs of each factor in each imputation to know which scores should be reversed
doms_sum["sgn", , ] <- apply( doms_sum["nms", , ] , 2 , function(x) startsWith( x , "-") ) %>%
  t %>% as.data.frame %>% mutate( across( everything() , ~ ifelse( . == T , -1 , 1 ) ) ) %>% t

# get rid of the minus sign in labels table
doms_sum["nms", , ] <- doms_sum["nms", , ] %>%
  t %>% as.data.frame %>% mutate( across( everything() , ~ gsub( "-" , "" , . ) ) ) %>% t

# list all the domains
doms <- c(
  "proc_spd", # loaded on primarily by PST, the first factor in 82% data sets
  "epis_mem", # loaded on primarily by RAVLT, the second factor in 79% data sets
  "verb_wm", # loaded on primarily by DS, the third factor in 62% data sets
  "visp_mem", # loaded on primarily by FP, the fourth factor in 45% data sets
  "set_shift", # loaded on primarily by TMT and RAVLT-B, the fifth factor in 28% data sets
  "anxiety", # loaded on primarily by STAI, the sixth factor in 60% data sets
  "visp_wm" # loaded on primarily by SS, the seventh factor in 49% data sets
)

# switch signs where appropriate in EFA loadings and scores, and rename and sort columns
for ( i in 1:imp ) {
  for ( j in c("loadings","scores","Vaccounted") ) {
    # multiply by a diagonal matrix of 1 and -1
    if( j %in% c("loadings","scores") ) efa[[i]][[nf-2]][[j]] <- efa[[i]][[nf-2]][[j]] %*% diag(doms_sum["sgn", , i ] )
    # rename the columns
    colnames( efa[[i]][[5]][[j]] ) <- doms_sum["nms", , i ]
    # reorder the columns such that they are in the same order for each imputation
    efa[[i]][[nf-2]][[j]] <- efa[[i]][[nf-2]][[j]][, doms]
  }
}

# save all EFA models for post-processing
saveRDS( object = efa, file = "models/efa.rds" )


# ----------- pre-processing for longitudinal analyses  -----------

# before merging compute scaling values for DRS-2, BDI-II, LEDD, age and time
scl <- list(
  M = list(
    drs = mean( d1$drs_tot , na.rm = T ), # 136.90
    bdi = mean( d1$bdi , na.rm = T ), # 10.95
    led = mean( d1$ledd_mg , na.rm = T ), # 1197.21
    age = mean( d1$age_ass_y, na.rm = T ) # 59.63
  ),
  SD = list(
    drs = sd( d1$drs_tot , na.rm = T ), # 7.72
    bdi = sd( d1$bdi , na.rm = T ), # 7.19
    led = sd( d1$ledd_mg , na.rm = T ), # 687.20
    age = sd( d1$age_ass_y, na.rm = T ) # 8.27
  ),
  # add median time of pre-surgery assessment for GLMMs intercepts
  Md = list(
    time = -median( d1[d1$ass_type == "pre", ]$time_y , na.rm = T ) # 0.30
  )
)

# merge longitudinal d1 with baseline factor scores (joining by id)
d3 <- lapply(
  1:imp,
  function(i) d1 %>%
    # first prepare a pre-surgery df with id and factor scores for each patient
    left_join( cbind.data.frame(id = d2$id, efa[[i]][[nf-2]]$scores ) , by = "id" ) %>%
    filter( complete.cases(drs_tot) ) %>% # get rid of three dummy rows due to more than one stimulation parameter (no DRS-2)
    # scale all DRS-2, BDI-II, LEDD and time already
    mutate(
      time = time_y + scl$Md$time,
      drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
      bdi = ( bdi - scl$M$bdi ) / scl$SD$bdi,
      led = ( ledd_mg - scl$M$led ) / scl$SD$led,
      age = ( age_ass_y - scl$M$age ) / scl$SD$age,
      sex = as.factor( sex ), # for better estimation of BDI in the second (covariate) model
      cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ) # right censoring for DRS == 144
    ) %>%
    # keep only variables of interest
    select(
      id, time, drs, cens_drs, bdi, led, age, sex, # outcomes, demographics, clinics
      proc_spd, epis_mem, verb_wm, visp_mem, set_shift, anxiety, visp_wm # pre-surgery cognition
    )
)

# loop across all imputations to get means and SDs of the pre-surgery cognitive domains,
# then transform the pre-surgery cognition in each data set to (pre-surgery) zero mean, unit SD variables
for( i in doms ) {
  for ( j in 1:imp ) {
    # calculate scaling values
    scl$M[[i]][[j]] <- d3[[j]][[i]] %>% mean( na.rm = T )
    scl$SD[[i]][[j]] <- d3[[j]][[i]] %>% sd( na.rm = T )
    # scale in the jth imputed data set
    d3[[j]][[i]] <- case_when(
      # all but anxiety measures will be inverse such that parameters
      # can be interpreted as effect of deficit in said measure
      i == "anxiety" ~ ( d3[[j]][[i]] - scl$M[[i]][[j]] ) / scl$SD[[i]][[j]],
      i != "anxiety" ~ -( d3[[j]][[i]] - scl$M[[i]][[j]] ) / scl$SD[[i]][[j]]
    )
  }
}

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .95 # adapt_delta parameter


# ----------- base model with correlated random-effects  -----------

# set-up the linear model
f0.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time | id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf

# set-up priors
p0 <- c(
  # fixed effects
  prior( normal(0.3, .1), class = Intercept ),
  prior( normal(-.2, .1), class = b, coef = time ),
  prior( normal(0, .1), class = b, coef = proc_spd ),
  prior( normal(0, .1), class = b, coef = epis_mem ),
  prior( normal(0, .1), class = b, coef = verb_wm ),
  prior( normal(0, .1), class = b, coef = visp_mem ),
  prior( normal(0, .1), class = b, coef = set_shift ),
  prior( normal(0, .1), class = b, coef = anxiety ),
  prior( normal(0, .1), class = b, coef = visp_wm ),
  prior( normal(0, .1), class = b, coef = time:proc_spd ),
  prior( normal(0, .1), class = b, coef = time:epis_mem ),
  prior( normal(0, .1), class = b, coef = time:verb_wm ),
  prior( normal(0, .1), class = b, coef = time:visp_mem ),
  prior( normal(0, .1), class = b, coef = time:set_shift ),
  prior( normal(0, .1), class = b, coef = time:anxiety ),
  prior( normal(0, .1), class = b, coef = time:visp_wm ),
  # random effects
  prior( normal(0, .1), class = sd, coef = Intercept, group = id ),
  prior( normal(0, .1), class = sd, coef = time, group = id ),
  prior( lkj(2), class = cor ),
  # other distributional parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu )
)

# fit the model with Student response function
m <- list(
  m0_base = brm_multiple(
    formula = f0.drs, family = student(), prior = p0, data = d3,
    sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
    control = list( adapt_delta = ad ), file = "models/m0_base.rds"
  )
)


# ----------- primary model with uncorrelated random-effects  -----------

# set-up the linear model
f1.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time || id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf

# set-up priors
p1 <- p0[ -which(p0$class == "cor"), ]

# fit the model with Student response function
m$m1_nocov <- brm_multiple(
  formula = f1.drs, family = student(), prior = p1, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m1_nocov.rds"
)


# ----------- sensitivity check with flat priors  -----------

# the same linear model as for m1
f2.drs <- f1.drs

# will be using brms default priors
p2 <- NULL

# fit the model with Student response function
m$m2_flat_priors <- brm_multiple(
  formula = f2.drs, family = student(), prior = p2, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m2_flat_priors.rds"
)


# ----------- model with covariates  -----------

# set contrast for sex as a factor
for( i in 1:imp ) contrasts(d3[[i]]$sex) <- -contr.sum(2)/2 # female = -0.5, male = 0.5

# set-up the linear model for outcome
f3.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + age + mi(bdi) + mi(led) + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time || id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf + student()

# set-up linear models for covariates with missing values
f3.bdi <- bf( bdi | mi() ~ 1 + time + sex + age + mi(led) + (1 + time || id) ) + gaussian()
f3.led <- bf( led | mi() ~ t2(time) + (1 | id) ) + gaussian()

# set-up priors
p3 <- c(
  # DRS-2
  prior( normal(0.3, .1), class = Intercept, resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time, resp = drs ),
  prior( normal(0, .1), class = b, coef = age, resp = drs ),
  prior( normal(0, .1), class = b, coef = mibdi, resp = drs ),
  prior( normal(0, .1), class = b, coef = miled, resp = drs ),
  prior( normal(0, .1), class = b, coef = proc_spd, resp = drs ),
  prior( normal(0, .1), class = b, coef = epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = set_shift, resp = drs ),
  prior( normal(0, .1), class = b, coef = anxiety, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:proc_spd, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:set_shift, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:anxiety, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_wm, resp = drs ),
  prior( normal(0, .1), class = sd, coef = Intercept, group = id, resp = drs ),
  prior( normal(0, .1), class = sd, coef = time, group = id, resp = drs ),
  prior( exponential(1), class = sigma, resp = drs ),
  prior( gamma(2, 0.1), class = nu, resp = drs ),
  # BDI-II
  prior( normal(.6, .5), class = Intercept, resp = bdi ),
  prior( normal(0, .5), class = b, coef = time, resp = bdi ),
  prior( normal(0, .5), class = b, coef = sex1, resp = bdi ),
  prior( normal(0, .5), class = b, coef = miled, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  # LEDD
  prior( normal(0, 100), class = Intercept, resp = led ),
  prior( normal(0, 100), class = b, coef = t2time_1, resp = led ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = led ),
  prior( exponential(1), class = sigma, resp = led )
)

# fit the model as defined above
m$m3_wcov <- brm_multiple(
  formula = f3.drs + f3.bdi + f3.led, prior = p3, data = d3,
  sample_prior = F, seed = s, chains = ch, iter = it, warmup = wu, # not saving priors to spare some memory (and because I don`t need them for this model)
  control = list( adapt_delta = .99 ), file = "models/m3_wcov.rds"
)


# ----------- extract stan code -----------

# for each model extract raw stan code and save it for sharing
for ( i in names(m) ) write.table( x = stancode(m[[i]]),
                                   file = paste0(getwd(), "/models/", i, "_stancode.txt")
                                   )


# ----------- soft model checking -----------

# print the highest Rhat to get an idea whether chains converged
sapply( names(m) , function(i) max(m[[i]]$rhats ) )

# keep only the primary model for next calculations to spare time
m <- m$m1_nocov

# clean environment to spare RAM
rm( list = ls()[ !( ls() %in% c( "d3", "imp", "m" ) ) ] )
gc()

# compute PSIS-LOO for each imputation in the primary model
# first read loo if there is already something compited
if ( file.exists("models/dbs_longCOG_psis-loo.rds") ) l <- readRDS( "models/dbs_longCOG_psis-loo.rds" )

# if no PSIS-LOO was already computed, start from a scratch
if (!exists("l") ) {
  l <- list()
  for ( i in 1:imp ) {
    l[[i]] <- loo( m , newdata = d3[[i]], pointwise = T )
    saveRDS( l, "models/dbs_longCOG_psis-loo.rds" ) # save after each iteration
    print(i) # print the number of the last data set with PSIS-LOO
  }
  # otherwise continue from the next data set after the last one with already computed PSIS-LOO
} else {
  for ( i in 1:(imp-1) ) {
    i = length(l)
    if ( i < imp ) {
      l[[i+1]] <- loo( m , newdata = d3[[i+1]], pointwise = T )
      saveRDS( l, "models/dbs_longCOG_psis-loo.rds" )
      print(i+1) # print the number of the last data set with PSIS-LOO
    }
  }
}
