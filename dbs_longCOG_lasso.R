# All analyses reported in the article (Mana et al., in review) ran in R version 4.2.0 (2022-04-22),
# on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.6.

# I used the following versions of packages employed: dplyr_1.0.9, missMDA_1.18, psych_2.2.5, brms_2.17.0.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "missMDA", # for imputation
  "psych", # for EFA
  "brms" # for Bayesian model fitting / interface with Stan
)

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set some values for later
imp = 100 # number of multiple imputations to account for missing pre-surgery data
s = 87542 # seed for reproducibility

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# read the data set and prepare subsets for individual analyses
d0 <- read.csv( "data/20220508_dbs_longCOG_data.csv" , sep = "," )
d1 <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data
d2 <- d1[ d1$ass_type == "pre" , ] # only pre-surgery assessments of included patients


# ---- pre-surgery cognitive profile  ----

# for EFA keep only id and cognitive tests in d2
d2 <- d2[ , c( 2, which(names(d2) == "tmt_a"):which(names(d2) == "fp_dr"), which(names(d2) %in% paste0("staix",1:2)) ) ]

# change names such that they get correct label in post-processing
colnames(d2)[-1] <- paste0( "b_", colnames(d2)[-1] )
tests <- colnames(d2)[-1]

# log-transform reaction times before the analysis
for ( i in c( paste0("b_tmt_", c("a","b")), paste0("b_pst_", c("d","w","c")) ) ) d2[[i]] <- log( d2[[i]] )

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
# loop through all 100 data sets and  three to eight latent factors for each imputation
efa <- lapply( 1:imp, function(i)
  lapply( 3:8, function(j) fa( d2.imp$res.MI[[i]], # one-by-one use each imputed data 
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
doms <- c("proc_spd", # loaded on primarily by PST, the first factor in 82% data sets
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


# ---- pre-processing for longitudinal analyses  ----

# before merging compute scaling values for DRS-2, BDI-II, LEDD, age and time
scl <- list( M = list( drs = mean( d1$drs_tot , na.rm = T ), # 136.90
                       bdi = mean( d1$bdi , na.rm = T ), # 10.95
                       led = mean( d1$ledd_mg , na.rm = T ), # 1197.21
                       age = mean( d1$age_ass_y, na.rm = T ) # 59.63
                       ),
             SD = list( drs = sd( d1$drs_tot , na.rm = T ), # 7.72
                        bdi = sd( d1$bdi , na.rm = T ), # 7.19
                        led = sd( d1$ledd_mg , na.rm = T ), # 687.20
                        age = sd( d1$age_ass_y, na.rm = T ) # 8.27
                        ),
             # add median time of pre-surgery assessment for GLMMs intercepts s
             Md = list(
               time = -median( d1[d1$ass_type == "pre", ]$time_y , na.rm = T ) # 0.30
               )
             )

# merge longitudinal d1 with baseline factor scores (joining by id)
d3 <- lapply( 1:imp, function(i) d1 %>%
                # first prepare a pre-surgery df with id and factor scores for each patient
                left_join( cbind.data.frame( id = d2$id, efa[[i]][[nf-2]]$scores ) , by = "id" ) %>%
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
                        proc_spd, epis_mem, verb_wm, visp_mem, set_shift, anxiety, visp_wm # pre-surgery cognition
                ) %>%
                # add raw test scores for model comparisons
                left_join( cbind.data.frame( id = d2$id, d2.imp$res.MI[[i]] ), by = "id" )
              )

# loop across all imputations to get means and SDs of the pre-surgery cognitive domains and tests,
# then transform the pre-surgery cognition in each data set to (pre-surgery) zero mean, unit SD variables
for ( i in 1:imp ) {
  # start with pre-processing the cognitive domains
  for ( j in doms ) {
    # calculate scaling values
    scl$M[[j]][[i]] <- efa[[i]][[nf-2]]$scores[,j] %>% mean()
    scl$SD[[j]][[i]] <- efa[[i]][[nf-2]]$scores[,j] %>% sd()
    # scale in the jth imputed data set
    d3[[i]][[j]] <- case_when(
      # all but anxiety measures will be inverse such that parameters
      # can be interpreted as effect of deficit in said measure
      j == "anxiety" ~ ( d3[[i]][[j]] - scl$M[[j]][[i]] ) / scl$SD[[j]][[i]],
      j != "anxiety" ~ ( scl$M[[j]][[i]] - d3[[i]][[j]] ) / scl$SD[[j]][[i]]
    )
  }
  # next pre-process single cognitive tests
  for ( j in tests ) {
    # calculate scaling values
    scl$M[[j]][[i]] <- d2.imp$res.MI[[i]][,j] %>% mean()
    scl$SD[[j]][[i]] <- d2.imp$res.MI[[i]][,j] %>% sd()
    # scale in the jth imputed data set
    if ( j %in% c( paste0("b_tmt_",c("a","b")), paste0("b_pst_",c("d","w","c")), paste0("b_staix",1:2) ) ) {
      # all but reaction speed and anxiety measures will be inversed such that parameters
      # can be interpreted as effect of deficit in said measure
      d3[[i]][[j]] <- ( d3[[i]][[j]] - scl$M[[j]][[i]] ) / scl$SD[[j]][[i]]
    } else d3[[i]][[j]] <- ( scl$M[[j]][[i]] - d3[[i]][[j]] ) / scl$SD[[j]][[i]]
  }
}

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .99 # adapt_delta parameter

# save the data and scaling values for post-processing
saveRDS( list(d1 = d1, d3 = d3, scl = scl), "data/long_dat.rds" )


# ---- models set-up ---- 

# set-up the linear models
f <- list(
  m1_lasso_doms = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", doms, sep = " * ", collapse = " + "), " + (1 + time || id)" ) %>% as.formula() %>% bf(),
  m2_lasso_tests = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", tests, sep = " * ", collapse = " + "), " + (1 + time || id)" ) %>% as.formula() %>% bf(),
  m3_flat_doms = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", doms, sep = " * ", collapse = " + "), " + (1 + time || id)" ) %>% as.formula() %>% bf()
)

# set-up priors for the model m1 (cognitive domains lasso)
p <- list(
  m1_lasso_doms = c(
    # fixed effects
    prior( normal(0.3, .1), class = Intercept ),
    prior( lasso(1), class = b ),
    # random effects
    prior( normal(0, .1), class = sd, coef = Intercept, group = id ),
    prior( normal(0, .1), class = sd, coef = time, group = id ),
    # other distributional parameters
    prior( exponential(1), class = sigma ),
    prior( gamma(2, 0.1), class = nu )
    )
  )

# add priors for models m2 (cognitive tests lasso) and m3 (congnitive domains flat priors)
p$m2_lasso_tests <- p$m1_doms # the cognitive tests lasso model have identical priors to m1
p$m3_flat_doms <- NULL # default priors for the cognitive domains flat model


# ---- primary models fitting ----

# prepare a list for single models
m <- list()

# conduct model fitting
for ( i in names(f) ) m[[i]] <- brm_multiple( formula = f[[i]], family = student(), prior = p[[i]],
                                              data = d3, sample_prior = T, seed = s, chains = ch,
                                              iter = it, warmup = wu, control = list( adapt_delta = ad ),
                                              file = paste0( "models/",i,".rds" ),
                                              save_model = paste0("models/",i,".stan")
                                              )


# ----------- soft model checking -----------

# print the highest Rhat to get an idea whether chains converged
sapply( names(m) , function(i) max(m[[i]]$rhats ) )


# ---- PSIS-LOO for model model comparisons ----

# clean the environment
rm( list = ls()[ !( ls() %in% c( "d3", "imp", "m" ) ) ] )
m$m3_flat_doms <- NULL
gc()

# prepare a list for PSIS-LOO estimates
l <- list()

# fill the l list with PSIS-LOO of m1 and m2 (which are to be compared)
for ( i in names(m)[1:2] ) {
  for ( j in 1:imp ) {
    l[[i]][[j]] <- loo( m[[i]], newdata = d3[[j]])
    saveRDS( l, "models/lasso_psis_loo.rds" ) # save after each iteration to keep it safe
  }
}
