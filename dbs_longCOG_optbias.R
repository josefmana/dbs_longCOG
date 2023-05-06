# A neat small script that checks for overoptimistic bias in the best case scenario for choosing significant predictors of an
# outcome from a number of potential predictors via two-stage procedure that first pre-selects potential predictors via simple correlations
# and then estimates effect sizes via multiple regressions compared to Bayesian Lasso.

# list packages to be used
pkgs <- c( "rstudioapi", # setting working directory via RStudio API
           "conflicted", # conflict resolutions
           "tidyverse", "dplyr", "reshape2", # data wrangling
           "MASS", # multivariate normal distribution
           "brms", "bayestestR", # Bayesian model fitting and summaries
           "ggplot2", "patchwork" # plotting
           )

# load or install each of the packages as needed
for (i in pkgs) {
  if (i %in% rownames (installed.packages()) == F) install.packages(i) # install if it ain't installed yet
  if ( i %in% names (sessionInfo() $otherPkgs) == F ) library(i, character.only = T ) # load if it ain't loaded yet
}

# check for conflicts among packages and set preferences
conflict_scout()
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# set working directory (works in RStudio only)
setwd( dirname(getSourceEditorContext()$path) )

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_classic(base_size = 18) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# formats in which figures are to be saved
forms <- c(".jpg",".png",".tiff")


# ---- prepare data-generating function ----

# prepare a function on varying number of subjects, number of predictors, varying pre-selection threshold and effect sizes
# that prints a set of p-values (for two-step method) and their Bayesian equivalents (for the Bayesian Lasso), one for each
# potential predictor; takes the following inputs:
# N = number of subjects
# n_prds = number of predictors (integer)
# eff = list of effect sizes (vector)
# cor = covariance matrix (identical with correlation in case of standard normals)
# thres = p-value used as a threshold for pre-selection in the two-step method (numeric)
sim <- function( N = 126, n_prds = 7, eff = rep(0,7), cor = diag( rep(1,7) ), thres = .1 ) {
  
  # simulate "n_prds" predictors from multivariate standard normal distribution
  t <- mvrnorm( N, rep(0,n_prds), cor )
  
  # generate the outcome
  if ( length(eff) == 1 ) o <- rnorm( N, sum(t[ , 1] * eff), 1 ) # exactly one non-zero effect
  else o <- rnorm( N, rowSums(t[ , 1:length(eff) ] * eff), 1 ) # none or more than one non-zero effects
  
  # pre-select predictors
  sel <- sapply( 1:n_prds, function(i) ifelse( cor.test( t[ , i], o )$p.value < thres, 1, 0 ) )
  
  # compute final multiple linear regression for the two-step procedure
  # use only pre-selected variables
  if( sum(sel) == 0 ) reg <- lm( o ~ 1 ) # intercept only if no variable was selected
  else reg <- lm( o ~ t[ , which(sel == 1) ] )
  
  # extract p-values
  if ( sum(sel) == 0 ) p <- NULL # there ain't no p-value if no variable was pre-selected
  else p <- as.vector( summary(reg)$coefficients[ 2:(sum(sel)+1), 4 ] )
  
  # fit the Bayesian Lasso
  lasso <- brm( as.formula( paste0( "o ~ ", paste( paste0("t",1:n_prds), collapse = " + " ) ) ), # linear model
                family = gaussian(), prior = prior( lasso(1), class = b ), # likelihood and priors
                data =  cbind.data.frame(o, t) %>% `colnames<-`( c("o", paste0("t",1:n_prds) ) ) # observed variables
                )
  
  # extract probability of direction and recalculate to a p-value equivalent
  pd <- pd_to_p( p_direction(lasso)$pd, direction = "two-sided" )[ 2:(n_prds+1) ]
  
  # sum-up the results and print them
  res <- data.frame( var = paste0("V", 1:n_prds), two_step = rep(1, n_prds), lasso = pd ) # prepare ones (for those that were not tested, technically p = 1, easier to work with than NAs)
  res$two_step[ which(sel == 1) ] <- p
  return(res) # return the results
  
}


# ---- set-up covariance matrixes ----

# first prepare a function for no-covariance cases
nocov_fun <- function(K = 7) return( diag( rep(1,K) ) ) # a function for no correlation with varying number of measures
nocov <- nocov_fun(K = 23) # no covariance for K = 23 test measures (as many as we use in the current study)

# putative covariation (correlation in this case) structure of predictors in our study
# based on grouping of tests to common clusters with ~50 % shared variability (i.e., r = .7, R2 = .49) 
yocov <- matrix( c(
  1, .7, rep(0,21), .7, 1, rep(0,21), # State-Trait Anxiety Inventory
  rep(0,2), 1, .7, rep(0,19), rep(0,2), .7, 1, rep(0,19), # Trail Making Test
  rep(0,4), 1, rep(.7,2), rep(0,16), rep(0,4), .7, 1, .7, rep(0,16), rep(0,4), rep(.7,2), 1, rep(0,16), # Digit Span
  rep(0,7), 1, .7, rep(0,14), rep(0,7), .7, 1, rep(0,14), # Spatial Span
  rep(0,9), 1, rep(0,13), # Tower of London
  rep(0,10), 1, rep(.7,2), rep(0,10), rep(0,10), .7, 1, .7, rep(0,10), rep(0,10), rep(.7,2), 1, rep(0,10), # Stroop test
  rep(0,13), 1, .7, rep(0,8), rep(0,13), .7, 1, rep(0,8), # verbal fluency
  rep(0,15), 1, rep(0,7), # Similarities (abstraction)
  rep(0,16), 1, rep(.7,4), rep(0,2), rep(0,16), .7, 1, rep(.7,3), rep(0,2), rep(0,16), rep(.7,2), 1, rep(.7,2), rep(0,2), rep(0,16), rep(.7,3), 1, .7, rep(0,2), rep(0,16), rep(.7,4), 1, rep(0,2), # RAVLT 
  rep(0,21), 1, .7, rep(0,21), .7, 1 # Family Pictures test
), ncol = 23, byrow = T )

# plot the correlation matrix
melt(yocov) %>%
  rename( "Pearson's r" = "value" ) %>%
  ggplot( aes(x = Var1, y = rev(Var2), fill = `Pearson's r` ) ) +
  geom_tile( color = "grey" ) +
  scale_fill_gradient2( low = "white", high = cbPal[7], limit = c(0,1), midpoint = .33 ) +
  labs( x = "Row variable", y = "Column variable" ) +
  scale_x_continuous( breaks = 1:23, labels = 1:23 ) +
  scale_y_continuous( breaks = 1:23, labels = 23:1 ) +
  guides( fill = guide_colourbar( barwidth = 1, barheight = 18 ) ) +
  coord_fixed()

# save it
for ( i in forms ) ggsave( paste0("figures/Fig S1 simulation correlation matrix",i ), dpi = 300, width = 13.7, height = 7.79 )


# ---- simulating the null ----

# prepare a table that shows how false positive error rates depend on the threshold for case with N = 126 (the sample size of our
# study) and K = 23 tests
# first read d0 if there is already something computed
if ( file.exists("models/sims.rds") ) for ( i in names(readRDS("models/sims.RDS")) ) assign( i, readRDS("models/sims.RDS")[[i]] )

# if it ain't there, start from a scratch
if ( !exists("d0") ) {
  
  d0 <- list() # prepare a list
  
  for ( i in c("nocov","yocov") ) {
    
    # simulate the data
    d0[[i]] <- lapply( 1:1e2,
                       function(j)
                         # add. a column denoting the number of the simulation (nominal variable) and simulate
                         cbind.data.frame( data.frame( sim = rep(j,23) ), sim( N = 126, n_prds = 23, eff = rep(0,23), cor = get(i), thres = .2 ) )
                       
                       # after simulating, tidy up the table
                       ) %>% do.call( rbind.data.frame, . ) %>% add_column( cov = i, .before = 1 )
    
    # save after each iteration
    saveRDS( list( d0 = d0 ), "models/sims.rds" )
    
  }
}


# ---- simulating effects ----

# as the code is fitting few hundreds of Stan models, it crashes every so often, so we will put some insurances in place

# fit one to five medium-size predictors (r = .3) in two settings: a lot (K = 23) and medium amount (K = 7) tests
# if it there are no simulations of this type (see l. 135), start from a scratch
if ( !exists("d1") ) {
  
  d1 <- list() # prepare a list
  
  # loop through the possibilities
  for ( i in c("lot","med") ) {
    
    d1[[i]] <- list()
    K <- case_when( i == "lot" ~ 23, i == "med" ~ 7 ) # number of potential predictors
    
    # loop through 1-5 to-be true/false positives
    for ( j in 1:5 ) {
      
      # simulate the data
      d1[[i]][[j]] <- lapply( 1:1e2,
                              function(l)
                                # add. a column denoting the number of the simulation (nominal variable) and simulate
                                cbind.data.frame( data.frame( sim = rep(l,K) ), sim( N = 126, n_prds = K, eff = rep(0,K), cor = nocov_fun(K=K), thres = .2 ) )

                         # after simulating, tidy up the table
                         ) %>%
        do.call( rbind.data.frame, . ) %>%
        add_column( true_positives = j, .before = 1 ) %>%
        add_column( no_preds = K, .before = 1 )
      
      # save after each iteration
      saveRDS( list( d0 = d0, d1 = d1 ), "models/sims.rds" )
      
      }
    }
  # otherwise continue from the next data set after the last one with already computed simulations
} else {
  
  # loop through all the levels, this time no need for creating new lists as they already exist (read from "models/sims.RDS")
  for ( i in c("lot","med") ) {
    
    if( !exists( i, where = d1 ) ) d1[[i]] <- list()
    
    K <- case_when( i == "lot" ~ 23, i == "med" ~ 10 ) # number of potential predictors
    l <- length(d1[[i]]) # number of already simulated true positives
    
    for ( j in 1:5) {
      
      # compute only settings that are not computed yet
      if ( j > l ) {
        
        # simulate the data
        d1[[i]][[j]] <- lapply( 1:1e2,
                                function(k)
                                  # add. a column denoting the number of the simulation (nominal variable) and simulate
                                  cbind.data.frame( data.frame( sim = rep(k,K) ), sim( N = 126, n_prds = K, eff = rep(0,K), cor = nocov_fun(K=K), thres = .2 ) )
                                
                                # after simulating, tidy up the table
                                ) %>%
          do.call( rbind.data.frame, . ) %>%
          add_column( true_positives = j, .before = 1 ) %>%
          add_column( no_preds = K, .before = 1 )
        
        # save after each iteration
        saveRDS( list( d0 = d0, d1 = d1 ), "models/sims.rds" )
          
      }
    }
  }
}

# finish by simulating two kinds of effects with covaried predictors where we assume total medium effect of variables that
# represent TMT (3:4), and Stroop test (11:13) inspired by Kim et al. (2014)
# (i) each variable has a proportional weight with all clusters linking to 0.3
# (ii) exactly one variable from each cluster has the total 0.3 weight
if ( !exists("d2") ) {
  
  d2 <- list() # prepare a list
  
  # calculate the simulations
  d2$cluster <- lapply( 1:1e2, function(i) cbind.data.frame( data.frame( sim = rep(i,23) ), sim( N = 126, n_prds = 23, eff = c( rep(0,2), rep(.15,2),rep(0,6),rep(.1,3),rep(0,10) ), cor = yocov, thres = .2 ) ) ) %>% do.call( rbind.data.frame, . )
  d2$test <- lapply( 1:1e2, function(i) cbind.data.frame( data.frame( sim = rep(i,23) ), sim( N = 126, n_prds = 23, eff = c( rep(0,2), .3, rep(0,7),.3,rep(0,12) ), cor = yocov, thres = .2 ) ) ) %>% do.call( rbind.data.frame, . )
    
}

# save it all
saveRDS( list( d0 = d0, d1 = d1, d2 = d2 ), "models/sims.RDS" )


# ---- prepare summaries ----

# summary of false positives conditional on null hypothesis
t0 <- do.call( rbind.data.frame, d0 ) %>%
  # re-format such that we have a column for all variables to be pushed to ggplot
  pivot_longer( cols = c("two_step","lasso"), names_to = "meth", values_to = "p" ) %>% # prepare a single column of p-values
  mutate( sig = ifelse(p < .05, 1, 0 ) ) %>% # flag findings significant on 5% level
  pivot_wider( id_cols = c("cov","sim","meth"), values_from = sig, names_from = var ) %>% # pivot again such that each potential predictor has its own column with indicator of a "hit" per simulation/method pairs
  mutate( false_positives = rowSums( across( starts_with("V") ) ) ) %>% # calculate number of "hits"/false positives across variables for each simulation/method pair
  select( cov, meth, false_positives ) %>% # keep only variables of interest for the plot
  table() %>% as.data.frame() %>% # prepare a table of frequencies of false hits per method
  # prepare names
  mutate(
    Method = case_when( meth == "two_step" ~ "Two-step procedure", meth == "lasso" ~ "Bayesian Lasso" ),
    Covariance = case_when( cov == "nocov" ~ "Independent predictors", cov == "yocov" ~ "Covaried predictors" ) %>% factor( levels = c("Independent predictors","Covaried predictors"), ordered = T )
  )

# summary of false negatives under no covariance
t1 <- lapply( names(d1), function(i)
  
  # looping through levels of nocov simulations with true predictors present
  lapply( 1:length(d1[[i]]), function(j)
    
    # for each combination of d1 individually count false positives and false negatives
    d1[[i]][[j]] %>%
      pivot_longer( cols = c("two_step","lasso"), names_to = "meth", values_to = "p" ) %>% # prepare a single column of p-values
      mutate( sig = ifelse(p < .05, 1, 0 ) ) %>% # flag findings significant on 5% level
      pivot_wider( id_cols = c("no_preds","true_positives","sim","meth"), values_from = sig, names_from = var ) %>%
      mutate(
        false_positives = rowSums( across( paste0( "V", (j+1):case_when( i == "lot" ~ 23, i == "med" ~ 7 ) ) ) ),
        false_negatives = j - rowSums( across( paste0( "V", 1:j ) ) )
      )
    
  ) %>% do.call( rbind.data.frame, . ) %>% select( -starts_with("V") ) # keep only the outcomes so that it's easier to bind the data sets

) %>% do.call( rbind.data.frame, . ) # collapse to a single long data sets

  

# ---- plotting ----

# plot false positives rates conditional on null hypothesis
t0 %>%
  # shift frequencies by one such that there is some visual representation of zero as well
  ggplot( aes(x = false_positives, y = Freq+1, fill = Method ) ) +
  geom_bar( stat = "identity", position = position_dodge( width = .5), width = .5 ) + # bars
  geom_text( aes(label = Freq), vjust = -1.0, size = 4, position = position_dodge(width = .5 ) ) + # counts
  scale_y_continuous( limits = c(0,105), breaks = seq(1,101,25), labels = seq(0,100,25), name = "Count" ) +
  scale_x_discrete( name = "False positives per one hundred simulations" ) +
  scale_fill_manual( values = cbPal[c(2,1)] ) + # colors
  theme( legend.position = "bottom" ) +
  facet_wrap( ~Covariance ) # prepare subplots

# save it
for ( i in forms ) ggsave( paste0("figures/Fig S2 simulation false positives under null",i ), dpi = 300, width = 12.8, height = 7.79 )


