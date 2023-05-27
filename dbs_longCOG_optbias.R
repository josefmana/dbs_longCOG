# This is a script to check false error rates in the best case scenario for choosing significant predictors of a longitudinal outcome
# from a number of potential predictors via two-stage procedure that first pre-selects potential predictors via univariate regression
# and then continues to multiple regressions compared to Bayesian Lasso.

# list packages to be used
pkgs <- c( "rstudioapi", # setting working directory via RStudio API
           "tidyverse", "dplyr", "reshape2", # data wrangling
           "MASS", # multivariate normal distribution
           "lmerTest", "brms", "bayestestR", # model fitting and summaries
           "ggplot2", # plotting
           "english" # changing numbers to words
           )

# load or install each of the packages as needed
for (i in pkgs) {
  if (i %in% rownames (installed.packages()) == F) install.packages(i) # install if it ain't installed yet
  if ( i %in% names (sessionInfo() $otherPkgs) == F ) library(i, character.only = T ) # load if it ain't loaded yet
}

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

# read the structure of our data set
struct <- read.csv( "data/data_struct.csv", sep = "," )

# use all parallel CPU cores for computing
options( mc.cores = parallel::detectCores() )


# ---- prepare data-generating function ----

# prepare a function on varying number of subjects, number of predictors, varying pre-selection threshold and effect sizes
# that prints a set of p-values (for two-step method) and their Bayesian equivalents (for the Bayesian Lasso), one for each
# potential predictor; takes the following inputs:
# b0 = population-level pre-surgery intercept
# b1 = population-level annual cognitive decline
# n_prds = number of predictors (integer)
# eff = list of effect sizes (vector)
# cor = covariance matrix (identical with correlation in case of standard normals)
# thres = p-value used as a threshold for pre-selection in the two-step method (numeric)
# struct = longitudinal data structure
sim <- function( b0 = 0, b1 = -.2, n_prds = 7, eff = rep(0,7), cor = diag( rep(1,7) ), thres = .1, struct = struct ) {
  
  # extract number of patients and their IDs
  ID <- unique(struct$id)
  N <- length(ID)
  
  # simulate patient-specific predictors/variation
  t <- mvrnorm( N, rep(0,n_prds), cor ) %>% `rownames<-`( ID ) # "n_prds" predictors from multivariate standard normal distribution
  u <- mvrnorm( N, c(0,0), diag( c(.05^2,.05^2) ) ) %>% `row.names<-`( ID ) # patient-specific intercepts and slopes
  
  # prepare a data set combining pre-defined data structure (struct) and generative regression parameters
  d <- struct %>% mutate( b0 = b0, b1 = b1, # population-level slope
                          c0 = sapply( 1:nrow(struct), function(i) u[ struct$id[i], 1 ] ), # patient-specific intercept
                          c1 = sapply( 1:nrow(struct), function(i) u[ struct$id[i], 2 ] ) # patient-specific slope
                          )
  
  # add patients' cognitive profile
  for ( i in 1:ncol(t) ) d[ , paste0("p",i) ] <- sapply( 1:nrow(d), function(j) t[ d$id[j], i ] )
  
  # generate the outcome
  if ( length(eff) == 1 ) d$o <- with( d, rnorm( N, ( b0 + c0 )  + ( b1 + c1 ) * time + d[ , "p1"] * eff[1] * time, 1 ) ) # exactly one non-zero effect
  else d$o <- with( d, rnorm( N, ( b0 + c0 )  + ( b1 + c1 ) * time + rowSums( sapply( 1:n_prds, function(i) d[ , paste0("p",i) ] * eff[i] * d$time ) ), 1 ) ) # none or more than one non-zero effects
  
  # pre-select predictors
  sel <- sapply( paste0("p",1:n_prds), function(i) ifelse( summary( lmer( as.formula( paste0( "o ~ 1 + time * ", i, "+ (1 + time | id)") ), data = d ) )[["coefficients"]][3,5] < thres, 1, 0 ) )
  
  # compute final multiple linear regression for the two-step procedure
  # use only pre-selected variables
  if( sum(sel) == 0 ) reg <- lmer( o ~ 1 + time + (1 + time | id ), data = d ) # intercept only if no variable was selected
  else reg <- lmer( as.formula( paste0( "o ~ 1 + ", paste("time", names(sel)[sel == 1], sep = " * ", collapse = " + "), " + (1 + time | id )" ) ), data = d )
  
  # extract p-values
  if ( sum(sel) == 0 ) p <- NULL # there ain't no p-value if no variable was pre-selected
  else p <- as.vector( summary(reg)$coefficients[ paste0("time:",names(sel)[sel == 1]) , 5 ] )
  
  # fit the Bayesian Lasso
  lasso <- brm( as.formula( paste0( "o ~ 1 + ", paste( "time", paste0("p",1:n_prds), sep = " * ", collapse = " + " ), " + (1 + time | id)" ) ), # linear model
                family = gaussian(), prior = prior( lasso(1), class = b ), # likelihood and priors
                data = d # observed variables
                )
  
  # extract probability of direction and recalculate to a p-value equivalent
  pd <- pd_to_p( with( p_direction(lasso), pd[ Parameter %in% paste0("b_time:p",1:n_prds)] ), direction = "two-sided" )
  
  # sum-up the results and print them
  res <- data.frame( var = paste0("p", 1:n_prds), two_step = rep(1, n_prds), lasso = pd ) # prepare ones (for those that were not tested, technically p = 1, easier to work with than NAs)
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

# short description will be here

# as the code is fitting few hundreds of Stan models, it crashes every so often, so we will put some insurances in place
if ( file.exists("models/sims.rds") ) for ( i in names(readRDS("models/sims.RDS")) ) assign( i, readRDS("models/sims.RDS")[[i]] )

# if it ain't there, start from a scratch
if ( !exists("d0") ) {
  
  d0 <- list() # prepare a list
  
  # loop through two diffetent time slopes (b1 = 0 = null effect and b1 = -0.3 = small to moderate decline)
  for ( i in 1:2 ) {
    
    b <- case_when( i == 1 ~ 0, i == 2 ~ -.3 ) # prepare the slope
    d0[[i]] <- list() # prepare a list
    
    # loop through covariance structures for test-level inference
    for ( j in c("nocov","yocov") ) {
      
      # simulate the data
      d0[[i]][[paste0(j,23)]] <- lapply(
        # loop one hundred times
        1:1e2, function(k)
          # add a column denoting the number of the simulation (nominal variable) and simulate
          cbind.data.frame(
            data.frame( sim = rep(k,23) ), # column labeling single simulations 
            sim( b0 = 0, b1 = b, n_prds = 23, eff = rep(0,23), cor = get(j), thres = .2, struct = struct ) # synthetic data analysis results
          )
      ) %>%
      # after simulating, tidy up the table
      do.call( rbind.data.frame, . ) %>%
        add_column( decl = b, .before = 1 ) %>%
        add_column( type = paste0(j,23), .before = 1 )
      
      # save after each iteration
      saveRDS( list( d0 = d0 ), "models/sims.rds" )
      
    }
    
    # add null model of the factor scores (i.e., seven independent predictors)
    d0[[i]]$nocov7 <- lapply(
      # loop one hundred times
      1:1e2, function(j)
        # add a column denoting the number of the simulation (nominal variable) and simulate
        cbind.data.frame(
          data.frame( sim = rep(j,7) ), # column labeling single simulations
          sim( b0 = 0, b1 = b, n_prds = 7, eff = rep(0,7), cor = nocov_fun(7), thres = .2, struct = struct ) # synthetic data analysis results
        )
    ) %>%
    # after simulating, tidy up the table
    do.call( rbind.data.frame, . ) %>%
      add_column( decl = -(i-1)/10, .before = 1 ) %>%
      add_column( type = "nocov7", .before = 1 )
    
    # save once again
    saveRDS( list( d0 = d0 ), "models/sims.rds" )
    
  }

# otherwise continue from the next data set after the last one with already computed simulations
} else {
  
  # loop through several time slopes, this time no need for creating new lists as they already exist (read from "models/sims.RDS")
  for ( i in 1:2 ) {
    
    l <- length(d0[[i]]) # number of already simulated data sets
    
    # compute only settings that are not computed yet
    if ( i > l ) {
      
      b <- case_when( i == 1 ~ 0, i == 2 ~ -.3 ) # prepare the slope
      d0[[i]] <- list() # prepare a list
      
      # loop through covariance structures for test-level inference
      for ( j in c("nocov","yocov") ) {
        
        # simulate the data
        d0[[i]][[paste0(j,23)]] <- lapply(
          # loop one hundred times
          1:1e2, function(k)
            # add. a column denoting the number of the simulation (nominal variable) and simulate
            cbind.data.frame(
              data.frame( sim = rep(k,23) ), # column labeling single simulations 
              sim( b0 = 0, b1 = b, n_prds = 23, eff = rep(0,23), cor = get(j), thres = .2, struct = struct ) # synthetic data analysis results
            )
        ) %>%
        # after simulating, tidy up the table
        do.call( rbind.data.frame, . ) %>%
          add_column( decl = b, .before = 1 ) %>%
          add_column( type = paste0(j,23), .before = 1 )
        
        # save after each iteration
        saveRDS( list( d0 = d0 ), "models/sims.rds" )
        
      }
      
      # add null model of the factor scores (i.e., seven independent predictors)
      d0[[i]]$nocov7 <- lapply(
        # loop one hundred times
        1:1e2, function(j)
          # add a column denoting the number of the simulation (nominal variable) and simulate
          cbind.data.frame(
            data.frame( sim = rep(j,7) ), # column labeling single simulations
            sim( b0 = 0, b1 = b, n_prds = 7, eff = rep(0,7), cor = nocov_fun(7), thres = .2, struct = struct ) # synthetic data analysis results
          )
      ) %>%
      # after simulating, tidy up the table
      do.call( rbind.data.frame, . ) %>%
        add_column( decl = b, .before = 1 ) %>%
        add_column( type = "nocov7", .before = 1 )
      
      # save once again
      saveRDS( list( d0 = d0 ), "models/sims.rds" )
      
    }
  }
}




# ---- prepare summaries ----

# summary of false positives conditional on null hypothesis
t0 <- lapply(1:2, function(i)
  do.call( rbind.data.frame, d0[[i]] ) %>%
    # re-format such that we have a column for all variables to be pushed to ggplot
    pivot_longer( cols = c("two_step","lasso"), names_to = "meth", values_to = "p" ) %>% # prepare a single column of p-values
    mutate( sig = ifelse(p < .05, 1, 0 ) ) %>% # flag findings significant on 5% level
    pivot_wider( id_cols = c("decl","type","sim","meth"), values_from = sig, names_from = var ) %>% # pivot again such that each potential predictor has its own column with indicator of a "hit" per simulation/method pairs
    mutate( false_positives = rowSums( across( starts_with("p") ), na.rm = T ) ) %>% # calculate number of "hits"/false positives across variables for each simulation/method pair
    dplyr::select( decl, type, meth, false_positives ) %>% # keep only variables of interest for the plot
    table() %>% as.data.frame() %>% # prepare a table of frequencies of false hits per method
    # prepare names
    mutate(
      Method = case_when( meth == "two_step" ~ "Two-step procedure", meth == "lasso" ~ "Bayesian Lasso" ),
      Covariance = case_when( type == "nocov23" ~ "Independent predictors (k = 23)",
                              type == "yocov23" ~ "Covaried predictors (k = 23)",
                              type == "nocov7" ~ "Independent predictors (k = 7)") %>%
        factor( levels = c( paste0("Independent predictors (k = ", c(7,23),")" ),"Covaried predictors (k = 23)"), ordered = T )
    )
) %>% do.call( rbind.data.frame, . ) %>%
  mutate( Decline = paste0( ifelse( decl == "-0.3", "Small to moderate", "No"), " decline (b1 = ", decl, ")" ) )



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
  facet_grid( Decline ~ Covariance ) # prepare subplots

# save it
for ( i in forms ) ggsave( paste0("figures/Fig S2 simulation false positives under null",i ), dpi = 300, width = 13.1, height = 7.79 )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/sims.txt" )
