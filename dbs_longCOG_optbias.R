# A neat small script that checks for overoptimistic bias in the best case scenario for choosing significant predictors of an
# outcome from ten independent predictors via two-stage (same data) procedure that first pre-selects potential predictors via
# simple correlations and then estimates effect sizes via multiple regressions.

# list packages to be used
pkgs <- c("rstudioapi", "dplyr", "tidyverse", "MASS", "brms", "bayestestR" )

# load or install each of the packages as needed
for (i in pkgs) {
  if (i %in% rownames (installed.packages()) == F) install.packages(i) # install if it ain't installed yet
  if ( i %in% names (sessionInfo() $otherPkgs) == F ) library(i, character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- prepare data-generating function ----

# prepare a function with varying number of subjects, number of predictors, varying pre-selection threshold and effect size
# N = number of subjects
# n_prds = number of independent predictors (integer)
# eff = list of effect sizes (vector)
# cor = covariance (identical with correlation in case of standard normals) matrix
# thres = p-value used a threshold for pre-selection (numeric)
# meth = either two-step procedure ("two-step"), multiple regression ("multi") or Bayesian lasso ("lasso")
sim <- function( N = 126, n_prds = 10, eff = rep(0,10), cor = diag( rep(1,10) ), thres = .1, meth = "two-step" ) {
  
  # generate data and pre-select variables
  t <- mvrnorm( N, rep(0,n_prds), cor ) # simulate "n_prds" preedictors from multivariate standard normal distribution
  o <- rnorm( N, rowSums(t[ , 1:length(eff) ] * eff) ) # generate the outcome
  sel <- sapply( 1:n_prds, function(i) ifelse( cor.test( t[,i], o )$p.value < thres, 1, 0 ) ) # pre-select variables
  
  # compute the final multiple linear regression
  if( meth == "two-step" ) {
    
    # if the two-step procedure was selected, used only the preselected variables
    if( sum(sel) == 0 ) reg <- lm( o ~ 1 ) # intercept only if no variable was selected
    else reg <- lm( o ~ t[ , which(sel == 1) ] )
    
    # extract final p-values
    if ( sum(sel) == 0 ) p <- NULL
    else p <- as.vector( summary(reg)$coefficients[ 2:(sum(sel)+1), 4 ] ) # raw p-values
    
  } else if( meth == "multi" ) {
    
    # otherwise regress the outcome on all pre-surgery variables
    reg <- lm( o ~ t )
    p <- as.vector( summary(reg)$coefficients[ -1 , 4 ] )
    
  } else if ( meth == "lasso" ) {
    
    # finally a Bayesian lasso
    reg <- brm( as.formula( paste0( "o ~ ", paste( paste0("t",1:n_prds), collapse = " + " ) ) ),
                prior = prior( lasso(1), class = b ),
                data =  cbind.data.frame(o, t) %>% `colnames<-`( c("o", paste0("t",1:n_prds) ) )
                )
    
    # extract probability of direction and recalculate to a p-value equivalent
    p <- pd_to_p( p_direction(reg)$pd, direction = "two-sided" )[ 2:(n_prds+1) ]
    
  }
  
  # sum-up the results and print them
  res <- rep(1, n_prds) # prepare ones (for those that were not tested, techniqually p = 1, easier to work with compare to NAs)
  if( meth == "two-step" ) res[ which(sel == 1) ] <- p
  else res <- p # add p-value
  return(res) # return the results
  
}

# prepare two types of covariance matrixes
nocor <- function(K = 10) return( diag( rep(1,K) ) )
nocov <- nocor(K = 27)
yocov <- matrix( c(
  1,rep(.7,2),rep(0,24), .7,1,.7,rep(0,24), rep(.7,2),1,rep(0,24), # Trail Making Test
  rep(0,3),1,rep(0,23), # BNT
  rep(0,4),1,rep(.7,2),rep(0,20), rep(0,4),.7,1,.7,rep(0,20), rep(0,4),rep(.7,2),1,rep(0,20), # verbal memory
  rep(0,7),1,.7,rep(0,18), rep(0,7),.7,1,rep(0,18), # nonverbal memory
  rep(0,9),1,rep(0,17), # RCFT
  rep(0,10),1,.7,rep(.6,2),rep(.5,2),rep(0,11), rep(0,10),.7,1,rep(.6,2),rep(.5,2),rep(0,11), # Stroop color
  rep(0,10),rep(.6,2),1,.7,rep(.6,2),rep(0,11), rep(0,10),rep(.6,2),.7,1,rep(.6,2),rep(0,11), # Stroop words
  rep(0,10),rep(.5,2),rep(.6,2),1,.7,rep(0,11), rep(0,10),rep(.5,2),rep(.6,2),.7,1,rep(0,11), # Stroop color-word
  rep(0,16),1,.7,rep(0,9), rep(0,16),.7,1,rep(0,9), # verbl fluency
  rep(0,18),1,rep(0,8), # BDI
  rep(0,19),1,rep(.7,2),rep(0,5), rep(0,19),.7,1,.7,rep(0,5), rep(0,19),rep(.7,2),1,rep(0,5), # age and duration
  rep(0,22),1,rep(0,4), # education
  rep(0,23),1,rep(.7,2),0, rep(0,23),.7,1,.7,0, rep(0,23),rep(.7,2),1,0, # MDS-UPDRS III
  rep(0,26),1 # LEDD
), ncol = 27, byrow = T )


# ---- simulating the null ----

# prepare a table that shows how false positive error rates depend on the threshold for case with N = 103 (as in Kim et al., 2014)
# and K = 10 tests
d0 <-list()

# loop through covariance structures and methods of estimation
for ( i in c("nocov","yocov") ) {
  for ( j in c("two-step","multi","lasso") ) {
    
    d0[[i]][[j]] <- sapply( 1:1e2, function(k)
      sim( N = 103, n_prds = 27, eff = rep(0,27), cor = get(i), thres = .2, meth = j )
      ) %>% t() %>% as.data.frame()
    
  }
}

# summarize false negatives
fn0 <- sapply( names(d0), function(i) sapply( names(d0[[i]]), function(j) ifelse( d0[[i]][[j]] < .05, 1, 0 ) %>% rowSums() %>% table() ) )

# the results seem to be almost invariant to threshold used, overall, in ca 2/5 of cases there is at least one false positive

# Under the idealized conditions of this simulation, there is about 11 % chance of getting Kim et al.'s (2014) results of
# identifying (exactly) three significant predictors even when the null is true


# ---- simulating effects ----

# simulate a few thousands of random data sets with null effects
d1 <- lapply( c("two-step","multi","lasso"), function(i)
  sapply( 1:1e2, function(j) sim( N = 103, n_prds = 27, eff = rep(0.3,5), thres = .2, type = i ) ) %>% t() %>% as.data.frame()
)

# summarize false negatives
sapply( 1:3, function(i) ifelse( ( d1[[i]][ ,6:27] > .05 | is.na(d1[[i]][ ,6:27]) ), 0, 1 ) %>% rowSums() %>% table() )

# summarize true positives
sapply( 1:3, function(i) ifelse( ( d1[[i]][ ,1:5] > .05 | is.na(d1[[i]][ ,1:5]) ), 0, 1 ) %>% rowSums() %>% table() )
