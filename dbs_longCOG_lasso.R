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


# ---- pre-processing  ----

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

# set-up priors (using the same priors for both models)
p <- c(
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
for ( i in names(f) ) m[[i]] <- brm_multiple( formula = f[[i]], family = student(), prior = p,
                                              data = df, sample_prior = T, seed = s, chains = ch,
                                              iter = it, warmup = wu, control = list( adapt_delta = ad ),
                                              file = paste0( "models/",i,".rds" ),
                                              save_model = paste0("models/",i,".stan")
                                              )



# ---- model comparisons ----

# clean the environment
rm( list = ls()[ !( ls() %in% c( "d", "df", "imp", "m", "scl", "doms", "tests", "s", "cbPal", "var_nms" ) ) ] )
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

# compute loo comparisons
l_comp <- lapply( 1:imp, function(i) loo_compare( l$m1_lasso_doms[[i]], l$m2_lasso_tests[[i]] ) )

# prepare a table for loo comparisons
t_comp <- sapply(
  1:imp, function(i)
    cbind(
      # negative = m1_lasso_doms wins, positive = m2_lasso_tests wins
      if ( l$m1_lasso_doms[[i]]$estimates[1,1] == l_comp[[i]][,"elpd_loo"][1] ) l_comp[[i]][2,"elpd_diff"] else -l_comp[[i]][2,"elpd_diff"],
      l_comp[[i]][2,"se_diff"]
      )
  ) %>% t() %>% as.data.frame() %>%
  `colnames<-`( c("elpd_dif","se_dif") ) %>%
  rownames_to_column( "dataset" ) %>%
  mutate( ci_low = qnorm( .025, elpd_dif, se_dif ),
          ci_upp = qnorm( .975, elpd_dif, se_dif ),
          sig_dif = ifelse( ci_low > 0 | ci_upp < 0, "*", "" ),
          winner = ifelse( elpd_dif < 0, "domains", "tests" )
          )

# save the table (the "G" in "Tab G." means "github only", i.e. the table won't be presented in the article)
write.table( t_comp %>% mutate_if( is.numeric, ~ round(.,2) %>% sprintf("%.2f",.) ),
             file = "tables/Tab G1 model comparisons.csv", sep = ",", row.names = F )


# ---- soft model checking ----

# check the highest Rhat for chains convergence and the highest Pareto-k for influential outliers
cbind.data.frame(
  Rhat = sapply( names(m) , function(i) max(m[[i]]$rhats, na.rm = T ) ) %>% round(3),
  Pareto_k = sapply( names(l) , function(i) sapply( 1:imp, function(j) max( l[[i]][[j]]$diagnostics$pareto_k ) ) ) %>% apply(., 2, max ) %>% round(3)
)

# print all Pareto-k values for both models
f.pk1 <- sapply( names(l), function(i) sapply( 1:imp, function(j) l[[i]][[j]]$diagnostics$pareto_k ) ) %>%
  as.data.frame() %>%
  pivot_longer( cols = everything(), values_to = "Pareto-k", names_to = "Predictors:",
                names_transform = function(x) { ifelse( x == "m1_lasso_doms", "cognitive functions", "cognitive tests" ) }
                ) %>%
  ggplot( aes(x = `Pareto-k`, fill = `Predictors:`, color = NULL ) ) +
  geom_histogram( position = "identity" ) +
  scale_fill_manual( values = alpha( cbPal[c(7,6)], alpha = .35 ) ) +
  scale_y_continuous( name = "Frequency" ) +
  geom_vline( xintercept = 0.5, linetype = "dashed", color = "red" ) +
  theme( legend.position = "bottom" )

# print the high tail of all Pareto-k values
f.pk2 <- sapply( names(l), function(i) sapply( 1:imp, function(j) l[[i]][[j]]$diagnostics$pareto_k ) ) %>%
  as.data.frame() %>%
  pivot_longer( cols = everything(), values_to = "Pareto-k", names_to = "Predictors:",
                names_transform = function(x) { ifelse( x == "m1_lasso_doms", "cognitive functions", "cognitive tests" ) }
                ) %>%
  ggplot( aes(x = `Pareto-k`, fill = `Predictors:`, color = NULL ) ) +
  geom_histogram( position = "identity" ) +
  scale_fill_manual( values = alpha( cbPal[6], alpha = .35 ) ) + # because there are only "cognitive tests" observations with k > 0.5
  geom_vline( xintercept = 0.5, linetype = "dashed", color = "red" ) +
  scale_y_continuous( name = "Frequency" ) +
  scale_x_continuous( limits = c(.5,.66) ) +
  theme_classic( base_size = 12 ) +
  theme( legend.position = "none" )

# create a single plot with an inset
f.pk1 + annotation_custom( ggplotGrob(f.pk2), xmin = .3, xmax = .7, ymin = 3.2e3, ymax = 7.2e3 )

# save as Fig. S6
ggsave( "figures/Fig S6 psis-loo pareto-ks.jpeg", dpi = 600, width = 9.64, height = 6.54 )


# ---- posterior predictive checks ----

# simulate values for each patient each six months from 2 years before to 12 years after surgery
n_seq = 24

# re-format id in d4 to a factor
d$id <- factor( d$id, levels = unique(d$id) )

# prepare data sets for prediction for each subject in the data set
d_seq <- expand.grid( seq( from = -2, to = 12, length.out = n_seq ), levels(d$id) ) %>%
  `colnames<-` ( c("time_y","id") ) %>%
  mutate( time = time_y + scl$Md$time )

# add cognitive domains and tests columns
for ( i in c(doms,tests) ) d_seq[[i]] <- NA

# add subjects' median (w.r.t. imputations) cognitive profile
for ( i in 1:nrow(d_seq) ) {
  for ( j in c(doms,tests) ) {
    d_seq[i, j] <- sapply( 1:imp,
                           function(k)
                             df[[k]][ df[[k]]$id == d_seq[i,"id"] & df[[k]]$time < (0+scl$Md$time) , j ]
                           ) %>% median()
  }
}

# prepare a list for predictions stratified by subjects chunks
# calculating prediction of single data points based on fixed-effects, random-effects and remaining distributional (residual) parameters
ppred <- list()
gc()

# add predictions to ppred
for (i in names(m) ) {
  for ( j in unique(d_seq$id) ) ppred[[i]][[j]] <- d_seq[ d_seq$id == j, ] %>%
      add_predicted_draws( m[[i]], seed = s ) %>%
      mutate(drs = scl$M$drs + scl$SD$drs * .prediction) %>%
      median_hdi(.width = .95) %>%
      mutate(drs.upper = ifelse(drs.upper > 144, 144, drs.upper) ) # manual censoring
} # took 6558.316 sec in total

# save the original prediction file as .rds
saveRDS( ppred, "data/ppred.rds" )

# collapse predictions to a single table for each model
for ( i in names(ppred) ) ppred[[i]] <- do.call( rbind.data.frame, ppred[[i]] ) %>%
  add_column( `Predicted by:` = ifelse(i == "m1_lasso_doms", "cognitive functions", "cognitive tests") )

# collapse one more to get a single prediction data set for both models
ppred <- do.call( rbind.data.frame, ppred ) %>% as.data.frame() %>% `rownames<-`( 1:nrow(.) )

# re-code id in the data set and ppreds such that the they are anonymized
all( levels(d$id) == levels(ppred$id) ) # TRUE, continue
N <-  levels(d$id) %>% length() # extract the number of subjects
levels(ppred$id) <- paste0( "S", sprintf("%.3d", 1:N) )
levels(d$id) <- paste0( "S", sprintf("%.3d", 1:N) )

# save the original prediction file as .rds
write.table( ppred, "data/ppred.csv", sep = ",", row.names = F )

# plot predictions and observed data for each subject separately
d %>% mutate( drs = drs_tot ) %>%
  filter( complete.cases(drs) ) %>%
  # plotting proper
  ggplot( aes(x = time_y, y = drs) ) +
  # add prediction lines and 95% compatibility intervals
  geom_line( data = ppred, size = 1, aes(x = time_y, y = drs, group = `Predicted by:`, color = `Predicted by:` ) ) +
  geom_ribbon( data = ppred, alpha = .2, aes(x = time_y, y = drs, ymin = drs.lower, ymax = drs.upper, group = `Predicted by:`, fill = `Predicted by:` ) ) +
  # add observed points
  geom_point( color = "black", size = 1.25 ) +
  # finish the plot with the last cosmetic changes
  scale_color_manual( values = cbPal[c(2,3)] ) +
  scale_fill_manual( values = cbPal[c(2,3)] ) +
  facet_wrap( ~id, nrow =  14 ) + # arrange to a 14 x 9 grid
  labs( x = "Time after surgery (Years)", y = "DRS-2 (0-144 points)") +
  theme_classic( base_size = 7 ) +
  theme( legend.position = "bottom" )

# save it as Fig. S5
ggsave( "figures/Fig S5 posterior predictive check.jpeg", dpi = 600, width = 7.77, height = 11.19 )


# ---- models' posteriors ----

# extract and summarize posteriors of the "fixed-effects"
post <- lapply( names(m), function(i) m[[i]] %>%
                  spread_draws( `b_.*`, regex = T ) %>% # extract all "fixed-effect" parameters
                  select( contains("b_") ) %>% # get rid of the info about chains
                  # re-scale to DRS-2 raw scale
                  mutate( b_Intercept = ( b_Intercept * scl$SD$drs ) + scl$M$drs,
                          across( !b_Intercept, function(x) { x * scl$SD$drs } )
                          ) %>%
                  # calculate median posterior, 95% PPI and probability of being negative
                  apply( . , 2 , function(x) { c(b = median(x), PPI = hdi(x, .width = .95), pd = sum(x<0)/length(x) ) } ) %>%
                  # relocate time such that it is after the baseline effects
                  as.data.frame() %>%  relocate( b_time, .after = ifelse( i == "m1_lasso_doms", "b_visp_wm", "b_sc_staix2") ) %>%
                  # tidy up
                  t() %>% as.data.frame() %>% rownames_to_column( var = "Parameter" )
                ) %>% `names<-` ( names(m) )

# the plots
f <- lapply( names(post), function(i) post[[i]] %>%
               # select only the estimands (i.e., the interaction terms and the global time effect for comparison)
               slice( which( grepl("time", Parameter) ) ) %>%
               # rename the Parameter column according to the var_nms file such that the predictors are
               # named in accordance with the in-text math model (i.e., delta[predictor] )
               mutate( Parameter = sub( "time:", "", Parameter ) %>% sapply( ., function(x) var_nms[x, ] ) ) %>%
               # plotting proper
               ggplot( aes( x = reorder(Parameter, b, decreasing = T), y = b, ymin = PPI1, ymax = PPI2 ) ) +
               geom_pointrange( shape = 21, size = 4, fatten = 2,
                                fill = alpha( case_when( i == "m1_lasso_doms" ~ cbPal[7], i == "m2_lasso_tests" ~ cbPal[6] ), alpha = 1 ),
                                color = alpha( case_when( i == "m1_lasso_doms" ~ cbPal[7], i == "m2_lasso_tests" ~ cbPal[6] ), alpha = .3 )
                                ) +
               geom_hline( yintercept = 0, size = 1, linetype = "dashed", color = "black" ) +
               labs( x = NULL, y = "DRS-2 (points per year)") +
               scale_y_continuous( limits = c(-1,1) ) +
               scale_x_discrete( labels = function(x) parse( text = paste0( "delta[", x, "]" ) ) ) + # add predictor parameter names
               theme( axis.text = element_text( size = 16 ) ) +
               coord_flip()
               )

# put the plots side-to-side
( f[[2]] | f[[1]] ) + plot_annotation( tag_levels = "a" ) & theme( plot.tag = element_text(face = "bold") )

# save it
ggsave( "figures/Fig 3 estimands.jpg" , dpi = 600, width = 1.3*9.64, height = 11.7 )

# prepare the Tab. S2 (cog. tests posterior summary) and Tab. S3 (cog. functions posterior summary)
t <- lapply( names(post), function(i) post[[i]] %>%
               # tidy the table columns 
               mutate( Parameter = paste0( "\t", var_nms[Parameter,1] ),
                       b = sprintf( "%.2f", round(b,2 ) ),
                       `95% PPI` = paste0( "[", sprintf("%.2f",round(PPI1,2)),", ",sprintf("%.2f",round(PPI2,2)),"]" ),
                       `Pr(b < 0)` = sprintf( "%.3f", round(pd,3) )
                       ) %>%
               # add rows delineating types of parameters
               add_row( Parameter = "Global intercept (α)", .before = 1 ) %>%
               add_row( Parameter = "Baseline correlates (β) ", .before = 3 ) %>%
               add_row( Parameter = "Time-dependent effects (𝛿)", .before = ifelse( i == "m1_lasso_doms", 11, 27 ) ) %>%
               # keep only variables of interest
               select( Parameter, b, `95% PPI`, `Pr(b < 0)` )
             ) %>% `names<-`( c("S3","S2") )

# save the tables as .csv
for ( i in names(t) ) write.table( t[[i]], sep = ";", row.names = F, na = "", quote = F,
                                   file = paste0("tables/Tab ",i," summary of posteriors (cognitive ",
                                                 ifelse(i=="S3","functions","tests"),").csv" )
                                   )


# ---- stats for in-text reporting ----

# summary of models comparisons
table(t_comp[ , c("sig_dif","winner") ] )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/lasso.txt" )
