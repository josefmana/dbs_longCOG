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
theme_set( theme_classic(base_size = 24) )

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


# ---- comparing estimands ----

# extract posteriors from all models
draws <- list()

# loop through the models to fill-in the draws list
for ( i in c("m1_lasso_doms","m2_lasso_tests","m3_doms_cov","m4_tests_cov") ) {
  m <- readRDS( paste0(getwd(), "/models/", i, ".rds") )
  draws[[i]] <- as_draws_df(m) %>%
    `colnames<-` ( gsub("_drs","",names(.) ) ) %>% # rename parameters in covariate models to appropriate format
    select( starts_with( c("b_","bsp_") ) & !starts_with( c("b_bdi","b_led","bsp_bdi") ) ) # keep only fixed-effects
  rm(m)
  gc()
}

# pull all posteriors into a single file
post <- lapply( names(draws), function(i)
  
  # prepare a column with model labels
  draws[[i]] %>%
    # back-transform posteriors to raw DRS-2 (or DRS-2 per year) scale
    mutate( b_Intercept = ( b_Intercept * scl$SD$drs ) + scl$M$drs,
            across( !b_Intercept, function(x) { x * scl$SD$drs } )
            ) %>%
    # calculate median posterior and 95% PPI
    apply( . , 2 , function(x) { c( b = median(x), PPI = hdi(x, .width = .95) ) } ) %>% t() %>%
    as.data.frame() %>% rownames_to_column( var = "Parameter" ) %>% add_column( `Predicted by:` = i )
  
  ) %>%
  
  # pull both "cognitive functions" and "cognitive tests" models together
  do.call( rbind.data.frame, . ) %>%
  
  # make final adjustments
  mutate(
    # re-code model to a nicer and consistent form (will serve to color density plots)
    `Predicted by:` = case_when(
      grepl("m1_lasso_doms", `Predicted by:`) ~ "Cognitive functions",
      grepl("m2_lasso_tests", `Predicted by:`) ~ "Cognitive tests",
      grepl("m3_doms_cov", `Predicted by:`) ~ "Cognitive functions (with covariates)",
      grepl("m4_tests_cov", `Predicted by:`) ~ "Cognitive tests (with covariates)"
      ),
    # add labels for effects grouping (will serve to split facets of the plot)
    Group = ifelse( Parameter == "b_Intercept", "alpha",
                    ifelse( grepl( "time|bdi|led|age", Parameter ), "delta", "beta" )
                    )
    )

# prepare a list for subplots of Fig. S5 (cognitive functions) and Fig. S6 (cognitive tests)
fig <- list()

# loop through type of predictor (functions vs. tests)
for ( i in c("functions","tests") ) fig[[i]] <- lapply(
  
  # loop through the parameter type
  # (need to use lapply instead of for loop here to ensure the correct Greek symbols are printed,
  # the for loop always printed the last one used, i.e. delta, in all plots for some reason)
  c("alpha","beta","delta"), function(j)
    
    # select only relevant parameters
    post[ post$Group == j & grepl(i,post$`Predicted by:`), ] %>%
    
    # prepare the labels
    mutate(
      title = case_when( j == "alpha" ~ "Global intercept",
                         j == "beta" ~ "Baseline correlates",
                         j == "delta" ~ "Time-dependent effects"),
      Parameter = case_when(
        j == "delta" ~ sub( "time:", "", Parameter ) %>% sapply( ., function(x) var_nms[x, ] ),
        j != "delta" ~ Parameter %>% sapply( ., function(x) var_nms[x, ] ) )
    ) %>%
    
    # plotting proper
    ggplot( aes( x = reorder(Parameter, b, decreasing = T), y = b, ymin = PPI1, ymax = PPI2, fill = `Predicted by:` ) ) +
    geom_pointrange( shape = 21, size = 3, fatten = 2, position = position_dodge( width = .7 ), color = cbPal[1] ) +
    geom_hline( yintercept = 0, size = 1.5, linetype = "dashed", color = "red" ) +
    
    # scale the axes and colors
    scale_x_discrete( labels = ifelse( j == "alpha", parse( text = "alpha"), function(x) parse( text = paste0( j, "[", x, "]" ) ) ) ) +
    scale_y_continuous( limits = case_when( j == "alpha" ~ c( 139,141 ) ) ) +
    scale_fill_manual( values =  case_when( i == "functions" ~ cbPal[c(7,2)], i == "tests" ~ cbPal[c(6,3)] ) ) +
    
    # label and flip axes, position legend and add title
    labs( x = NULL, y = ifelse( j == "delta", "DRS-2 (points per year)", "DRS-2" ) ) +
    theme( axis.text = element_text( size = 24 ), legend.position = "bottom" ) +
    facet_wrap( . ~ title) +
    coord_flip()
    
) %>% `names<-`(c("alpha","beta","delta") ) # after the lapply loop is finished, name each plot accordingly

# arrange and print the cognitive tests plot
with( fig$tests, (alpha / beta / delta ) + plot_layout( heights = c(1,23,27), guides = "collect") & theme( legend.position = "bottom") )
ggsave( "figures/Fig S3 cognitiove tests covariate check.jpg" , dpi = 600, width = 1.2*9.64, height = 4.2*6.54 )

# arrange and print the cognitive functions plot
with( fig$functions, (alpha / beta / delta ) + plot_layout( heights = c(1,7,11), guides = "collect") & theme( legend.position = "bottom") )
ggsave( "figures/Fig S4 cognitive functions covariate check.jpg" , dpi = 600, width = 1.2*9.64, height = 3.3*6.54 )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/covar.txt" )
