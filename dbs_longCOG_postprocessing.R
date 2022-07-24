# All analyses reported in the article (Mana et al., in review) ran in R version 4.2.0 (2022-04-22),
# on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.4.

# I used the following versions of packages employed: dplyr_1.0.9, tidyverse_1.3.1, psych_2.2.5,
# brms_2.17.0, tidybayes_3.0.2, ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "tidyverse", # for pivot_longer
  "psych", # for EFA
  "brms", # for Bayesian model fitting / interface with Stan
  "tidybayes", # for posteriors manipulations
  "ggplot2", # for plotting
  "patchwork" # for ggplots manipulations
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

# set ggplot theme
theme_set( theme_classic(base_size = 25) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# set number of multiple imputations used to account for missing pre-surgery data
imp = 100

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# read the raw data set as well as imputed baseline data sets
d1 <- readRDS("data/long_dat.rds")$d1
d3 <- readRDS("data/long_dat.rds")$d3
scl <- readRDS("data/long_dat.rds")$scl

# read a file with mapping variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")

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

# ----------- baseline cognitive profile post-processing  -----------

# read all EFA models
efa <- readRDS( "models/efa.rds")

# prepare a table for performance indexes of each solution for each imputed data set
fat <- array(
  # create an empty 6 (factor solutions) x 5 (performance indexes) x 100 (imputations) array
  data = NA, dim = c( length( efa[[1]] ), 5, imp ),
  # add dimension names
  dimnames = list(
    paste0( 1:length(efa[[1]])+2, "-factor" ), # number of factors
    c("TLI", "RMSEA", "RMSEA_90_CI_low", "RMSEA_90_CI_upp", "var_account"), # performance indexes
    1:imp # number of imputed data set
  )
)

# loop through the array and gather all performance indexes
for ( i in 1:imp ) {
  for ( j in 1:length(efa[[i]]) ) {
    fat[ j , , i ] <- c(
      efa[[i]][[j]]$TLI, # Tucker-Lewis index, > 0.9 is considered good
      efa[[i]][[j]]$RMSEA[1], # root-mean square error of approximation, < 0.08 is considered good
      efa[[i]][[j]]$RMSEA[2], # 90% CI RMSEA, lower boundary
      efa[[i]][[j]]$RMSEA[3], # 90% CI RMSEA, upper boundary (ideally should be < 0.08)
      max( efa[[i]][[j]]$Vaccounted["Cumulative Var",] ) # total variance "accounted for" by all included factors (i.e., the highest one from "Cumulative Var")
    )
  }
}

# ----------- Tab S1 summary of performance index of factor analyses  -----------

# summarize the fat table
t.s1 <- data.frame( Model = paste0( 3:8, "-factor"), TLI = NA, RMSEA = NA, RMSEA_90_CI_upp = NA, var_account = NA )

# fill-in summary of each model into t.s1
for ( i in t.s1$Model ) {
  for ( j in names(t.s1)[-1] ) {
    t.s1[ which(t.s1$Model == i) , j ] <- paste0(
      sprintf( "%.3f" , round( mean( fat[ i, j, ] ), 3 ) ), " (",
      sprintf( "%.3f" , round( sd( fat[ i, j, ] ), 3 ) ), ")"
    )
  }
}

# prepare a list for RMSEA and TLI frequencies across all imputations
fat_perc <- list(
  # upper RMSEA lower than 0.08
  `upper bound RMSEA < 0.08 (%)` = t( fat[ , "RMSEA_90_CI_upp" , ] ) %>%
    as.data.frame %>%
    mutate( across( everything() , ~ replace( . , . > .08 , NA) ) ),
  # TLI higher than 0.9
  `TLI > 0.90 (%)` = t( fat[ , "TLI" , ] ) %>%
    as.data.frame %>%
    mutate( across( everything() , ~ replace( . , . < .9 , NA) ) )
)

# calculate the percentages
for ( i in names(fat_perc) ) fat_perc[[i]] <- sapply(
  names(fat_perc[[i]]), function(j) 100 * sum( complete.cases( fat_perc[[i]][[j]] ) ) / imp
)

# bind the percentages to the FA summary table (t.s1)
t.s1 <- t.s1 %>% left_join(
  do.call( cbind.data.frame , fat_perc ) %>%
    # prepare a column "model" and bind the percentages to the t.s1
    rownames_to_column( var = "Model")
)

# change column names that are not fancy enough
colnames(t.s1)[4:5] <- c( "RMSEA 90% CI (upper bound)", "Total variance accounted for")

# save as csv
write.table( t.s1, file = "tables/Tab S1 factor analysis performance indexes.csv", sep = ",", row.names = F, quote = F)


# ----------- Fig S2 factor analyses performance indexes  -----------

# set-up a list to contain Fig S2 component figures
f.s2 <- list()

# loop through TLI and upper 90% CI RMSEA
for ( i in c("TLI","RMSEA_90_CI_upp") ) {
  f.s2[[i]] <- fat[ , i , ] %>%
    # format the table for plotting
    t %>% as.data.frame %>%
    pivot_longer( everything() , names_to = "Model" , values_to = i ) %>%
    rename( "index" = i ) %>%
    # plotting proper
    ggplot( aes( x = index , fill = Model) ) +
    geom_density( alpha = .4 , color = NA ) +
    scale_fill_manual( values = cbPal[3:8] ) + # use colorblind-friendly palette
    geom_vline( # add a vertical line depicting good performance heuristics
      linetype = "dashed", size = 1.2, xintercept = case_when( i == "TLI" ~ .9, i == "RMSEA_90_CI_upp" ~ .08 )
    ) +
    labs( y = "Density", x = case_when( i == "TLI" ~ "TLI", i == "RMSEA_90_CI_upp" ~ "RMSEA (upper 90% CI)" ) ) +
    scale_x_continuous( limits = case_when( i == "TLI" ~ c(0.6,1.05), i == "RMSEA_90_CI_upp" ~ c(0.03, 0.12) ),
                        breaks = case_when( i == "TLI" ~ seq(.6,1,.1), i == "RMSEA_90_CI_upp" ~ seq(.04,.12,.02) ),
                        labels = case_when(
                          i == "TLI" ~ sprintf( "%.1f", round( seq(.6, 1, .1), 1) ),
                          i == "RMSEA_90_CI_upp" ~ sprintf( "%.2f", round( seq(.04,.12,.02), 2) )
                        ))
}

# arrange Fig S2 for printing
f.s2$TLI / f.s2$RMSEA_90_CI_upp + plot_layout( guides = "collect" ) +
  plot_annotation( tag_levels = "a" ) & theme( plot.tag = element_text(face = "bold") )

# save as Fig S2
ggsave( "figures/Fig S2 factor analysis performance indexes.jpg", height = 1.25 * 6.07, width = 1.25 * 11.5, dpi = 300 )


# ----------- Tab 2 factor loadings  -----------

# the 7-factor solution was chosen
nf = 7

# prepare an array for loading matrices of each imputed EFA
loads <- array(
  # create an empty 25 ( 23 tests + 2 variance accounted) x 7 (factors) x 100 (imputations) array
  data = NA, dim = c(25, nf, imp),
  dimnames = list( c( rownames(efa[[1]][[nf-2]]$loadings) , "Proportion Var", "Cumulative Var" ), doms, 1:imp )
)

# fill-in all loadings from efa objects prepared above
for( i in 1:imp ) loads[ , , i ] <- efa[[i]][[nf-2]]$loadings %>%
  as.data.frame %>%
  bind_rows( efa[[i]][[nf-2]]$Vaccounted[ "Proportion Var", ] ) %>%
  bind_rows(
    apply( efa[[i]][[nf-2]]$Vaccounted[ "Proportion Var", ] %>% t, 1 , cumsum ) %>% t %>% as.data.frame()
  ) %>% as.matrix()

# write a summary of loadings across all imputed data sets
t2 <- matrix(
  data = NA, nrow = nrow(loads[ , , 1]), ncol = ncol(loads[ , , 1]),
  dimnames = list( rownames(loads[ , , 1 ]) , colnames(loads[ , , 1]) )
) %>% as.data.frame()

# fill-in averages and SDs across all imputations 
for ( i in rownames(t2) ) {
  t2[ i , ] <- paste0(
    sprintf( "%.2f" , round( loads[i, , ] %>% t %>% colMeans , 2 ) ), " (", # mean
    sprintf( "%.2f" , round( loads[i, , ] %>% t %>% apply( . , 2 , sd ), 2) ), ")" # SD
  )
}

# make rownames to a column
t2 <- t2 %>% rownames_to_column( var = "Test" )

# rename tests (rows) and factors (columns) such that they are publication-ready
for ( i in 1:(nrow(t2)-2) ) t2[i,"Test"] <- var_nms[ t2[i,]$Test,  ]
for ( i in 2:ncol(t2) ) colnames(t2)[i] <- var_nms[ colnames(t2)[i],  ]

# save as csv
write.table( t2, file = "tables/Tab 2 factor loadings.csv", sep = ",", row.names = F, quote = F )


# ----------- Fig S4 prior predictive check  -----------

# including prior predictive check to "post-processing" to keep the "statistical modelling" file clean
# set-up the linear model (doing PPC for the primary model only, ie, m1_nocov,  the rest should be "centered"
# around it as they are the same model with added parameters with zero-centered priors)
f.drs <- paste0( "drs | cens(cens_drs) ~ 1 + ", paste( "time", doms, sep = " * " , collapse = " + " ),
                 " + (1 + time || id)" ) %>% as.formula %>% bf

# set-up priors
p <- c(
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
  # other distributional parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu )
)

# fit the model sampling only from the prior
m <- brm_multiple(
  formula = f.drs, family = student(), prior = p, data = d3, sample_prior = "only",
  seed = 87542, chains = 4, iter = 200, warmup = 100, control = list( adapt_delta = .95 ), 
  file = "models/m1_prior_prediction.rds", save_model = "models/m1_prior_prediction.stan"
)

# prepare data at which PPC will be evaluated
d_seq <- expand.grid( time_y = c(-2,12), id = d1$id, proc_spd = 0, epis_mem = 0, verb_wm = 0,
                      visp_mem = 0, set_shift = 0, anxiety = 0, visp_wm = 0) %>%
  mutate( time = time_y + scl$Md$time) %>%
  add_epred_draws(m) # add predicted mean response (DRS-2) 

# randomly choose which prior draws to visualize
set.seed(87542) # seed for reproducibility
it <- sample( 1:max( unique( d_seq$.draw ) ), size = 25, replace = F)

# plot twenty-five PPCs to a 5x5 grid
d_seq %>%
  slice( which(.draw %in% it) ) %>%
  mutate( drs = .epred * scl$SD$drs + scl$M$drs ) %>%
  ggplot( aes( x = time_y, y = drs, group = id) ) +
    geom_line( alpha = .1 ) +
    geom_hline( yintercept = 144 ) +
    geom_hline( yintercept = 138, linetype = "dashed" ) +
    labs( x = "Time (years)", y = "DRS-2 prior prediction" ) +
    scale_y_continuous( breaks = c(50,100,150), labels = c(50,100,150) ) +
    scale_x_continuous( limits = c(-2,12) ) +
    facet_wrap( ~ .draw ) +
    theme_bw( base_size = 18 ) +
    theme( strip.background = element_blank(),
           strip.text.x = element_blank(),
           panel.grid = element_blank()
           )

# save the ensuing plot as Fig S4
ggsave( "figures/Fig S4 prior predictive check.jpg", dpi = 600, width = 10.1, height = 11.2 )


# ----------- extract MCMC draws -----------

# remove unneeded objects
rm( list = ls()[ !( ls() %in% c( "cbPal", "d1", "d3", "doms", "imp", "scl", "var_nms" ) ) ] )

# extract posteriors from all models
draws <- list()

# loop through the models to fill-in the draws list
for ( i in c("m0_base","m1_nocov","m2_flat_priors","m3_wcov") ) {
  m <- readRDS( paste0(getwd(), "/models/", i, ".rds") )
  draws[[i]] <- posterior::as_draws_df(m)
  rm(m)
  gc()
}

# now prepare list with posteriors from each model
# make a list of variable names for all fixed-effects
pars <- list(
  intercept = names(draws$m1_nocov)[1], # global intercept
  base = names(draws$m1_nocov)[3:9], # baseline correlates
  time = names(draws$m1_nocov)[c(2,10:16)] # time-dependent effects
)

# prepare a list for posteriors of relevant fixed-effects
post <- list()

# loop through models to fill-in posteriors
for ( i in names(draws) ) post[[i]] <- draws[[i]] %>%
  `colnames<-` ( gsub("_drs","",names(.) ) ) %>% # rename parameters in m3_wcov to appropriate format
  select( which( names(.) %in% unlist(pars) ) )

# pull all posteriors into a single file
post <- do.call( rbind.data.frame , post ) %>%
  # prepare a column with model labels
  rownames_to_column( var = "Model" ) %>%
  # pivot such that all parameters are in a single long column
  pivot_longer(
    cols = unlist(pars, use.names = F), values_to = "DRS-2", names_to = "Parameter"
  ) %>%
  # make final adjustments
  mutate(
    # re-code model to a nicer and consistent form (will serve to color density plots)
    Model = case_when(
      grepl("m0_base", Model) ~ "Correlated varying effects",
      grepl("m1_nocov", Model) ~ "Uncorrelated varying effects",
      grepl("m2_flat_priors", Model) ~ "Flat priors",
      grepl("m3_wcov", Model) ~ "Covariates"
    ),
    # add labels for effects grouping (will serve to split facets of the plot)
    Group = case_when(
      Parameter %in% pars$intercept ~ "intercept",
      Parameter %in% pars$base ~ "base",
      Parameter %in% pars$time ~ "time"
    ),
    # back-transform posteriors to raw DRS-2 (or DRS-2 per year) scale
    `DRS-2` = case_when(
      Parameter == "b_Intercept" ~ `DRS-2` * scl$SD$drs + scl$M$drs, # intercept
      Parameter != "b_Intercept" ~ `DRS-2` * scl$SD$drs # slopes
    )
  )

# change character columns to ordered factors such that ggplot draws them in correct order
# need to reverse the order by rev(.) because ggplot counts from bottom on the y-axis
post$Model <- factor( post$Model , levels = unique(post$Model)[c(2,1,3,4)]  , ordered = T )
post$Parameter <- factor( post$Parameter , levels = rev( unlist(pars, use.names = F) ) , ordered = T )
post$Group <- factor( post$Group , levels = rev( unique(post$Group) ) , ordered = T )


# ----------- Fig 2 posterior primary model posteriors  -----------

# prepare list for facets of Fig. 2 (m1_nocov posteriors)
f2 <- list()

# loop through the three parameter groups
for ( i in rev( levels(post$Group) ) ) {
  f2[[i]] <- post[ post$Group == i & post$Model == "Uncorrelated varying effects" , ] %>%
    mutate(
      title = case_when( i == "intercept" ~ "Global intercept",
                         i == "base" ~ "Baseline correlates",
                         i == "time" ~ "Time-dependent effects")
    ) %>%
    ggplot( aes(y = Parameter, x = `DRS-2`, fill = stat(x < 0) ) ) +
    stat_slab( geom = "slab", size = 1,  alpha = .4 ) +
    coord_cartesian(
      ylim = case_when( i == "intercept" ~ c(1.3,1.5), i == "base" ~ c(1.3,7.5), i == "time" ~ c(1.3,8.5) )
    ) +
    geom_vline( xintercept = 0, linetype = "dashed" ) + # add line showing zero effect
    scale_fill_manual( values = cbPal[1:2]) +
    scale_x_continuous(
      limits = case_when( i == "intercept" ~ c(138, 143) ),
      name = case_when( i == "intercept" ~ "" , i == "base" ~ "DRS-2", i == "time" ~ "DRS-2 (points per year)")
    ) +
    scale_y_discrete( name = NULL, labels = rev( var_nms[ pars[[i]] , ] ) ) +
    facet_wrap( . ~ title) +
    theme(
      legend.position = "None",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 20),
      axis.text.y = element_text( vjust = -1.5 ),
      axis.ticks.y = element_blank()
    )
}

# put them together
( ( f2$intercept / f2$base + plot_layout( heights = c(1,7) ) ) | f2$time )

# save as Fig 2
ggsave( "figures/Fig 2 model posteriors.jpg", height = 1.3 * 8.53, width = 1.6 * 9.05, dpi = 300 )


# ----------- Fig S6 posterior comparison across models  -----------

# plot it one by one per group
f.s6 <- list()

# loop through the three parameter groups
for ( i in rev( levels(post$Group) ) ) f.s6[[i]] <- post[ post$Group == i , ] %>%
  mutate(
    title = case_when( i == "intercept" ~ "Global intercept",
                       i == "base" ~ "Baseline correlates",
                       i == "time" ~ "Time-dependent effects")
  ) %>%
  ggplot( aes(y = Parameter, x = `DRS-2`, fill = Model) ) +
  stat_slab( geom = "slab", size = 1,  alpha = .4 ) +
  coord_cartesian(
    ylim = case_when( i == "intercept" ~ c(1.3,1.5), i == "base" ~ c(1.3,7.5), i == "time" ~ c(1.3,8.5) )
  ) +
  geom_vline( xintercept = 0, linetype = "dashed" ) +
  scale_fill_manual( values = cbPal[c(7,8,6,4)]) +
  scale_x_continuous(
    limits = case_when( i == "intercept" ~ c(138, 143) ),
    name = case_when( i == "intercept" ~ "" , i == "base" ~ "DRS-2", i == "time" ~ "DRS-2 (points per year)")
  ) +
  scale_y_discrete( name = NULL, labels = rev( var_nms[ pars[[i]] , ] ) ) +
  facet_wrap( . ~ title) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text( vjust = -1.5 ),
    axis.ticks.y = element_blank()
  )

# put the facets for each parameter type to a single plot
( ( f.s6$intercept / f.s6$base + plot_layout( heights = c(1,7) ) ) | f.s6$time )  +
  plot_layout( guides = "collect" ) & theme( legend.position = "bottom" )

# save as Fig S6
ggsave( "figures/Fig S6 posteriors across models.jpg", height = 1.3 * 8.53, width = 1.6 * 9.05, dpi = 300 )


# ----------- Fig S5 prior vs posterior comparison in the primary model  -----------

# keep only m1_nocov posteriors and add priors of m1_nocov for Fig S5
samples <- list(
  # extract priors
  Prior = draws$m1_nocov %>% select( starts_with("prior_"), -contains("sd"), -prior_nu, -prior_sigma ) %>%
    `colnames<-` ( gsub("prior_","",names(.) ) ) %>% # rename such that columns are named as in pars
    rename( "b_Intercept" = "Intercept" ), # rename the Intercept manually
  # extract posteriors
  Posterior = draws$m1_nocov %>% select( which( names(.) %in% unlist(pars) ) )
)

# pull all samples from m1_nocov into a single file
samples <- do.call( rbind.data.frame , samples ) %>%
  # prepare a column with model labels
  rownames_to_column( var = "Distribution" ) %>%
  # pivot such that all parameters are in a single long column
  pivot_longer(
    cols = unlist(pars, use.names = F), values_to = "DRS-2", names_to = "Parameter"
  ) %>%
  # make final adjustments
  mutate(
    # re-code distribution (will serve to color density plots)
    Distribution = sub( "\\..*", "", Distribution),
    # add labels for effects grouping (will serve to split facets of the plot)
    Group = case_when(
      Parameter %in% pars$intercept ~ "intercept",
      Parameter %in% pars$base ~ "base",
      Parameter %in% pars$time ~ "time"
    ),
    # back-transform posteriors to raw DRS-2 (or DRS-2 per year) scale
    `DRS-2` = case_when(
      Parameter == "b_Intercept" ~ `DRS-2` * scl$SD$drs + scl$M$drs, # intercept
      Parameter != "b_Intercept" ~ `DRS-2` * scl$SD$drs # slopes
    )
  )

# change character columns to ordered factors such that ggplot plots them in a correct order
samples$Distribution <- factor( samples$Distribution , levels = rev( unique(samples$Distribution) ), ordered = T )
samples$Parameter <- factor( samples$Parameter , levels = rev( unlist(pars, use.names = F) ) , ordered = T )
samples$Group <- factor( samples$Group , levels = rev( unique(samples$Group) ) , ordered = T )

# Prepare part Fig. S5 (m1_nocov priors vs posteriors)
f.s5 <- list()

# loop through all three parameter groups
for ( i in rev( levels(samples$Group) ) ) {
  f.s5[[i]] <- samples[ samples$Group == i , ] %>%
    mutate(
      title = case_when( i == "intercept" ~ "Global intercept",
                         i == "base" ~ "Baseline correlates",
                         i == "time" ~ "Time-dependent effects")
    ) %>%
    ggplot( aes(y = Parameter, x = `DRS-2`, color = Distribution ) ) +
    stat_slab( geom = "slab", linetype = "solid", size = 2, fill = NA,  ) +
    coord_cartesian(
      ylim = case_when( i == "intercept" ~ c(1.3,1.5),
                        i == "base" ~ c(1.3,7.5),
                        i == "time" ~ c(1.3,8.5) )
    ) +
    geom_vline( xintercept = 0, linetype = "dashed" ) +
    scale_color_manual( values = cbPal[c(6,1)] ) +
    scale_x_continuous(
      limits = case_when( i == "intercept" ~ c(135, 143) ),
      name = case_when( i == "intercept" ~ "" , i == "base" ~ "DRS-2", i == "time" ~ "DRS-2 (points per year)")
    ) +
    scale_y_discrete( name = NULL, labels = rev( var_nms[ pars[[i]] , ] ) ) +
    facet_wrap( . ~ title) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 20),
      axis.text.y = element_text( vjust = -1.5 ),
      axis.ticks.y = element_blank()
    )
}

# put the facets of Fig S5 together
( ( f.s5$intercept / f.s5$base + plot_layout( heights = c(1,7) ) ) | f.s5$time )  +
  plot_layout( guides = "collect" ) & theme( legend.position = "bottom" )

# save as Fig S5
ggsave( "figures/Fig S5 posteriors vs priors.jpg", height = 1.25 * 8.53, width = 1.5 * 9.05, dpi = 300 )


# ----------- Tab S2 summary of posteriors -----------

# prepare a table for m1_nocov (the primary model) posteriors summary
t.s2 <- data.frame( par = as.character( unlist(pars) ), b = NA, PPI = NA, pd = NA ) %>% column_to_rownames("par")

# fill-in with appropriate estimates from draws$m1_nocov
for ( i in rownames(t.s2) ) t.s2[i,] <- c(
  sprintf( "%.2f", round( median(draws$m1_nocov[[i]] * scl$SD$drs + if(i == "b_Intercept") scl$M$drs else 0 ), 2 ) ),
  paste0( "[", sprintf( "%.2f", round( hdi(draws$m1_nocov[[i]] * scl$SD$drs + if(i == "b_Intercept") scl$M$drs else 0 ), 2 )[,1] ),
          ", ", sprintf( "%.2f", round( hdi(draws$m1_nocov[[i]] * scl$SD$drs + if(i == "b_Intercept") scl$M$drs else 0 ), 2 )[,2] ),
          "]" ),
  sprintf( "%.3f", round( sum( draws$m1_nocov[[i]] < 0 ) / nrow( draws$m1_nocov) , 3 ) )
)

# change some of the var_nms to a more appropriate format for Tab S2
var_nms[ rownames(t.s2)[2:8], ] <- paste0( var_nms[sub( "b_", "", rownames(t.s2)[2:8] ),], " (", var_nms[ rownames(t.s2)[2:8], ], ")" )

# make last adjustments to the table
t.s2 <- t.s2 %>%
  rename( "95% PPI" = "PPI", "Pr(b < 0)" = "pd" ) %>%
  add_column( Parameter = paste0( "\t", var_nms[rownames(t.s2), ] ), .before = "b") %>%
  add_row( Parameter = "Global intercept", .before = 1 ) %>%
  add_row( Parameter = "Baseline correlates", .before = 3 ) %>%
  add_row( Parameter = "Time-dependent effects", .before = 11 )

# save as csv
write.table( t.s2, file = "tables/Tab S2 summary of parameters posteriors.csv", sep = ";", row.names = F, na = "", quote = F )

# remove objects used to summarize posteriors
rm( list = c( "draws", "f.s5", "f.s6", "f2", "m", "pars", "post", "samples", "t.s2" ) )
gc()

# ----------- prepare posterior predictions -----------

# prepare posterior predictions for Processing Speed and Episodic Memory for Fig 5 in a multistep procedure
# using only the primary model to be reported in the main text
m <- readRDS( "models/m1_nocov.rds" )

# write down a list of predictors effect of which we're going to visualize in Fig. 4
prds <- c( "proc_spd" , "epis_mem" )

# prepare a raw data set for visualizations (for Fig 4 and Fig S3)
d4 <- d1 %>%
  select( id, time_y , drs_tot ) %>%
  filter( complete.cases(drs_tot) ) %>%
  rename( "drs" = "drs_tot" ) %>%
  # add median of all pre-surgery cognitive variables across imputations
  mutate(
    proc_spd = sapply( 1:imp , function(i) d3[[i]]$proc_spd ) %>% apply( . , 1 , median ),
    epis_mem = sapply( 1:imp , function(i) d3[[i]]$epis_mem ) %>% apply( . , 1 , median ),
    verb_wm = sapply( 1:imp , function(i) d3[[i]]$verb_wm ) %>% apply( . , 1 , median ),
    visp_mem = sapply( 1:imp , function(i) d3[[i]]$visp_mem ) %>% apply( . , 1 , median ),
    set_shift = sapply( 1:imp , function(i) d3[[i]]$set_shift ) %>% apply( . , 1 , median ),
    anxiety = sapply( 1:imp , function(i) d3[[i]]$anxiety ) %>% apply( . , 1 , median ),
    visp_wm = sapply( 1:imp , function(i) d3[[i]]$visp_wm ) %>% apply( . , 1 , median )
  )

# get quantile groups for each patients according to Processing Speed and Episodic Memory (for Fig 4)
# as well as all other pre-surgery cognitive variables (for Tab 2)
d5 <- d4[ d4$time_y < 0 , ] %>%
  mutate(
    # re-coding such that 1 means the lowermost quantile and 4 means the uppermost quantile
    proc_spd_pent = -ntile(proc_spd, 5) + 6,
    epis_mem_pent = -ntile(epis_mem, 5) + 6,
    verb_wm_pent = -ntile(verb_wm, 5) + 6,
    visp_mem_pent = -ntile(visp_mem, 5) + 6,
    set_shift_pent = -ntile(set_shift, 5) + 6,
    anxiety_pent = -ntile(anxiety, 5) + 6,
    visp_wm_pent = -ntile(visp_wm, 5) + 6
  )

# add these groups to the longitudinal data set (d4)
d4 <- d4 %>% left_join( d5 %>% select( id, contains("_pent") ) , by = "id" )


# ----------- Fig 3 posterior predictions -----------

# prepare a prediction dummy data set
d_seq <- list(
  epred_fix = list(), # prediction of the expectation (epred) based on fixed-effects only
  epred_all = list() # prediction of the expectation based on both fixed- and random-effects
)

# write down how many predictions (time-points) per predictor group/quantile to calculate
n_seq = 30

# loop through predictors
for ( i in names(d_seq) ) {
  for ( j in prds ) d_seq[[i]][[j]] <- expand.grid(
    as.vector( by( d4[[j]], d4[[ paste0(j, "_pent") ]] , median ) ),
    seq( from = -2, to = 12, length.out = n_seq )
  ) %>% `colnames<-` ( c(j, "time_y") ) %>%
      # add all variables needed
      mutate(
        time = time_y + scl$Md$time,
        proc_spd = if( j == "proc_spd") proc_spd else 0,
        epis_mem = if( j == "epis_mem") epis_mem else 0,
        grp = rep( 1:5, n_seq ), # in the expand.grid above, need to write predictor first, time second, otherwise this is incorrect
        verb_wm = 0, visp_mem = 0, set_shift = 0, anxiety = 0, visp_wm = 0, id = "sub000"
      )
}

# add predictions of the expectation (e_pred) to d_seq
for ( i in names(d_seq) ) {
  for ( j in prds ) d_seq[[i]][[j]] <- d_seq[[i]][[j]] %>%
      add_epred_draws( m , allow_new_levels = T, seed = 1, if (i == "epred_fix") { re_formula = NA } ) %>%
      mutate(.epred = scl$M$drs + scl$SD$drs * .epred) %>%
      median_hdi(.width = .95)
}

# prepare variable names for the plot
nms <- list( proc_spd = "Processing speed" , epis_mem = "Episodic memory" )

# prepare Fig 3
f3 <- list()

# fill-in plots for Processing speed and Episodic memory
for( i in prds ) {
  f3[[i]] <- d4[ , c( "time_y","id","drs") ] %>%
    mutate( grp = d4[ , paste0(i, "_pent") ] ) %>%
    ggplot( aes(x = time_y, y = drs, group = id) ) +
      geom_hline( yintercept = 139, linetype = "dotted", size = 1.5, alpha = 1 ) +
      geom_line( size = .5, alpha = .33 ) +
      geom_point( size = 4, alpha = .33 ) +
      geom_ribbon(data = d_seq$epred_all[[i]], alpha = .1,
                  aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, fill = grp )
                  ) +
      geom_ribbon(data = d_seq$epred_fix[[i]], alpha = .3,
                  aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, fill = grp )
      ) +
      geom_line(data = d_seq$epred_fix[[i]], size = 2.5,
                aes( x = time_y, y = .epred, color = grp )
      ) +
      scale_color_gradient( low = "red", high = "grey66" ) +
      scale_y_continuous(name = "DRS-2", limits = c(80,153), breaks = seq(80,140,20), labels = seq(80,140,20) ) +
      scale_x_continuous(name = "Time from surgery (years)", limits = c(-3,12), breaks = seq(-2,12,2), labels = seq(-2,12,2) ) +
      facet_wrap( ~ grp, nrow = 1,
                  labeller = as_labeller( c(`1` = "first pentile",
                                            `2` = "second pentile",
                                            `3` = "third pentile",
                                            `4` = "fourth pentile",
                                            `5` = "fifth pentile")
                  )) +
      ggtitle( nms[[i]] ) +
      theme( legend.position = "none", plot.title = element_text(hjust = .5) )
}

# arrange Fig 3 subplots for saving
( f3$proc_spd / f3$epis_mem ) + plot_annotation( tag_levels = "a" ) & theme( plot.tag = element_text(face = "bold") )

# save as Fig 3
ggsave( "figures/Fig 3 posterior predictions.jpg", width = 1.75 * 10.1, height = 1.25 * 11.2, dpi = 300 )


# ----------- Tab 3 posterior predictions -----------

# prediction dummy data set for a new patient's expectation prediction stratified by different levels of
# each of pre-surgery cognitive profile components
d_seq <- list()

# fill-in with values to predict for each of the pre-surgery cognitive variables
for ( i in doms ) d_seq[[i]] <- expand.grid(
  as.vector( by( d4[[i]], d4[[ paste0(i, "_pent") ]] , median ) ), c(0,1)
) %>% `colnames<-` ( c(i, "time_y") ) %>%
  # add all variables needed
  mutate(
    time = time_y + scl$Md$time,
    id = "sub000", # new (unobserved) patient
    proc_spd = if( i == "proc_spd") proc_spd else 0,
    epis_mem = if( i == "epis_mem") epis_mem else 0,
    verb_wm = if( i == "verb_wm") verb_wm else 0,
    visp_mem = if( i == "visp_mem") visp_mem else 0,
    set_shift = if( i == "set_shift") set_shift else 0,
    anxiety = if( i == "anxiety") anxiety else 0,
    visp_wm = if( i == "visp_wm") visp_wm else 0
  )

# compute prediction of the expectation for a new patient (i.e., based on both fixed- and random-effects)
for ( i in doms ) d_seq[[i]] <- d_seq[[i]] %>%
  posterior_epred( object = m , allow_new_levels = T, seed = 1 , re_formula = NULL ) %>%
  # calculate contrasts ( (post - pre) / 1 year ) for each pentile
  as.data.frame %>% mutate( `1st pentile` = ( V6 - V1 ) * scl$SD$drs,
                            `2nd pentile` = ( V7 - V2 ) * scl$SD$drs,
                            `3rd pentile` = ( V8 - V3 ) * scl$SD$drs,
                            `4th pentile` = ( V9 - V4 ) * scl$SD$drs,
                            `5th pentile` = ( V10 - V5 ) * scl$SD$drs) %>%
  # keep only the contrast for next calculations
  select( contains("pentile") ) %>%
  # calculate median contrast and 95% PPI for each pentile
  median_hdi( .width = .95 ) %>% select( !starts_with(".") ) %>%
  matrix( nrow = 5 , byrow = T ) %>% as.data.frame  %>%
  # extract in the form "Median [95% PPI]" for each pentile
  mutate( `Md [95% PPI]` = paste0( sprintf("%.2f", round( as.numeric(V1), 2) ), " [",
                                   sprintf("%.2f", round( as.numeric(V2), 2) ), ", ",
                                   sprintf("%.2f", round(as.numeric(V3), 2) ), "]" )
  ) %>% select( `Md [95% PPI]` )
  
# put it all together into a single table
t3 <- do.call( cbind.data.frame , d_seq ) %>%
  `colnames<-` ( var_nms[ doms , ] ) %>%
  `rownames<-` ( paste0( c("1st","2nd","3rd","4th","5th") , " pentile") ) %>%
  rownames_to_column( var = "stratum" )

# save the table as csv
write.table( t3, file = "tables/Tab 3 posterior predictions.csv", sep = ";", row.names = F, na = "", quote = F )


# ----------- Fig S3 per patient posterior predictions -----------

# simulate values for each patient each half year from 2 years before to 12 years after surgery
n_seq = 24

# re-format id in d4 to a factor
d4$id <- factor( d4$id, levels = unique(d4$id) )

# prepare data sets for prediction for each subject in the data set
d_seq <- expand.grid( seq( from = -2, to = 12, length.out = n_seq ), levels(d4$id) ) %>%
  `colnames<-` ( c("time_y","id") ) %>%
  mutate(time = time_y + scl$Md$time,
         proc_spd = NA, epis_mem = NA, verb_wm = NA, visp_mem = NA, set_shift = NA, anxiety = NA, visp_wm = NA
         )

# add subjects' median (w.r.t. imputations) cognitive profile
for ( i in levels(d4$id) ) d_seq[ d_seq$id == i, 4:10 ] <- d4[ d4$time_y < 0 & d4$id == i, colnames(d4) %in% doms ]

# get chunks of subjects to compute predictions for
# twenty-one chunks by six subjects so that there`s enough memory for each chunk
chunks <- split( unique(d4$id) , ceiling( seq_along(unique(d4$id)) / 6 ) )

# prepare a list for predictions stratified by subjects chunks
# calculating prediction of single data points based on fixed-effects, random-effects and remaining distributional (residual) parameters
preds <- list()

# add predictions to preds
for ( i in 1:length(chunks) ) preds[[i]] <- d_seq[ d_seq$id %in% chunks[[i]] , ] %>%
  add_predicted_draws( m, seed = 1) %>%
  mutate(drs = scl$M$drs + scl$SD$drs * .prediction) %>%
  median_hdi(.width = .95) %>%
  mutate(drs.upper = ifelse(drs.upper > 144, 144, drs.upper) ) # manual censoring

# collapse predictions to a single table
preds <- do.call( rbind.data.frame , preds )

# re-code id in d4 and preds such that the they are anonymized in the figure
all( levels(d4$id) == levels(preds$id) ) # TRUE, continue
levels(d4$id) <- paste0( "S", sprintf( "%.3d" , 1:126) )
levels(preds$id) <- paste0( "S", sprintf( "%.3d" , 1:126) )

# plot predictions and observed data for each subject separately
d4 %>%
  ggplot( aes(x = time_y, y = drs) ) +
    # add prediction lines and 95% compatibility intervals
    geom_line( data = preds, size = 2, aes(x = time_y, y = drs, group = id, color = cbPal[7] ) ) +
    geom_ribbon( data = preds, alpha = .1, aes(x = time_y, y = drs,
                                               ymin = drs.lower, ymax = drs.upper,
                                               group = id, fill = cbPal[7] ) ) +
    # add observed points
    geom_point( color = "black", size = 4 ) +
    # finish the plot with the last cosmetic changes
    facet_wrap( ~id, nrow =  14 ) + # arrange to a 14 x 9 grid
    labs( x = "Time after surgery (Years)", y = "DRS-2 (0-144 points)") +
    theme( legend.position = "none" )

# save as Fig S3
ggsave( "figures/Fig S3 per-patient posterior predictions.jpg", width = 2 * 10.1, height = 6 * 5.39, dpi = 300 )


# ----------- summarize PSIS-LOO -----------

# read the file
l <- readRDS("models/dbs_longCOG_psis-loo.rds")

# compute maximal Pareto-Ks
t.loo <- data.frame( `pareto-k > 0.5` = rep(NA,imp), `pareto-k > 0.7` = NA, `pareto-k > 1.0` = NA, `max pareto-k` = NA )

# fill-in the table
for ( i in 1:imp ) t.loo[i,] <- c(
  sum( l[[i]]$diagnostics$pareto_k > 0.5 ),
  sum( l[[i]]$diagnostics$pareto_k > 0.7 ),
  sum( l[[i]]$diagnostics$pareto_k > 1.0 ),
  max( l[[i]]$diagnostics$pareto_k)
)

# sum
table( t.loo$pareto.k...0.5 ) # one case with Pareto-k > 0.5
table( t.loo$pareto.k...0.7 )
table( t.loo$pareto.k...1.0 )

# look at the value of models with Pareto-k > 0.5
t.loo[ t.loo$pareto.k...0.5 > .5 , ] # model no. 54, max(pareto-k) == 0.53


# ----------- linear vs non-linear fit -----------

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .99 # adapt_delta parameter
s = 87543 # seed for reproducibility

# set-up the linear model (doing PPC for the primary model only, ie, m1_nocov,  the rest should be "centered"
# around it as they are the same model with added parameters with zero-centered priors)
f4.drs <- list( linear = bf( drs | cens(cens_drs) ~ time + (1 + time || id) ),
                spline = bf( drs | cens(cens_drs) ~ t2(time) + (1 + time || id) )
)

# set-up priors (using brms default non-informative to allow for information from data to prevail)
p4 <- NULL

# fit the models
m4 <- list()
for ( i in names(f4.drs) ) m4[[i]] <- brm(
  formula = f4.drs[[i]], family = student(), prior = p4, data = d3[[69]], # using only drs and time, so imputed datasets in d3 are equivalent to scaled d1
  seed = s, iter = it, warmup = wu, chains = ch, control = list( adapt_delta = ad ),
  file = paste0( getwd(), "/models/m4_time_only_", i, ".rds"),
  save_model = paste0( getwd(), "/models/m4_time_only_", i, ".stan")
)

# find the largest Rhat to get idea about convergence
sapply( names(m4), function(i) max( rhat(m4[[i]]) ) )

# add PSIS-LOO to both model for influential variables check and model comparisons
for ( i in names(m4) ) m4[[i]] <- add_criterion( m4[[i]] , criterion = c("loo","waic") )

# plot pareto-ks
par( mfrow = c(1,2) ) # plot them next to each other
for ( i in names(m4) ) plot( loo(m4[[i]]) )
par( mfrow = c(1,1) ) # return graphic options to default

# compare the models via PSIS-LOO
loo( m4$linear, m4$spline, moment_match = T )


# ----------- Fig 4 linear vs non-linear fit -----------

# prepare data to be predicted
d_seq <- data.frame( time_y = seq(-2,12,length.out = 50), id = NA ) %>%
  mutate( time = time_y + scl$Md$time )

# add predictions of expectation (epreds, based on fixed-effects only) from both linear and non-linear models
preds <- lapply( names(m4) , function(i)
  d_seq %>%
    add_epred_draws( m4[[i]], re_formula = NA ) %>%
    mutate( .epred = .epred * scl$SD$drs + scl$M$drs ) %>%
    median_hdi( .width = .95 ) %>%
    add_column( Model = factor( ifelse(i == "linear", "Linear", "Smooths"), levels = c("Smooths","Linear"), ordered = T ) )
  )

# collapse the linear and spline predictions for plotting purposes to a single file
preds <- do.call( rbind.data.frame, preds )

# plot linear and spline fits over each other
preds %>%
  ggplot( aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, color = Model, fill = Model) ) +
    geom_ribbon( alpha = .2 , color = NA ) +
    geom_line( size = 2.5 , alpha = .75 ) +
    scale_y_continuous(name = "DRS-2", limits = c(120,145), breaks = seq(120,150,10), labels = seq(120,150,10) ) +
    scale_x_continuous(name = "Time from surgery (years)", limits = c(-2,12), breaks = seq(-2,12,2), labels = seq(-2,12,2) ) +
    scale_color_manual( values = cbPal[c(3,2)] ) +
    scale_fill_manual( values = cbPal[c(3,2)] ) +
    theme( legend.position = c(0.11,0.13) )

# save as Fig 4
ggsave( "figures/Fig 4 linear vs non-linear fit.jpg", width = 1.25 * 10.1, height = 1.25 * 5.39, dpi = 300 )
