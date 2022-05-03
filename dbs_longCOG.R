# All analyses reported in the article (Mana et al., in review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 19043).

# I used the following versions of packages employed: dplyr_1.0.7, tidyverse_1.3.0, DiagrammeR_1.0.9, DiagrammeRsvg_0.1,
# rsvg_2.3.1, missMDA_1.18, psych_2.2.3, brms_2.16.3, psych_2.2.3, tidybayes_2.3.1, ggplot2_3.3.3 and patchwork_1.1.1

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "tidyverse", # for pivot_longer
  "DiagrammeR", # for flowcharts
  "DiagrammeRsvg", # for saving plots crated in DiagrammeR
  "rsvg", # for saving plots crated in DiagrammeR
  "missMDA", # for imputation
  "psych", # for EFA
  "brms", # for Bayesian model fitting / interface with Stan
  "loo", # for PSIS-LOO based operations
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

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
s = 87542 # seed for reproducibility
ch = 4 # number of chains
it = 2500 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .99 # adapt_delta parameter

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present,
# prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler 
# and version.

# read the dataset and prepare subsets for individual analyses
# In this file, my strategy is to label the sequential version
# of dataset as d# (where # goes for a number of the data version)
d0 <- read.csv( "data/20220426_dbs_longCOG_data.csv" , sep = "," )
d1 <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data
d2 <- d1[ d1$ass_type == "pre" , ] # only pre-surgery assessments of included patients

# ----------- participant inclusion flowchart  -----------

# check that when selecting only rows containing pre-surgery assessment,
# there ain't no patient duplicated
isTRUE(
  all.equal(
    d0[ d0$ass_type == "pre" , ]$id,
    unique( d0[ d0$ass_type == "pre" , ]$id )
  )
) # TRUE

# create a table summarizing reasons for excluding patients
t <- table(d0[ d0$ass_type == "pre" , ]$why_excluded)
print(t)

# using numbers from t create an inclusion/exclusion flowchart
f1 <- " digraph {
  
  /// define nodes which will include numbers reflecting the inclusion/exclusion process
  node [ fontname = Calibri, fontsize = 24, shape = box, style = rounded , width = 2.5 , margin= 0.25 , penwidth = 2 ];
  
  /// box for all patients, i.e., sum(t) = 200 patients
  all_pats [ label =<
  <b>200 consecutive PD patients </b><br/><br/>Local database 2000-2020<br/>General University Hospital<br/>in Prague
  >];
         
  /// box for all STN-DBS patients, i.e., sum(t[c(1,4,5,6,7,8,9,11)])) = 173 patients
  all_stn [ label =<
  <b>173 patients </b><br/><br/>implanted with STN-DBS
  >];
  
  /// box for all included patients, i.e., t[4] = 126 patients
  all_incl [ label =<
  <b>126 patients </b><br/><br/>followed-up longitudinally
  >];
  
  /// nodes for exluded patients, specifying only these characteristics that will differ from the nodes above
  node [ fixedsize = T, width = 5.5, height = 2.2 ];
  
  /// box for non-STN-DBS patients, i.e., t[c(3,13,2,10,12)] with sum(t[c(3,13,2,10,12)]) = 27 patients
  excl_nostn [ label =<
  <b>27 patients excluded due to</b><br align = 'left'/><br align = 'left'/>
  12 GPi-DBS<br align = 'left'/>
  4 VIM-DBS<br align = 'left'/>
  4 duodopa<br align = 'left'/>
  5 rejected<br align = 'left'/>
  2 suspended<br align = 'left'/>
  >];
  
  /// box for STN-DBS excluded patients, i.e., t[c(1,7,9, 8, 6,11, 5)] with sum(t[c(1,7,9,8,6,11,5)]) = 57 patients
  excl_stn [ label =<
  <b>47 patients excluded due to</b><br align = 'left'/><br align = 'left'/>
  24 pre-surgery data missing<br align = 'left'/>
  18 follow-up data missing<br align = 'left'/>
  3 unilateral STN-DBS<br align = 'left'/>
  2 not speaking Czech<br align = 'left'/>
  >];
  
  /// dummy nodes for horizontally forking out of the vertical 'inclusion flow' to excluded sides
  node [ shape = rectangle, width = 0, height = 0, label = '', fill = black ];
  
  /// create directed edges in the inclusion (from dummy/fork vertically to inclusion boxes)
  /// and the exlusion (from forks horizontally to exclusion boxes) parts of the flowchart
  /// first make the arrows bigger
  edge [ arrowsize = 2.5, penwidth = 2.5 ]
  
  /// next specifiy paths
  fork1 -> all_stn; fork2 -> all_incl;
  
  /// for the horizontal paths use 'rank = same' to ensure they nodes are level
  { rank = same ; fork1 -> excl_nostn }
  { rank = same ; fork2 -> excl_stn }
  
  /// create non-directed edges from the inclusion boxes to dummy/fork boxes (vertically)
  edge [ dir = none, penwidth = 2.5 ]
  all_pats -> fork1; all_stn -> fork2;

  /// finally seperate dummy/fork nodes from exclusion boxes by some reasonable distance (in inches)
  nodesep = 1
  }" 

# save the flowchart as png
grViz(f1) %>% export_svg %>% charToRaw %>% rsvg_png("figures/Fig.1 inclusion-exclusion flowchart.png")


# ----------- sample description  -----------

# list all stimulation parameters
pars <- names(d1)[which(names(d1)=="current_r_mA"):which(names(d1)=="frequency_l_Hz")]

# prepare dataframe to fill-in
t1 <- data.frame(
  Md = rep( NA , length(pars) ), `Min-Max` = NA, M = NA, SD = NA, row.names = pars
)

# fill-in with statistics of all stimulation parameters
for ( i in pars ) {
  # median, rounded to hundredths
  t1[ i , "Md" ] <- sprintf(
    "%.2f" , round( median( d1[[i]], na.rm = T ) , digits = 2 )
  )
  # min-max, rounded to decimals
  t1[ i , "Min.Max" ] <- paste0(
    sprintf( "%.1f" , round( min( d1[[i]], na.rm = T ) , digits = 1 ) ), "-",
    sprintf( "%.1f" , round( max( d1[[i]], na.rm = T ) , digits = 1 ) )
  )
  # mean, rounded to hundredths
  t1[ i , "M" ] <- sprintf(
    "%.2f" , round( mean( d1[[i]], na.rm = T ) , digits = 2 )
  )
  # SD, rounded to hundredths
  t1[ i , "SD" ] <- sprintf(
    "%.2f" , round( sd( d1[[i]], na.rm = T ) , digits = 2 )
  )
}

# list all variables that will be included in Tab. 2 (baseline characteristics)
vars <- names(d2)[which(names(d2)=="age_stim_y"):which(names(d2)=="fp_dr")]

# prepare a dataframe to fill-in
t2 <- data.frame(
  N = rep( NA, length(vars) ), Md = NA, `Min-Max` = NA, M = NA, SD = NA, row.names = vars
)

# fill-in with statistics of all variables but sex (which is nominal)
for ( i in vars[-3] ) {
  # number of data points
  t2[ i , "N" ] <- sum( !is.na(d2[[i]]) )
  # median, rounded to integers
  t2[ i , "Md" ] <- sprintf(
    "%.0f" , round( median(d2[[i]], na.rm = T ) , digits = 0 )
  )
  # min-max, rounded to integers
  t2[ i , "Min.Max" ] <- paste0(
    sprintf( "%.0f" , round( min(d2[[i]], na.rm = T ) , digits = 0 ) ), "-",
    sprintf( "%.0f" , round( max(d2[[i]], na.rm = T ) , digits = 0 ) )
  )
  # mean, rounded to hundredths
  t2[ i , "M" ] <- sprintf(
    "%.2f" , round( mean(d2[[i]], na.rm = T ) , digits = 2 )
  )
  # SD, rounded to hundredths
  t2[ i , "SD" ] <- sprintf(
    "%.2f" , round( sd(d2[[i]], na.rm = T ) , digits = 2 )
  )
}

# add frequency of males to the table sex row
t2[ "sex" , c("N","Min.Max")] <- c(
  # number of entries
  sum( !is.na(d2$sex) ),
  # add frequency (percentage) to Min.Max
  paste0(
    # frequency
    table(d2$sex)["male"], " (",
    # percentage, rounded to integers
    sprintf(
      "%.0f" , round(
        100 * ( table(d2$sex)[ "male" ] / sum( table(d2$sex) ) ), 0
      )
    ), " %)"
  )
)

# prepare a histogram of the distribution of assessments across time (Fig. 2A)
f2 <- list()

# need to use the complete.cases command because for three patient I have duplicated
# rows due to more than one stimulation parameter
f2$hist <- d1[ complete.cases(d1$drs_tot) , ] %>% 
  ggplot( aes(x = time_y) ) +
  stat_bin( binwidth = .5, position = position_nudge(x = -.5*.5) ) + # creates bars
  stat_bin(
    binwidth = .5, geom = "text", aes(label = ..count..), vjust = -1.0,
    position = position_nudge(x = -.5*.5), size = 6
  ) + # add numbers
  labs(
    x = "Time from STN-DBS surgery (years)",
    y = "Number of Assessments"
  ) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0, 69),
    breaks = seq(0, 60, 10), labels = seq(0, 60, 10)
  ) +
  scale_x_continuous(
    limits = c(-2, 12), breaks = seq(-2, 12, 1), labels = seq(-2, 12, 1)
  )

# prepare a bin plot showing distribution of the number of assessments per patient (Fig.2B)
f2$bin <- table( d1[ complete.cases(d1$drs_tot) , ]$id ) %>%
  as.data.frame() %>%
  ggplot( aes(x = Freq) ) +
  geom_bar( width = .25 ) +
  geom_text(
    stat = "count",
    aes(label = ..count..),
    vjust = -1.0, size = 6
  ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0, 65) ) +
  labs(
    x = "Number of Assessments per Patient",
    y = "Number of Patients"
  )

# arrange Fig.2A and Fig.2B for printing
f2$hist / f2$bin + plot_annotation( tag_levels = "A" )

# save as Fig.2
ggsave( "figures/Fig.2 distribution of assessments.png" , height = 2.5 * 6.12 , width = 1.5 * 11.6 , dpi = "retina" )


# ----------- pre-surgery cognitive profile  -----------

# for EFA keep only id and cognitive tests in d2
d2 <- d2[ , c( 2, which(names(d2) == "tmt_a"):which(names(d2) == "fp_dr"), which(names(d2) %in% paste0("staix",1:2)) ) ]

# multiply impute missing values in d2
# first find out the optimal number of components for multiple combinations
nb <- estim_ncpPCA( d2[,-1] , ncp.max = 10 )

# impute via PCA-based multiple imputation (n = 100 imputations)
set.seed(s) # set seed for reproducility
d2.imp <- MIPCA( d2[ ,- 1] , ncp = nb$ncp )

# fit EFA to each imputed dataset with 3:8 factors
efa <- lapply(
  # loop through all 100 datasets
  1:length(d2.imp$res.MI), function(i)
    lapply(
      # loop through three to eight latent factors for each imputation
      3:8, function(j)
        fa(
          d2.imp$res.MI[[i]], # one-by-one use each imputed dataset
          nfactors = j, # fit 3-8 factor solutions
          rotate = "varimax", # rotate varimax to enforce orthogonality and for interpretation purposes
          scores = "regression" # compute regression scores for each patient
        )
    )
)

# prepare a table for performance indexes of each solution for each imputed dataset
fat <- array(
  # create n x m x l array
  dim = c(
    length( efa[[1]] ), # m = number of factor solutions I computed above
    6, # n = six dimensions for performance indexes
    length( d2.imp$res.MI ) # l = number of imputations
  ),
  # add dimension names
  dimnames = list(
    paste0( 1:length(efa[[1]])+2, "_factors" ), # number of factors
    c("TLI", "BIC", "RMSEA", "RMSEA_90_CI_low", "RMSEA_90_CI_upp", "var_account"), # performance indexes
    1:length(d2.imp$res.MI) # number of imputed dataset
  )
)

# loop through the array and gather all performance indexes
# loop through the array and gather all performance indexes
for ( i in 1:length(d2.imp$res.MI) ) {
  for ( j in 1:length(efa[[i]]) ) {
    fat[ j , , i ] <- c(
      round( efa[[i]][[j]]$TLI , 3 ), # Tucker-Lewis index, > 0.9 is considered good
      round( efa[[i]][[j]]$BIC , 3 ), # Bayesian Information Criterion, lower values are considered better
      round( efa[[i]][[j]]$RMSEA[1], 3 ), # root-mean square error of approximation, < 0.08 is considered good
      round( efa[[i]][[j]]$RMSEA[1], 3 ), # 90% CI RMSEA, lower boundary
      round( efa[[i]][[j]]$RMSEA[1], 3 ), # 90% CI RMSEA, upper boundary
      round(efa[[i]][[j]]$Vaccounted["Cumulative Var",j+2], 3) # total variance "accounted for" by all included factors
    )
  }
}

# print the upper 90% CI of RMSEA for all model/imputation pairs
# get rid of values > .08 to see only models with good fit
t(fat[ , "RMSEA_90_CI_upp" , ] ) %>%
  as.data.frame %>%
  mutate( across( everything() , ~ replace( . , . > .08 , "") ) )

# go through all six-factor solutions' loading matrices
# create a convenience function so that I don't go crazy immediately
print_load <- function(i) print( efa[[i]][[4]]$loadings, cutoff = .4, sort = T )
# weird ones - # 24, 31, 38, 45, 52, 63, 82
# include reverse scored factors - # 15, 19, 24, 38, 45, 46, 47, 63, 91, 98

# save the imputed data set for future control
saveRDS( d2.imp , "data/20220503_100_imputed_df.rds" )




# selecting the 6-factor solution due to the best BIC, reasonable stats w.r.t. competition and interpretable factors
# which are missing in models with better stats that include a one-item factor
# time to name cognitive domains inferred from pre-surgery tests
for ( i in c("loadings","scores","Vaccounted")) {
  colnames( efa[[4]][[i]] ) <- c(
    "ment_flex", # 'mental flexibility' is loaded on primarily by Stroop task and verbal fluencies
    "epis_mem", # 'episodic memory' is loaded on primarily by RAVLT
    "visp_wm", # 'visuo-spatial working memory/switching' is loaded on primarily by Spatial Span, TMTs and LNS
    "visp_mem", # 'visuo-spatial memory' is loaded on primarily by Family Picture Test
    "verb_wm", # 'verbal working memory' is loaded on primarily by Digit Spans and Letter-Number sequencing (LNS)
    "anxiety" # 'anxiety! is loaded on primarily by STAI
  )
}

# prepare Tab.3 with loadings and variance accounted for by the 6-factor model
tab3 <- unclass( efa[[4]]$loadings ) %>%
  as.data.frame() %>%
  # add variance accounted for
  bind_rows( as.data.frame(efa[[4]]$Vaccounted[ c("Proportion Var", "Cumulative Var"), ]) ) %>%
  # rename and round
  mutate(
    ment_flex = sprintf( "%.3f" , round( ment_flex , 3 ) ),
    epis_mem = sprintf( "%.3f" , round( epis_mem , 3 ) ),
    visp_wm = sprintf( "%.3f" , round( visp_wm , 3 ) ),
    visp_mem = sprintf( "%.3f" , round( visp_mem , 3 ) ),
    verb_wm = sprintf( "%.3f" , round( verb_wm , 3 ) ),
    anxiety = sprintf( "%.3f" , round( anxiety , 3 ) )
  ) %>%
  # give rownames their own variable/column
  rownames_to_column(
    var = "variable"
  )

# save all EFA models
saveRDS( object = efa, file = "models/efa.rds" )

# merge longitudinal d1 with baseline factor scores from efa[[4]] (joining by id)
d1 <- d1 %>%
  left_join( cbind.data.frame(id = d2$id, efa[[4]]$scores ) , by = "id" ) %>%
  filter( complete.cases(drs_tot) ) # get rid of three dummy rows due to more than one stimulation parameter (no DRS-2)

# ----------- pre-processing for longitudinal analyses  -----------

# save scaling values for variables to be included in longitudinal analyses
scl <- list(
  M = list(
    # longitudinal variables, i.e., DRS-2 (outcome), BDI and LEDD (confounders)
    drs = mean( d1$drs_tot , na.rm = T ), # 136.90
    bdi = mean( d1$bdi , na.rm = T ), # 10.95
    led = mean( d1$ledd_mg , na.rm = T ), # 1194.03
    # pre-surgery prognostic factors/cognitive profile (calculate M from pre-surgery assessments only), all zero
    ment_flex = mean( d1[d1$ass_type=="pre", ]$ment_flex , na.rm = T ),
    epis_mem = mean( d1[d1$ass_type=="pre", ]$epis_mem , na.rm = T ),
    visp_wm = mean( d1[d1$ass_type=="pre", ]$visp_wm , na.rm = T ),
    visp_mem = mean( d1[d1$ass_type=="pre", ]$visp_mem , na.rm = T ),
    verb_wm = mean( d1[d1$ass_type=="pre", ]$verb_wm , na.rm = T ),
    anxiety = mean( d1[d1$ass_type=="pre", ]$anxiety , na.rm = T )
  ),
  SD = list(
    # longitudinal variables, i.e., DRS-2 (outcome), BDI and LEDD (confounders)
    drs = sd( d1$drs_tot , na.rm = T ), # 7.72
    bdi = sd( d1$bdi , na.rm = T ), # 7.19
    led = sd( d1$ledd_mg , na.rm = T ), # 686.69
    # pre-surgery prognostic factors/cognitive profile (calculate SD from pre-surgery assessments only)
    ment_flex = sd( d1[d1$ass_type=="pre", ]$ment_flex , na.rm = T ), # 0.91
    epis_mem = sd( d1[d1$ass_type=="pre", ]$epis_mem , na.rm = T ), # 0.93
    visp_wm = sd( d1[d1$ass_type=="pre", ]$visp_wm , na.rm = T ), # 0.89
    visp_mem = sd( d1[d1$ass_type=="pre", ]$visp_mem , na.rm = T ), # 0.99
    verb_wm = sd( d1[d1$ass_type=="pre", ]$verb_wm , na.rm = T ), # 0.87
    anxiety = sd( d1[d1$ass_type=="pre", ]$anxiety , na.rm = T ) # 0.90
  ),
  # add median time of pre-surgery assessment for GLMMs intercepts
  Md = list(
    time = -median( d1[d1$ass_type == "pre", ]$time_y , na.rm = T ) # 0.30
  )
)

# prepare dataset for all longitudinal analyses
d3 <- d1 %>% mutate(
  # base variable for a longitudinal model
  time = time_y + scl$Md$time,
  drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
  bdi = ( bdi - scl$M$bdi ) / scl$SD$bdi,
  led = ( ledd_mg - scl$M$led ) / scl$SD$led,
  sex = as.factor( sex ), # for better estimation of BDI in the second and third models
  # censoring variables
  cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ), # right censoring for DRS == 144
  cens_bdi = ifelse( bdi == min(bdi, na.rm = T) , "left" , "none" ), # left censoring for BDI == 0
  # pre-surgery prognostic factors/cognitive profile, scale such that higher value imply deficit
  ment_flex = -( ment_flex - scl$M$ment_flex ) / scl$SD$ment_flex,
  epis_mem = -( epis_mem - scl$M$epis_mem ) / scl$SD$epis_mem,
  visp_wm = -( visp_wm - scl$M$visp_wm ) / scl$SD$visp_wm,
  visp_mem = -( visp_mem - scl$M$visp_mem ) / scl$SD$visp_mem,
  verb_wm = -( verb_wm - scl$M$verb_wm ) / scl$SD$verb_wm,
  anxiety = -( anxiety - scl$M$anxiety ) / scl$SD$anxiety
) %>% select(
  # keep only variables of interest
  id, time, drs, cens_drs, bdi, cens_bdi, sex, led, # outcomes, time and sex
  ment_flex, epis_mem, visp_wm, visp_mem, verb_wm, anxiety # pre-surgery prognostic factors/cognitive profile
)


# ----------- description models  -----------

# prepare a list for all models and PSIS-LOOs
m <- list()
l <- list()

# model with no covariates
# set-up the linear model
f1.drs <- bf( drs | cens(cens_drs) ~ 1 + time + (1 + time | id) )

# set-up priors
p1 <- c(
  # fixed-effects
  prior( normal(0.3, .1), class = Intercept ),
  prior( normal(-.2, .1), class = b, coef = time ),
  # random-effects
  prior( normal(0, .1), class = sd, coef = Intercept, group = id ),
  prior( normal(0, .1), class = sd, coef = time, group = id ),
  # variances
  prior( exponential(1), class = sd ),
  prior( exponential(1), class = sigma ),
  # covariances
  prior( lkj(2), class = cor ),
  # "degrees of freedom"/"normality parameter" nu
  prior( gamma(2, 0.1), class = nu )
)

# fit the total effect of time GLMM
m$desc$no_cov <- brm(
  formula = f1.drs, family = student(), prior = p1, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/descriptive.rds"
)

# add LOO and WAIC criteria
add_criterion( m$desc$no_cov , criterion = c("loo","waic") )

# save the PSIS-LOO
l$desc$no_cov <- loo(m$desc$no_cov)

# save the priors for documentation purposes
pr <- list( desc = list( no_cov = prior_summary(m$desc$no_cov) ) )

# some posterior checks
max( rhat(m$desc$no_cov) ) # 1.005
max( loo(m$desc$no_cov)$diagnostics$pareto_k ) # 0.57

# model with covariates
# set-up the linear model
f2.drs <- bf( drs | cens(cens_drs) ~ 1 + mi(bdi) + mi(led)  + time+ (1 + time | id) )
f2.bdi <- bf( bdi | cens(cens_bdi) + mi() ~ 1 + mi(led) + sex + time + (1 + time | id) )
f2.led <- bf( led | mi() ~ t2(time) + (1 | id) )

# set contrast for sex
contrasts(d3$sex) <- -contr.sum(2)/2 # female = -0.5, male = 0.5

# set-up priors
p2 <- c(
  # DRS-2
  prior( normal(0.3, .1), class = Intercept , resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time , resp = drs  ),
  prior( normal(0, .1), class = b, coef = mibdi, resp = drs ),
  prior( normal(0, .1), class = b, coef = miled, resp = drs ),
  prior( normal(0, .1), class = sd, coef = Intercept, group = id , resp = drs  ),
  prior( normal(0, .1), class = sd, coef = time, group = id , resp = drs  ),
  prior( exponential(1), class = sd , resp = drs ),
  prior( exponential(1), class = sigma , resp = drs ),
  prior( gamma(2, 0.1), class = nu , resp = drs ),
  # BDI-II
  prior( normal(.6, .5), class = Intercept, resp = bdi ),
  prior( normal(0, .5), class = b, coef = time, resp = bdi ),
  prior( normal(0, .5), class = b, coef = sex1, resp = bdi ),
  prior( normal(0, .5), class = b, coef = miled, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  prior( exponential(1), class = sd, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  prior( gamma(2, 0.1), class = nu, resp = bdi ),
  # LEDD
  # keeping default for fixed-effect on LEDD
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = led ),
  prior( exponential(1), class = sd, resp = led ),
  prior( exponential(1), class = sigma, resp = led ),
  prior( gamma(2, 0.1), class = nu, resp = led ),
  # group-level correlation matrix
  prior( lkj(2), class = cor )
)

# fit the model
m$desc$w_cov <- brm(
  formula = f2.drs + f2.bdi + f2.led + set_rescor(F),
  family = student(), prior = p2, data = d3, sample_prior = T,
  seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/descriptive_w_cov.rds"
)

# add LOO and WAIC criteria
add_criterion( m$desc$w_cov , criterion = c("loo","waic") , resp = "drs" )

# save the PSIS-LOO
l$desc$w_cov <- loo(m$desc$w_cov, resp = "drs")

# save the priors for documentation purposes
pr$desc$w_cov <- prior_summary(m$desc$w_cov)

# some posterior checks
max( rhat(m$desc$w_cov) ) # 1.011
max( loo(m$desc$w_cov, resp = "drs")$diagnostics$pareto_k ) # 0.60


# ----------- prognostic models  -----------

# write down the cognitive domains
doms <- colnames( efa[[4]]$scores )

# model with no covariates
# set-up the linear model
f3.drs <- bf(
  as.formula(
    paste0(
      "drs | cens(cens_drs) ~ 1 + ",
      paste( "time", doms, sep = " * " , collapse = " + " ),
      " + (1 + time | id)"
    )
  )
)

# set-up priors
p3 <- c(
  # fixed effects
  prior( normal(0.3, .1), class = Intercept ),
  prior( normal(-.2, .1), class = b, coef = time ),
  prior( normal(0, .1), class = b, coef = ment_flex ),
  prior( normal(0, .1), class = b, coef = epis_mem ),
  prior( normal(0, .1), class = b, coef = visp_wm ),
  prior( normal(0, .1), class = b, coef = visp_mem ),
  prior( normal(0, .1), class = b, coef = verb_wm ),
  prior( normal(0, .1), class = b, coef = anxiety ),
  prior( normal(0, .1), class = b, coef = time:ment_flex ),
  prior( normal(0, .1), class = b, coef = time:epis_mem ),
  prior( normal(0, .1), class = b, coef = time:visp_wm ),
  prior( normal(0, .1), class = b, coef = time:visp_mem ),
  prior( normal(0, .1), class = b, coef = time:verb_wm ),
  prior( normal(0, .1), class = b, coef = time:anxiety ),
  # random effects
  prior( normal(0, .1), class = sd, coef = Intercept, group = id ),
  prior( normal(0, .1), class = sd, coef = time, group = id ),
  prior( exponential(1), class = sd ),
  # other distribution parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu ),
  # correlation matrices
  prior( lkj(2), class = cor )
)

# fit the GLMM
m$prog$no_cov <- brm(
  formula = f3.drs, family = student(), prior = p3, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/prognostic.rds"
)

# add LOO and WAIC criteria
add_criterion( m$prog$no_cov , criterion = c("loo","waic") )

# save the PSIS-LOO
l$prog$no_cov <- loo(m$prog$no_cov)

# save the priors for documentation purposes
pr$prog$no_cov <- prior_summary(m$prog$no_cov)

# some posterior checks
max( rhat(m$prog$no_cov) ) # 1.004
max( loo(m$prog$no_cov)$diagnostics$pareto_k ) # 0.55

# model with covariates
# set-up the linear model
f4.drs <- bf(
  as.formula(
    paste0(
      "drs | cens(cens_drs) ~ 1 + mi(bdi) + mi(led) + ",
      paste( "time", doms, sep = " * " , collapse = " + " ),
      " + (1 + time | id)"
    )
  )
)

# linear models for BDI-II and LEDD remain the same as in m2
f4.bdi <- f2.bdi
f4.led <- f2.led

# set-up priors
p4 <- c(
  # DRS-2
  prior( normal(0.3, .1), class = Intercept , resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time , resp = drs  ),
  prior( normal(0, .1), class = b, coef = mibdi, resp = drs ),
  prior( normal(0, .1), class = b, coef = miled, resp = drs ),
  prior( normal(0, .1), class = b, coef = ment_flex, resp = drs ),
  prior( normal(0, .1), class = b, coef = epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_mem, resp = drs  ),
  prior( normal(0, .1), class = b, coef = verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = anxiety, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:ment_flex, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_mem, resp = drs  ),
  prior( normal(0, .1), class = b, coef = time:verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:anxiety, resp = drs ),
  prior( normal(0, .1), class = sd, coef = Intercept, group = id , resp = drs  ),
  prior( normal(0, .1), class = sd, coef = time, group = id , resp = drs  ),
  prior( exponential(1), class = sd , resp = drs ),
  prior( exponential(1), class = sigma , resp = drs ),
  prior( gamma(2, 0.1), class = nu , resp = drs ),
  # BDI-II
  prior( normal(.6, .5), class = Intercept, resp = bdi ),
  prior( normal(0, .5), class = b, coef = time, resp = bdi ),
  prior( normal(0, .5), class = b, coef = sex1, resp = bdi ),
  prior( normal(0, .5), class = b, coef = miled, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  prior( exponential(1), class = sd, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  prior( gamma(2, 0.1), class = nu, resp = bdi ),
  # LEDD
  # keeping default for fixed-effect on LEDD
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = led ),
  prior( exponential(1), class = sd, resp = led ),
  prior( exponential(1), class = sigma, resp = led ),
  prior( gamma(2, 0.1), class = nu, resp = led ),
  # group-level correlation matrix
  prior( lkj(2), class = cor )
)

# fit the "prognostic" model
m$prog$w_cov <- brm(
  formula = f4.drs + f4.bdi + f4.led + set_rescor(F),
  family = student(), prior = p4, data = d3, sample_prior = T,
  seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/prognostic_w_cov.rds"
)

# add LOO and WAIC criteria
add_criterion( m$prog$w_cov , criterion = c("loo","waic") , resp = "drs" )

# save the PSIS-LOO
l$prog$w_cov <- loo(m$prog$w_cov, resp = "drs")

# save the priors for documentation purposes
pr$prog$w_cov <- prior_summary(m$prog$w_cov)

# some posterior checks
max( rhat(m$prog$w_cov) ) # 1.060
max( loo(m$prog$w_cov, resp = "drs")$diagnostics$pareto_k ) # 0.46


# ----------- joint quality check & model comparisons  -----------

# compare the models via PSIS-LOO
# for now compare m1 to m3 and m2 to m4 (because the have different number of points for PSIS-LOO)
# prepare pointwise elpd matrixes
lpd_point <- lapply(
  names(l$desc), function(i)
    cbind(
      l$desc[[i]]$pointwise[,"elpd_loo"],
      l$prog[[i]]$pointwise[,"elpd_loo"]
    )
)

# create the table
comp <- lapply(
  names(l$desc) , function(i)
    # compare decriptive vs. prognostic models within no covariates/with covariates
    # cannot compare models with vs within covariates due to missing raw dat in covariates models
    loo_compare(
      l$desc[[i]], l$prog[[i]]
    )[ paste0( c("m$desc$","m$prog$") , i ) , c("elpd_diff","se_diff") ] %>%
    as.data.frame() %>%
    # change rownames to something intelligible
    rownames_to_column( var = "model" ) %>%
    mutate(
      model = ifelse(
        grepl("desc", model), "descriptive", "prognostic"
      )
    )
)

# give them names
names(lpd_point) <- names(l$desc)
names(comp) <- names(l$desc)

# add weights
for ( i in names(comp) ) {
  # pseudo-BMA weights
  comp[[i]][,"pseudobma"] <- rbind(
    pseudobma_weights( lpd_point[[i]], BB = F )[[1]],
    pseudobma_weights( lpd_point[[i]], BB = F )[[2]]
  )[,1]
  # pseudo-BMA+ weights (with Bayesian bootstrap)
  comp[[i]][,"pseudobma+"] <- rbind(
    pseudobma_weights( lpd_point[[i]], BB = T )[[1]],
    pseudobma_weights( lpd_point[[i]], BB = T )[[2]]
  )[,1]
  # stacking weigts
  comp[[i]][,"stacking"] <- rbind(
    stacking_weights( lpd_point[[i]] )[[1]],
    stacking_weights( lpd_point[[i]] )[[2]]
  )[,1]
}

# round the results
for ( i in names(comp) ) comp[[i]] <- comp[[i]] %>%
  mutate_if( is.numeric, round, 3 )

# prepare colors to use in graphs (rangi2 and a colorblind palette)
rangi2 <- "#8080FF"
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# loop through each model to create PSIS-LOO Pareto k plots
f.s1 <- list()
for ( i in names(l) ) {
  for ( j in names(l[[i]]) ) {
    f.s1[[i]][[j]] <- l[[i]][[j]]$diagnostics$pareto_k %>%
      as.data.frame() %>%
      rownames_to_column( "Data point") %>%
      mutate_if( is.character, as.numeric ) %>%
      rename( "Pareto shape k" = ".") %>%
      ggplot( aes(x = `Data point`, y = `Pareto shape k`) ) +
      geom_hline(
        yintercept = c(0, .5, .7) ,
        linetype = c("dotted", "dotdash", "dashed"),
        color = c( "black", "red", "red" )
      ) +
      geom_point( shape = 3, color = rangi2, size = 2, stroke = 1.5 ) +
      scale_x_continuous(
        breaks = seq(0,350,50) , labels = seq(0,350,50)
      ) +
      scale_y_continuous(
        limits = c(-.17,.7), breaks = seq(0,.7,.1),
        labels = sprintf( "%.1f", round( seq(0,.7,.1), 1) )
      )
  }
}

# arrange the plots
( f.s1$desc$no_cov + f.s1$desc$w_cov ) /
( f.s1$prog$no_cov + f.s1$prog$w_cov ) + 
plot_annotation( tag_levels = "A" )

# save the grid as Fig.S1
ggsave(
  "figures/Fig.s1 Pareto k diagnostics.png" ,
  height = 2 * 6.07 , width = 2 * 11.5 , dpi = "retina"
)

# ----------- visualization of models' posteriors  -----------

# prepare list for posteriors
post <- list(
  no_cov = list( desc = list(), prog = list() ),
  w_cov = list( desc = list(), prog = list() )
)

# first extract all parameters
for ( i in names(post) ) {
  for ( j in names(post[[i]]) ) {
    post[[i]][[j]]$pars <- m[[j]][[i]] %>%
      # fixed-effects, variance-covariance components,
      # other distributional parameters
      spread_draws(
        `b.*|sd.*|cor.*|sigma.*|nu.*` , regex = T
      )
    # keep only DRS related parameters in w_cov models
    if ( i == "w_cov" ) post[[i]][[j]]$pars <- post[[i]][[j]]$pars %>%
        select( contains("drs") )
  }
}

# seperate them by categories
# at the same time, back-transform to DRS-2 scale where appropriate
# always put it into a single column to make it ready for plotting
for ( i in names(post) ) {
  for ( j in names(post[[i]]) ) {
    # the global Intercept
    post[[i]][[j]]$intercept <- post[[i]][[j]]$pars %>% select( starts_with("b_") & contains("Intercept") )
    # other pre-surgery effects
    post[[i]][[j]]$preop <- post[[i]][[j]]$pars %>% select( starts_with("b_") & !contains("Intercept") & !contains("time") )
    # time-varying fixed-effects
    post[[i]][[j]]$time <- post[[i]][[j]]$pars %>% select( starts_with("b_") & contains("time") )
    # covariates
    # using a cheat code here, they're both imputed so they start with "bsp" instead of "b_"
    post[[i]][[j]]$cov <- post[[i]][[j]]$pars %>% select( starts_with("bsp") )
    # random-effects standard deviations
    post[[i]][[j]]$sds <- post[[i]][[j]]$pars %>% select( starts_with("sd_") )
    # random-effects correlations
    post[[i]][[j]]$cor <- post[[i]][[j]]$pars %>% select( starts_with("cor_") )
    # residual sigma
    post[[i]][[j]]$sigma <- post[[i]][[j]]$pars %>% select( contains("sigma") )
    # residual nu
    post[[i]][[j]]$nu <- post[[i]][[j]]$pars %>% select( contains("nu") )
    
    # after done with splitting, remove the original $pars to spare space
    post[[i]][[j]]$pars <- NULL
    
  }
}

# back-transform posteriors to DRS-2 scale where appropriate
# pivot longer to get the same format everywhere ready for plotting
for ( i in names(post) ) {
  for ( j in names(post[[i]]) ) {
    for ( k in names(post[[i]][[j]]) ) {
      # keep the null list alone and continue on
      if ( length( post[[i]][[j]][[k]] ) == 0 ) next
      # otherwise, pivot and back-transform (where appropriate)
      post[[i]][[j]][[k]] <- post[[i]][[j]][[k]] %>%
        # pivoting
        pivot_longer( cols = everything() ) %>%
        # renaming and back-transforming
        mutate(
          # renaming
          name = case_when(
            # posteriors from w_cov models
            i == "w_cov" & k == "cov" ~ sub( "bsp_drs_mi" , "" , name ),
            i == "w_cov" & k %in% c("intercept","preop","time") ~ sub( "b_drs_" , "" , name ),
            i == "w_cov" & k == "sds" ~ sub( "sd_id__drs_", "", name),
            # posteriors from no_cov models
            i == "no_cov" & k %in% c("intercept","preop","time") ~ sub( "b_" , "" , name),
            i == "no_cov" & k == "sds" ~ sub( "sd_id__", "", name),
            # residual parameters across models
            k == "cor" ~ "cor",
            k == "sigma" ~ "sigma",
            k == "nu" ~ "nu"
          ),
          # back-transforming
          value = case_when(
            k == "intercept" ~ value * scl$SD$drs + scl$M$drs,
            k %in% c("preop","time","cov","sds") ~ value * scl$SD$drs,
            k %in% c("cor","sigma","nu") ~ value
          )
        )
    }
  }
}




# 2022-05-01 preliminary plotting for github (well, for me of course)

# plotting, basic
data.frame(
  par = post$no_cov$desc$intercept$name,
  desc = post$no_cov$desc$intercept$value,
  prog = post$no_cov$prog$intercept$value
) %>% pivot_longer(
  cols = c("desc","prog"), values_to = "drs-2", names_to = "model"
) %>%
  ggplot( aes( x = `drs-2`, fill = model, color = model) ) +
  geom_density( size = 1 , alpha = .5 ) +
  scale_color_manual( values = cbPal[2:3]) +
  scale_fill_manual( values = cbPal[2:3])

# plotting, advanced
data.frame(
  par = post$no_cov$prog$time$name,
  no_cov = post$no_cov$prog$time$value,
  w_cov = post$w_cov$prog$time$value
) %>% pivot_longer(
  cols = c("no_cov","w_cov"), values_to = "drs-2", names_to = "model"
) %>%
  ggplot( aes( y = par, x = `drs-2`, fill = model, color = model ) ) +
  stat_halfeye( geom = "slab", slab_linetype = "solid" , slab_size = 1,  slab_alpha = .5 ) +
  geom_vline( xintercept = 0 , linetype = "dashed" ) +
  scale_color_manual( values = cbPal[2:3]) +
  scale_fill_manual( values = cbPal[2:3])

# an interesting one, looking how much did adding prognostic factors reduce patient-level variability
data.frame(
  par = post$no_cov$prog$sds$name,
  desc = post$no_cov$desc$sds$value,
  prog = post$no_cov$prog$sds$value
) %>% pivot_longer(
  cols = c("desc","prog"), values_to = "DRS-2 (SD)", names_to = "model"
) %>%
  ggplot( aes( y = par, x = `DRS-2 (SD)`, fill = model, color = model ) ) +
  stat_halfeye( geom = "slab", slab_linetype = "solid" , slab_size = 1,  slab_alpha = .6 ) +
  scale_color_manual( values = cbPal[2:3]) +
  scale_fill_manual( values = cbPal[2:3])
