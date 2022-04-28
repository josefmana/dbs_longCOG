# All analyses reported in the article (Mana et al., in review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 19043).

# I used the following versions of packages employed: dplyr_1.0.7, tidyverse_1.3.0, DiagrammeR_1.0.9, DiagrammeRsvg_0.1,
# rsvg_2.3.1, brms_2.16.3, psych_2.2.3, tidybayes_2.3.1, ggplot2_3.3.3 and patchwork_1.1.1

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "tidyverse", # for pivot_longer
  "DiagrammeR", # for flowcharts
  "DiagrammeRsvg", # for saving plots crated in DiagrammeR
  "rsvg", # for saving plots crated in DiagrammeR
  "brms", # for Bayesian model fitting / interface with Stan
  "psych", # for EFA
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
flowchart <- " digraph {
  
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
grViz(flowchart) %>% export_svg %>% charToRaw %>% rsvg_png("figures/Fig.1 inclusion-exclusion flowchart.png")


# ----------- directed acyclic graph -----------


# ----------- sample description  -----------

# list all stimulation parameters
pars <- names(d1)[which(names(d1)=="current_r_mA"):which(names(d1)=="frequency_l_Hz")]

# prepare dataframe to fill-in
tab1 <- data.frame(
  Md = rep( NA , length(pars) ), `Min-Max` = NA, M = NA, SD = NA, row.names = pars
)

# fill-in with statistics of all stimulation parameters
for ( i in pars ) {
  # median, rounded to hundredths
  tab1[ i , "Md" ] <- sprintf(
    "%.2f" , round( median( d1[[i]], na.rm = T ) , digits = 2 )
  )
  # min-max, rounded to decimals
  tab1[ i , "Min.Max" ] <- paste0(
    sprintf( "%.1f" , round( min( d1[[i]], na.rm = T ) , digits = 1 ) ), "-",
    sprintf( "%.1f" , round( max( d1[[i]], na.rm = T ) , digits = 1 ) )
  )
  # mean, rounded to hundredths
  tab1[ i , "M" ] <- sprintf(
    "%.2f" , round( mean( d1[[i]], na.rm = T ) , digits = 2 )
  )
  # SD, rounded to hundredths
  tab1[ i , "SD" ] <- sprintf(
    "%.2f" , round( sd( d1[[i]], na.rm = T ) , digits = 2 )
  )
}

# list all variables that will be included in Tab. 2 (baseline characteristics)
vars <- names(d2)[which(names(d2)=="age_stim_y"):which(names(d2)=="fp_dr")]

# prepare a dataframe to fill-in
tab2 <- data.frame(
  N = rep( NA, length(vars) ), Md = NA, `Min-Max` = NA, M = NA, SD = NA, row.names = vars
)

# fill-in with statistics of all variables but sex (which is nominal)
for ( i in vars[-3] ) {
  # number of data points
  tab2[ i , "N" ] <- sum( !is.na(d2[[i]]) )
  # median, rounded to integers
  tab2[ i , "Md" ] <- sprintf(
    "%.0f" , round( median(d2[[i]], na.rm = T ) , digits = 0 )
  )
  # min-max, rounded to integers
  tab2[ i , "Min.Max" ] <- paste0(
    sprintf( "%.0f" , round( min(d2[[i]], na.rm = T ) , digits = 0 ) ), "-",
    sprintf( "%.0f" , round( max(d2[[i]], na.rm = T ) , digits = 0 ) )
  )
  # mean, rounded to hundredths
  tab2[ i , "M" ] <- sprintf(
    "%.2f" , round( mean(d2[[i]], na.rm = T ) , digits = 2 )
  )
  # SD, rounded to hundredths
  tab2[ i , "SD" ] <- sprintf(
    "%.2f" , round( sd(d2[[i]], na.rm = T ) , digits = 2 )
  )
}

# add frequency of males to the table sex row
tab2[ "sex" , c("N","Min.Max")] <- c(
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
# need to use the complete.cases command because for three patient I have duplicated
# rows due to more than one stimulation parameter
f2.hist <- d1[ complete.cases(d1$drs_tot) , ] %>% 
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
f2.bin <- table( d1[ complete.cases(d1$drs_tot) , ]$id ) %>%
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
f2.hist / f2.bin + plot_annotation( tag_levels = "A" )

# save as Fig.2
ggsave( "figures/Fig.2 distribution of assessments.png" , height = 2 * 6.12 , width = 1.5 * 11.6 , dpi = "retina" )


# ----------- pre-surgery cognitive profile  -----------

# for EFA keep only id and cognitive tests in d2
d2 <- d2[ , c( 2, which(names(d2) == "tmt_a"):which(names(d2) == "fp_dr"), which(names(d2) %in% paste0("staix",1:2)) ) ]

# parallel test for number of factors
fa.parallel( d2[ , -1 ] ) # suggests 4-6 factors and 4 components

# fit EFA to get some meaningful latent cognitive factors from pre-surgery neuropsychology
efa <- lapply(
  X = 3:8, # fit 3-8 factor solutions
  FUN = function(i)
    fa(
      d2[ -1 ], # use all pre-surgery cognitive test, drop id
      nfactors = i, # fit 3-8 factor solutions
      rotate = "varimax", # rotate via Varimax (to enforce orthogonality)
      scores = "regression", # add regression scores inferred to each patient from the factor solutions
      missing = T , impute = "median" # impute missing values with a median, the weakest part of the whole analysis
    )
)

# create a table for performance indexes of each factor solution
fat <- data.frame(
  no_factors = 3 : (2+length(efa)), # identifiers
  TLI = NA, # Tucker-Lewis index, > 0.9 is considered good
  BIC = NA, # Bayesian Information Criterion, lower values (including negative) are considered better
  RMSEA = NA, # root-mean square error of approximation, lower is better, < 0.08 is considered good
  RMSEA_90_CI = NA, # 90% CI for RMSEA
  var_account = NA # total variance "accounted for" by all included factors
)

# fill-in the fat table
for (i in 1:length(efa)) {
  fat[i, ] <- c(
    i+2, # number of factors (because I fitted 3:8 factor solutions and this loop goes through i = 1:6, it's i+2)
    sprintf( "%.3f" , round(efa[[i]]$TLI, 3) ), # Tucker-Lewis index
    sprintf( "%.3f" , round(efa[[i]]$BIC, 3) ), # Bayesian Information Criterion
    sprintf( "%.3f" , round(efa[[i]]$RMSEA[1], 3) ), # root-mean square error of approximation
    paste0(
      "[", sprintf( "%.3f", round(efa[[i]]$RMSEA[2], 3) ), ", ",
      sprintf( "%.3f", round(efa[[i]]$RMSEA[3], 3) ), "]"
    ), # root-mean square error of approximation, 90% CI
    sprintf( "%.3f" , round(efa[[i]]$Vaccounted["Cumulative Var",i+2], 3) ) # total variance "accounted for"
  )
}

# 6-factor solution is leading by my inspection of fat table:
# 1) lowest BIC, 2) everything else solid w.r.t. competing  models
# However, judging by TLI and RMSEA, the fit ain't ideal for either model, a limitation to mention
# Check loading patterns of 5-8 factor solutions to see if they make sense
print( efa[[3]]$loadings, cutoff = .4, sort = T ) # reasonable loadings but really bad stats
print( efa[[4]]$loadings, cutoff = .4, sort = T ) # reasonable loadings, acceptable stats (w.r.t. competition)
print( efa[[5]]$loadings, cutoff = .4, sort = T ) # MR7 includes only one test with loading > 0.4 (Sim.)
print( efa[[6]]$loadings, cutoff = .4, sort = T ) # ok stats but MR7 includes only one test with loading > 0.4 (Sim.)

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
    anxiety = mean( d1[d1$ass_type=="pre", ]$anxiety , na.rm = T ),
    # pre-surgery prognostic factors/clinical variables
    dur_stim = mean( d1[d1$ass_type=="pre", ]$dur_stim_y , na.rm = T ), # 11.67
    age_stim = mean( d1[d1$ass_type=="pre", ]$age_stim_y , na.rm = T ), # 57.25
    led_pre = mean( d1[d1$ass_type=="pre", ]$ledd_pre_mg , na.rm = T ), # 1696.88
    updrsiii_pre = mean( d1[d1$ass_type=="pre", ]$mds_updrs_iii_med_off , na.rm = T ), # 45.79
    ldopa_resp = mean( d1[d1$ass_type=="pre", ]$ldopa_resp , na.rm = T ) # 52.64
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
    anxiety = sd( d1[d1$ass_type=="pre", ]$anxiety , na.rm = T ), # 0.90
    # pre-surgery prognostic factors/clinical variables
    dur_stim = sd( d1[d1$ass_type=="pre", ]$dur_stim_y , na.rm = T ), # 4.05
    age_stim = sd( d1[d1$ass_type=="pre", ]$age_stim_y , na.rm = T ), # 7.96
    led_pre = sd( d1[d1$ass_type=="pre", ]$ledd_pre_mg , na.rm = T ), # 672.33
    updrsiii_pre = sd( d1[d1$ass_type=="pre", ]$mds_updrs_iii_med_off , na.rm = T ), # 10.93
    ldopa_resp = sd( d1[d1$ass_type=="pre", ]$ldopa_resp , na.rm = T ) # 12.81
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
  #led = ( ledd_mg - scl$M$led ) / scl$SD$led,
  sex = as.factor( sex ), # for better estimation of BDI in the second and third models
  cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ), # right censoring for DRS == 144
  cens_bdi = ifelse( bdi == min(bdi, na.rm = T) , "left" , "none" ), # left censoring for BDI == 0
  # pre-surgery prognostic factors/cognitive profile, scale such that higher value imply deficit
  ment_flex = -( ment_flex - scl$M$ment_flex ) / scl$SD$ment_flex,
  epis_mem = -( epis_mem - scl$M$epis_mem ) / scl$SD$epis_mem,
  visp_wm = -( visp_wm - scl$M$visp_wm ) / scl$SD$visp_wm,
  visp_mem = -( visp_mem - scl$M$visp_mem ) / scl$SD$visp_mem,
  verb_wm = -( verb_wm - scl$M$verb_wm ) / scl$SD$verb_wm,
  anxiety = -( anxiety - scl$M$anxiety ) / scl$SD$anxiety,
  # pre-surgery prognostic factors/clinics
  dur_stim = ( dur_stim_y - scl$M$dur_stim ) / scl$SD$dur_stim,
  age_stim = ( age_stim_y - scl$M$age_stim ) / scl$SD$age_stim,
  led_pre = ( ledd_pre_mg - scl$M$led_pre ) / scl$SD$led_pre,
  updrsiii_pre = ( mds_updrs_iii_med_off - scl$M$updrsiii_pre ) / scl$SD$updrsiii_pre,
  ldopa_resp = -( ldopa_resp - scl$M$ldopa_resp ) / scl$SD$ldopa_resp
) %>% select(
  # keep only variables of interest
  id, time, drs, cens_drs, bdi, cens_bdi, sex, # outcomes, time and sex
  ment_flex, epis_mem, visp_wm, visp_mem, verb_wm, anxiety, # pre-surgery prognostic factors/cognitive profile
  dur_stim, age_stim, led_pre, updrsiii_pre, ldopa_resp # pre-surgery prognostic factors/clinics
)


# ----------- total effect of time  -----------

# prior predictive check

# set-up the linear model
f.drs1 <- bf( drs | cens(cens_drs) ~ 1 + time + (1 + time | id) )

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
m1 <- brm(
  formula = f.drs1, family = student(), prior = p1, data = d3,
  seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/total_effect.rds"
)

# add LOO and WAIC criteria
add_criterion( m1 , criterion = c("loo","waic") )

# save the priors for documentation purposes
pr <- list( glmm1 = prior_summary(m1) )

# some posterior checks
max( rhat(m1) ) # 1.007
max( loo(m1)$diagnostics$pareto_k ) # 0.53
plot( loo(m1) ) # three cases > 0.5, no case > 0.7

# ----------- direct effect of time  -----------

# or rather depression confounded via LEDD effect?

# prior predictive simulation

# set-up the linear models
f.drs2 <- bf( drs | cens(cens_drs) ~ 1 + time + mi(bdi) + (1 + time | p | id) )
f.bdi2 <- bf( bdi | cens(cens_bdi) + mi() ~ 1 + time + sex + (1 + time | p | id) )

# set contrast for sex
contrasts(d3$sex) <- -contr.sum(2)/2 # female = -0.5, male = 0.5

# set-up priors
p2 <- c(
  # fixed-effects
  prior( normal(0.3, .1), class = Intercept , resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time , resp = drs  ),
  prior( normal(0, .1), class = b, coef = mibdi, resp = drs ),
  prior( normal(.6, .5), class = Intercept, resp = bdi ),
  prior( normal(0, .5), class = b, coef = time, resp = bdi ),
  prior( normal(0, .5), class = b, coef = sex1, resp = bdi ),
  # random-effects
  prior( normal(0, .1), class = sd, coef = Intercept, group = id , resp = drs  ),
  prior( normal(0, .1), class = sd, coef = time, group = id , resp = drs  ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  # variances
  prior( exponential(1), class = sd , resp = drs  ),
  prior( exponential(1), class = sigma , resp = drs  ),
  prior( exponential(1), class = sd, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  # covariances
  prior( lkj(2), class = cor ),
  # "degrees of freedom"/"normality parameter" nu
  prior( gamma(2, 0.1), class = nu , resp = drs  ),
  prior( gamma(2, 0.1), class = nu, resp = bdi )
)

# fit the "direct" effect of time GLMM
m2 <- brm(
  formula = f.drs2 + f.bdi2 + set_rescor(F), family = student(), prior = p2, data = d3,
  seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/direct_effect.rds"
)

# add LOO and WAIC criteria
add_criterion( m2 , criterion = c("loo","waic") , resp = "drs" )
add_criterion( m2 , criterion = c("loo","waic") , resp = "bdi" , newdata = d3[ complete.cases(d3$bdi), ] )

# save the priors for documentation purposes
pr$glmm2 <- prior_summary(m2)

# some posterior checks
max( rhat(m2) ) # 1.008
max( loo(m2, resp = "drs")$diagnostics$pareto_k ) # 0.53
max( loo(m2, resp = "bdi", newdata = d3[ complete.cases(d3$bdi), ])$diagnostics$pareto_k ) # 0.68
par( mfrow = c(2,1) )
plot( loo(m2, resp = "drs") ) # four cases > 0.5, no case > 0.7
plot( loo(m2, resp = "bdi", newdata = d3[ complete.cases(d3$bdi), ]$diagnostics$pareto_k) ) # seven cases > 0.5


# ----------- prognostic pre-surgery cognitive profile  -----------

# prior predictive simulation

# set-up the linear models
f.drs3 <- bf(
  # covariates
  drs | cens(cens_drs) ~ 1 + time +  mi(bdi) +
    # demographic/clinic predictors
    age_stim*time + mi(dur_stim)*time + mi(led_pre)*time + mi(updrsiii_pre)*time + mi(ldopa_resp)*time +
    # pre-surgery cognitive profile predictors
    ment_flex*time + epis_mem*time + visp_wm*time + visp_mem*time + verb_wm*time + anxiety*time +
    # patient-specific effects
    (1 + time | p | id)
)

# predictors with missing values
f.bdi3 <- bf( bdi | cens(cens_bdi) + mi() ~ 1 + time + sex + (1 + time | p | id) )
f.dur3 <- bf( dur_stim | mi() ~ 1 + age_stim )
f.ledpre3 <- bf( led_pre | mi() ~ 1 + mi(dur_stim) )
f.updrsiii3 <- bf( updrsiii_pre | mi() ~ 1 + sex + age_stim + mi(dur_stim) )
f.ldoparesp3 <- bf( ldopa_resp | mi() ~ 1 )

# set-up priors
p3 <- c(
  # DRS-2
  prior( normal(0.3, .1), class = Intercept , resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time , resp = drs  ),
  prior( normal(0, .1), class = b, coef = mibdi, resp = drs ),
  prior( normal(0, .1), class = b, coef = age_stim, resp = drs ),
  prior( normal(0, .1), class = b, coef = midur_stim, resp = drs ),
  prior( normal(0, .1), class = b, coef = miled_pre, resp = drs ),
  prior( normal(0, .1), class = b, coef = miupdrsiii_pre, resp = drs ),
  prior( normal(0, .1), class = b, coef = mildopa_resp, resp = drs ),
  prior( normal(0, .1), class = b, coef = ment_flex, resp = drs ),
  prior( normal(0, .1), class = b, coef = epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = anxiety, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:age_stim, resp = drs ),
  prior( normal(0, .1), class = b, coef = midur_stim:time, resp = drs ),
  prior( normal(0, .1), class = b, coef = miled_pre:time, resp = drs ),
  prior( normal(0, .1), class = b, coef = miupdrsiii_pre:time, resp = drs ),
  prior( normal(0, .1), class = b, coef = mildopa_resp:time, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:ment_flex, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_mem, resp = drs ),
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
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  prior( exponential(1), class = sd, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  prior( gamma(2, 0.1), class = nu, resp = bdi ),
  # group-level correlation matrix
  prior( lkj(2), class = cor ),
  # pre-surgery disease duration
  prior( student_t(3, -0.2, 2.5) , class = Intercept, resp = durstim),
  prior( normal(0, .5), class = b, coef = age_stim, resp = durstim ),
  prior( exponential(1), class = sigma, resp = durstim ),
  prior( gamma(2, 0.1), class = nu, resp = durstim ),
  # pre-surgery LED
  prior( student_t(3, -0.1, 2.5), class = Intercept, resp = ledpre ),
  prior( normal(0, .5), class = b, coef = midur_stim, resp = ledpre ),
  prior( exponential(1), class = sigma, resp = ledpre ),
  prior( gamma(2, 0.1), class = nu, resp = ledpre ),
  # pre-surgery UPDRS-III (off medication)
  prior( student_t(3, 0, 2.5), class = Intercept, resp = updrsiiipre ),
  prior( normal(0, .5), class = b, coef = sex1, resp = updrsiiipre ),
  prior( normal(0, .5), class = b, coef = age_stim, resp = updrsiiipre ),
  prior( normal(0, .5), class = b, coef = midur_stim, resp = updrsiiipre ),
  prior( exponential(1), class = sigma, resp = updrsiiipre ),
  prior( gamma(2, 0.1), class = nu, resp = updrsiiipre ),
  # pre-surgery levodopa response
  prior( student_t(3, 0.2, 2.5), class = Intercept, resp = ldoparesp ),
  prior( exponential(1), class = sigma, resp = ldoparesp ),
  prior( gamma(2, 0.1), class = nu, resp = ldoparesp )
)

# fit the "prognostic" model
m3 <- brm(
  formula = f.drs3 + f.bdi3 + f.dur3 + f.ledpre3 + f.updrsiii3 + f.ldoparesp3 + set_rescor(F),
  family = student(), prior = p3, data = d3,
  seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/prognostic.rds"
)

# add LOO and WAIC criteria
add_criterion( m3 , criterion = c("loo","waic") , resp = "drs" )
add_criterion( m3 , criterion = c("loo","waic") , resp = "bdi" , newdata = d3[ complete.cases(d3$bdi), ] )

# save the priors for documentation purposes
pr$glmm3 <- prior_summary(m3)

# some posterior checks
max( rhat(m3) ) # 1.01
max( loo(m3, resp = "drs")$diagnostics$pareto_k ) # 0.47
max( loo(m3, resp = "bdi", newdata = d3[ complete.cases(d3$bdi), ])$diagnostics$pareto_k ) # 0.66
par( mfrow = c(2,1) )
plot( loo(m3, resp = "drs") ) # all cases < 0.5
plot( loo(m3, resp = "bdi", newdata = d3[ complete.cases(d3$bdi), ]$diagnostics$pareto_k) ) # nine cases > 0.5