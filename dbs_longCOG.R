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
  "rsvg", "DiagrammeRsvg", # for saving plots crated in DiagrammeR
  "brms", # for Bayesian model fitting / interface with STAN
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

# merge longitudinal d1 with baseline factor scores from efa[[4]] (joining by id)
d1 <- d1 %>% left_join( cbind.data.frame(id = d2$id, efa[[4]]$scores ) , by = "id" )

# ----------- total effect of time  -----------


# ----------- direct effect of time  -----------


# ----------- prognostic pre-surgery cognitive profile  -----------