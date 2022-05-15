# All analyses reported in the article (Mana et al., in review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 19043).

# I used the following versions of packages employed: dplyr_1.0.7, tidyverse_1.3.0, DiagrammeR_1.0.9, DiagrammeRsvg_0.1,
# rsvg_2.3.1, missMDA_1.18, psych_2.2.3, brms_2.16.3, loo_2.4.1, tidybayes_2.3.1, ggplot2_3.3.3 and patchwork_1.1.1

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
  #"loo", # for PSIS-LOO based operations
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

# set number of multiple imputations to account for missing pre-surgery data
imp = 100
s = 87542 # seed for reproducibility

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# read the data set and prepare subsets for individual analyses
d0 <- read.csv( "data/20220508_dbs_longCOG_data.csv" , sep = "," )
d1 <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data
d2 <- d1[ d1$ass_type == "pre" , ] # only pre-surgery assessments of included patients

# read a file with mapping variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 )


# ----------- participant inclusion flowchart  -----------

# check that when selecting only rows containing pre-surgery assessment there ain't no patient duplicated
isTRUE(
  all.equal(
    d0[ d0$ass_type == "pre" , ]$id,
    unique( d0[ d0$ass_type == "pre" , ]$id )
  )
) # TRUE

# create a table summarizing reasons for excluding patients
t0 <- table(d0[ d0$ass_type == "pre" , ]$why_excluded)
print(t0)

# using numbers from t0 create an inclusion/exclusion flowchart
f1 <- " digraph {
  
  /// define nodes which will include numbers reflecting the inclusion/exclusion process
  node [ fontname = Calibri, fontsize = 24, shape = box, style = rounded , width = 2.5 , margin= 0.25 , penwidth = 2 ];
  
  /// create a box for all patients, i.e., sum(t) = 200 patients
  all_pats [ label =<
  <b>200 consecutive PD patients </b><br/><br/>Local database 2000-2020<br/>General University Hospital<br/>in Prague
  >];
         
  /// create a box for all STN-DBS patients, i.e., sum(t[c(1,4,5,6,7,8,9,11)])) = 173 patients
  all_stn [ label =<
  <b>173 patients </b><br/><br/>implanted with STN-DBS
  >];
  
  /// create a box for all included patients, i.e., t[4] = 126 patients
  all_incl [ label =<
  <b>126 patients </b><br/><br/>followed-up longitudinally
  >];
  
  /// create nodes for exluded patients, specifying only these characteristics that will differ from the nodes above
  node [ fixedsize = T, width = 5.5, height = 2.2 ];
  
  /// create a box for non-STN-DBS patients, i.e., t[c(3,13,2,10,12)] with sum(t[c(3,13,2,10,12)]) = 27 patients
  excl_nostn [ label =<
  <b>27 patients excluded due to</b><br align = 'left'/><br align = 'left'/>
  12 GPi-DBS<br align = 'left'/>
  4 VIM-DBS<br align = 'left'/>
  4 duodopa<br align = 'left'/>
  5 rejected<br align = 'left'/>
  2 suspended<br align = 'left'/>
  >];
  
  /// create a box for STN-DBS excluded patients, i.e., t[c(1,7,9, 8, 6,11, 5)], sum(t[c(1,7,9,8,6,11,5)]) = 47 patients
  excl_stn [ label =<
  <b>47 patients excluded due to</b><br align = 'left'/><br align = 'left'/>
  24 pre-surgery data missing<br align = 'left'/>
  18 follow-up data missing<br align = 'left'/>
  3 unilateral STN-DBS<br align = 'left'/>
  2 not speaking Czech<br align = 'left'/>
  >];
  
  /// create dummy nodes for horizontally forking out of the vertical 'inclusion flow' to excluded sides
  node [ shape = rectangle, width = 0, height = 0, label = '', fill = black ];
  
  /// create directed edges in the inclusion (from dummy/fork vertically to inclusion boxes)
  /// and the exlusion (from forks horizontally to exclusion boxes) parts of the flowchart
  /// first make the arrows bigger
  edge [ arrowsize = 2.5, penwidth = 2.5 ]
  
  /// specifiy paths
  fork1 -> all_stn; fork2 -> all_incl;
  
  /// for the horizontal paths use 'rank = same' to ensure their nodes are level
  { rank = same ; fork1 -> excl_nostn }
  { rank = same ; fork2 -> excl_stn }
  
  /// create non-directed edges from the inclusion boxes to dummy/fork boxes (vertically)
  edge [ dir = none, penwidth = 2.5 ]
  all_pats -> fork1; all_stn -> fork2;

  /// seperate dummy/fork nodes from exclusion boxes by some reasonable distance (in inches)
  nodesep = 1
  }"

# save the flowchart as Fig 1
grViz(f1) %>% export_svg %>% charToRaw %>% rsvg_png("figures/Fig 1 inclusion-exclusion flowchart.png")


# ----------- sample description  -----------

# list all stimulation parameters
pars <- names(d1)[which(names(d1)=="current_r_mA"):which(names(d1)=="frequency_l_Hz")]

# prepare a data frame to be filled-in with stimulation parameters' summary
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

# list all variables that will be included in Tab 2 (baseline characteristics)
vars <- names(d2)[which(names(d2)=="age_stim_y"):which(names(d2)=="fp_dr")]

# prepare a data frame to be filled-in with baseline characteristics
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

# prepare a histogram of the distribution of assessments across time (Fig 2A)
f2 <- list()

# need to use the complete.cases command because three patient have duplicated rows
# due to more than one stimulation parameter
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

# prepare a bin plot showing distribution of the number of assessments per patient (Fig 2B)
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

# arrange Fig 2A and Fig 2B for printing
f2$hist / f2$bin + plot_annotation( tag_levels = "A" )

# save as Fig 2
ggsave( "figures/Fig 2 distribution of assessments.png" , height = 2.5 * 6.12 , width = 1.5 * 11.6 , dpi = "retina" )


# ----------- pre-surgery cognitive profile  -----------

# for EFA keep only id and cognitive tests in d2
d2 <- d2[ , c( 2, which(names(d2) == "tmt_a"):which(names(d2) == "fp_dr"), which(names(d2) %in% paste0("staix",1:2)) ) ]

# log-transform reaction times before the analysis
for ( i in c( paste0("tmt_", c("a","b")), paste0("pst_", c("d","w","c")) ) ) d2[[i]] <- log( d2[[i]] )

# find out the optimal number of components for multiple imputation
nb <- estim_ncpPCA( d2[,-1] , ncp.min = 0, ncp.max = 10 , nbsim = imp )

# impute via PCA-based multiple imputation (n = 100 imputations)
set.seed(s) # set seed for reproducibility
d2.imp <- MIPCA( d2[ ,- 1] , ncp = nb$ncp )

# calculate parallel test for each imputed data set 
set.seed(s) # set seed for reproducibility
p.test <- lapply( 1:100 , function(i) fa.parallel( d2.imp$res.MI[[i]] ) )

# look at the results of parallel tests
table( sapply( 1:100 , function(i) p.test[[i]]$nfact ) )

# fit EFA to each imputed data set with 3:8 factors
efa <- lapply(
  # loop through all 100 data sets
  1:imp, function(i)
    lapply(
      # loop through three to eight latent factors for each imputation
      3:8, function(j)
        fa(
          d2.imp$res.MI[[i]], # one-by-one use each imputed data set
          nfactors = j, # fit 3-8 factor solutions
          rotate = "varimax", # rotate varimax to enforce orthogonality and for interpretation purposes
          scores = "regression" # compute regression scores for each patient
        )
    )
)

# prepare a table for performance indexes of each solution for each imputed data set
fat <- array(
  # create an empty 6 (factor solutions) x 5 (performance indexes) x 100 (imputations) array
  data = NA, dim = c( length( efa[[1]] ), 5, imp ),
  # add dimension names
  dimnames = list(
    paste0( 1:length(efa[[1]])+2, "_factors" ), # number of factors
    c("TLI", "RMSEA", "RMSEA_90_CI_low", "RMSEA_90_CI_upp", "var_account"), # performance indexes
    1:imp # number of imputed data set
  )
)

# loop through the array and gather all performance indexes
for ( i in 1:imp ) {
  for ( j in 1:length(efa[[i]]) ) {
    fat[ j , , i ] <- c(
      round( efa[[i]][[j]]$TLI , 3 ), # Tucker-Lewis index, > 0.9 is considered good
      round( efa[[i]][[j]]$RMSEA[1], 3 ), # root-mean square error of approximation, < 0.08 is considered good
      round( efa[[i]][[j]]$RMSEA[2], 3 ), # 90% CI RMSEA, lower boundary
      round( efa[[i]][[j]]$RMSEA[3], 3 ), # 90% CI RMSEA, upper boundary (ideally should be < 0.08)
      round(efa[[i]][[j]]$Vaccounted["Cumulative Var",j+2], 3) # total variance "accounted for" by all included factors
    )
  }
}

# summarize the fat table
fat_sum <- data.frame( model = paste0( 3:8, "_factors"), TLI = NA, RMSEA = NA, RMSEA_90_CI_upp = NA, var_account = NA )

# fill-in summary of each model into fat_sum
for ( i in fat_sum$model ) {
  for ( j in names(fat_sum)[-1] ) {
    fat_sum[ which(fat_sum$model == i) , j ] <- paste0(
      sprintf( "%.3f" , round( mean( fat[ i, j, ] ), 3 ) ), " (",
      sprintf( "%.3f" , round( sd( fat[ i, j, ] ), 3 ) ), ")"
    )
  }
}

# prepare a list for RMSEA and TLI frequencies across all imputations
fat_perc <- list(
  # upper RMSEA lower than 0.08
  `RMSEA_upp < 0.08 (%)` = t( fat[ , "RMSEA_90_CI_upp" , ] ) %>%
    as.data.frame %>%
    mutate( across( everything() , ~ replace( . , . > .08 , NA) ) ),
  # TLI higher than 0.9
  `TLI > 0.90 (%)` = t( fat[ , "TLI" , ] ) %>%
    as.data.frame %>%
    mutate( across( everything() , ~ replace( . , . < .9 , NA) ) )
)

# calculate the percentages/frequencies
for ( i in names(fat_perc) ) fat_perc[[i]] <- sapply(
  names(fat_perc[[i]]), function(j) sum( complete.cases( fat_perc[[i]][[j]] ) )
)

# bind the percentages to the FA summary table (fat_sum)
fat_sum <- fat_sum %>% left_join(
  do.call( cbind.data.frame , fat_perc ) %>%
    # prepare a column "model" and bind the percentages to the fat_sum
    rownames_to_column( var = "model")
)

# visualize performance indexes across imputed data sets with density plots (Fig S1)
# set-up a list to contain Fig S1 component figures
f.s1 <- list()

# loop through TLI and upper 90% CI RMSEA
for ( i in c("TLI","RMSEA_90_CI_upp") ) {
  f.s1[[i]] <- fat[ , i , ] %>%
    # format the table for plotting
    t %>% as.data.frame %>%
    pivot_longer( everything() , names_to = "Model" , values_to = i ) %>%
    # mutate variables to an appropriate form
    mutate( Model = substr( gsub( "_", "-", Model ) , 1 , nchar(Model)-1 ) ) %>%
    rename( "index" = i ) %>%
    # plotting proper
    ggplot( aes( x = index , fill = Model) ) +
    geom_density( alpha = .4 , color = NA ) +
    scale_fill_manual( values = cbPal[3:8] ) + # use colorblind-friendly palette
    geom_vline( # add a vertical line depicting good performance heuristics
      linetype = "dashed", size = 1.2, xintercept = case_when( i == "TLI" ~ .9, i == "RMSEA_90_CI_upp" ~ .08 )
    ) +
    labs( y = "Density", x = case_when( i == "TLI" ~ "TLI", i == "RMSEA_90_CI_upp" ~ "RMSEA (upper 90% CI)" ) )
}

# arrange Fig S1 for printing
f.s1$TLI / f.s1$RMSEA_90_CI_upp + plot_layout( guides = "collect" ) + plot_annotation( tag_levels = "A" )

# save as Fig S1
ggsave(
  "figures/Fig S1 factor analysis performance indexes.png", height = 1.75 * 6.07, width = 1.75 * 11.5, dpi = "retina"
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
doms <- c(
  "proc_spd", # loaded on primarily by PST, the first factor in 83% data sets
  "epis_mem", # loaded on primarily by RAVLT, the second factor in 81% data sets
  "verb_wm", # loaded on primarily by DS, the third factor in 62% data sets
  "visp_mem", # loaded on primarily by FP, the fourth factor in 46% data sets
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


# write a summary of loadings across all 100 imputed data sets
t3 <- matrix(
  data = NA, nrow = nrow(loads[ , , 1]), ncol = ncol(loads[ , , 1]),
  dimnames = list( rownames(loads[ , , 1 ]) , colnames(loads[ , , 1]) )
) %>% as.data.frame()

# fill-in averages and SDs across all 100 imputations 
for ( i in rownames(t3) ) {
  t3[ i , ] <- paste0(
      sprintf( "%.2f" , round( loads[i, , ] %>% t %>% colMeans , 2 ) ), " (", # mean
      sprintf( "%.2f" , round( loads[i, , ] %>% t %>% apply( . , 2 , sd ), 2) ), ")" # SD
    )
}

# visualize loadings (i.e.,median loadings across imputations)
# prepare labels for the factors
lab <- as_labeller( c( "proc_spd" = "Processing\nspeed",
                       "epis_mem" = "Episodic\nmemory",
                       "verb_wm" = "Verbal\nworking memory",
                       "visp_mem" = "Visuospatial\nmemory",
                       "set_shift" = "Set shifting",
                       "anxiety" = "Anxiety",
                       "visp_wm" = "Spatial\nworking memory"
                       )
                    )

# create a median factor loading plot as Fig 3
# start by calculating median loading for each test/factor pairs across imputations
f3 <- sapply (
  rownames(t3)[ 1:(nrow(t3)-2) ],
  function(i) loads[ i , , ] %>% t %>% apply( . , 2 , median )
) %>% as.data.frame %>%
  # change to an appropriate format
  mutate( f.order = 1:nf ) %>% # prepare an order variable to sort tests on y-axis
  rownames_to_column( var = "Factor" )  %>% # make factors from row names to explicit variables
  pivot_longer( # change to a long format
    rownames(t3)[ 1:(nrow(t3)-2) ] , values_to = "Loading" , names_to = "Test"
  ) %>%
  mutate( t.order = rep(1:(nrow(t3)-2), nf ) ) %>% # prepare an order variable to sort factors
  # plot it
  ggplot( aes( x = reorder(Test, -t.order) , y = abs(Loading) , fill = Loading) ) +
    facet_wrap( ~ reorder(Factor, f.order) , nrow = 1 , labeller = lab ) + # place the factors to separate facets
    geom_bar( stat = "identity" ) + # create the bars
    coord_flip() + # flip axes so that tests are on the y-axis and loadings on x-axis
    scale_fill_gradient2( # define the fill color gradient: orange = positive, blue = negative
      name = "Loading", high = cbPal[2], mid = "white", low = cbPal[6], midpoint = 0, guide = F
    ) +
    scale_x_discrete( name = "Test" , labels = rev( var_nms[ rownames(t3)[ 1:(nrow(t3)-2) ] , ] )  ) +
    scale_y_continuous(
      name = "Loading Strength", limits = c(-0.05,1.05),
      labels = sprintf( "%.1f" , round( seq(0, 1, .4) , 2 ) ), breaks = seq(0, 1, .4) 
    ) +
  theme( panel.grid.major = element_line() ) # add grid lines for easier orientation

# save as Fig 3
ggsave( "figures/Fig 3 factor loadings.png", height = 1.5*11.8, width = 1.8*10.8, dpi = "retina" )

# save all EFA models
saveRDS( object = efa, file = "models/efa.rds" )

# ----------- pre-processing for longitudinal analyses  -----------

# before merging compute scaling values for DRS-2, BDI-II, LEDD, age and time
scl <- list(
  M = list(
    drs = mean( d1$drs_tot , na.rm = T ), # 136.90
    bdi = mean( d1$bdi , na.rm = T ), # 10.95
    led = mean( d1$ledd_mg , na.rm = T ), # 1197.21
    age = mean( d1$age_ass_y, na.rm = T ) # 59.63
  ),
  SD = list(
    drs = sd( d1$drs_tot , na.rm = T ), # 7.72
    bdi = sd( d1$bdi , na.rm = T ), # 7.19
    led = sd( d1$ledd_mg , na.rm = T ), # 687.20
    age = sd( d1$age_ass_y, na.rm = T ) # 8.27
  ),
  # add median time of pre-surgery assessment for GLMMs intercepts
  Md = list(
    time = -median( d1[d1$ass_type == "pre", ]$time_y , na.rm = T ) # 0.30
  )
)

# merge longitudinal d1 with baseline factor scores (joining by id)
d3 <- lapply(
  1:imp,
  function(i) d1 %>%
    # first prepare a pre-surgery df with id and factor scores for each patient
    left_join( cbind.data.frame(id = d2$id, efa[[i]][[nf-2]]$scores ) , by = "id" ) %>%
    filter( complete.cases(drs_tot) ) %>% # get rid of three dummy rows due to more than one stimulation parameter (no DRS-2)
    # scale all DRS-2, BDI-II, LEDD and time already
    mutate(
      time = time_y + scl$Md$time,
      drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
      bdi = ( bdi - scl$M$bdi ) / scl$SD$bdi,
      led = ( ledd_mg - scl$M$led ) / scl$SD$led,
      age = ( age_ass_y - scl$M$age ) / scl$SD$age,
      sex = as.factor( sex ), # for better estimation of BDI in the second (covariate) model
      cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ), # right censoring for DRS == 144
    ) %>%
    # keep only variables of interest
    select(
      id, time, drs, cens_drs, bdi, led, age, sex, # outcomes, demographics, clinics
      proc_spd, epis_mem, verb_wm, visp_mem, set_shift, anxiety, visp_wm # pre-surgery cognition
    )
)

# loop across all imputations to get means and SDs of the pre-surgery cognitive domains,
# then transform the pre-surgery cognition in each data set to (pre-surgery) zero mean, unit SD variables
for( i in doms ) {
  for ( j in 1:imp ) {
    # calculate scaling values
    scl$M[[i]][[j]] <- d3[[j]][[i]] %>% mean( na.rm = T )
    scl$SD[[i]][[j]] <- d3[[j]][[i]] %>% sd( na.rm = T )
    # scale in the jth imputed data set
    d3[[j]][[i]] <- case_when(
      # all but anxiety measures will be inverse such that parameters
      # can be interpreted as effect of deficit in said measure
      i == "anxiety" ~ ( d3[[j]][[i]] - scl$M[[i]][[j]] ) / scl$SD[[i]][[j]],
      i != "anxiety" ~ -( d3[[j]][[i]] - scl$M[[i]][[j]] ) / scl$SD[[i]][[j]]
    )
  }
}

# clean environment for space (will need it)
rm( list = ls()[ !( ls() %in% c( "cbPal", "d1", "d3", "doms", "s", "scl", "var_nms" ) ) ] )
gc() # garbage collection to ease RAM a bit if needed

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .99 # adapt_delta parameter


# ----------- no covariates prognostic model  -----------

# set-up the linear model
f1.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time || id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf

# set-up priors
p1 <- c(
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
  prior( exponential(10), class = sd, coef = Intercept, group = id ),
  prior( exponential(10), class = sd, coef = time, group = id ),
  prior( exponential(10), class = sd ),
  # other distributional parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu )
)

# fit the model with Student response function
brm_multiple(
  formula = f1.drs, family = student(), prior = p1, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m1_nocov.rds"
) # not creating an object to spare RAM


# ----------- prognostic model with flat priors  -----------

# the same linear model as for m1
f2.drs <- f1.drs

# will be using brms default priors
p2 <- NULL

# fit the model with Student response function
brm_multiple(
  formula = f2.drs, family = student(), prior = p2, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m2_flat_priors.rds"
)







# ----------- prognostic model with covariates -----------

# set-up the linear model
f2.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + age + mi(bdi) + mi(led) + ", # outcome and covariates
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time || id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf

# set-up linear models for covariates with missing values
f2.bdi <- bf( bdi | mi() ~ 1 + mi(led) + sex + time + (1 + time | id) )
f2.led <- bf( led | mi() ~ t2(time) + (1 | id) )

# set-up priors
p2 <- c(
  # fixed effects
  prior( normal(0.3, .1), class = Intercept, resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time, resp = drs ),
  prior( normal(0, .1), class = b, coef = proc_spd, resp = drs ),
  prior( normal(0, .1), class = b, coef = epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = set_shift, resp = drs ),
  prior( normal(0, .1), class = b, coef = anxiety, resp = drs ),
  prior( normal(0, .1), class = b, coef = visp_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:proc_spd, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:epis_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:verb_wm, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_mem, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:set_shift, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:anxiety, resp = drs ),
  prior( normal(0, .1), class = b, coef = time:visp_wm, resp = drs ),
  # random effects
  prior( exponential(10), class = sd, coef = Intercept, group = id ),
  prior( exponential(10), class = sd, coef = time, group = id ),
  prior( exponential(10), class = sd ),
  # other distributional parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu )
)

# fit the model with Student response function
m1 <- brm_multiple(
  formula = f1.drs, family = student(), prior = p1, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m1_nocov.rds"
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
  # stacking weights
  comp[[i]][,"stacking"] <- rbind(
    stacking_weights( lpd_point[[i]] )[[1]],
    stacking_weights( lpd_point[[i]] )[[2]]
  )[,1]
}

# round the results
for ( i in names(comp) ) comp[[i]] <- comp[[i]] %>%
  mutate_if( is.numeric, round, 3 )

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
