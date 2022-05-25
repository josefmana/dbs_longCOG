# All analyses reported in the article (Mana et al., in review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 22000).

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
imp = 20
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
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")


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
d2.imp <- MIPCA( d2[ ,- 1] , ncp = nb$ncp , nboot = imp )

# calculate parallel test for each imputed data set 
set.seed(s) # set seed for reproducibility
p.test <- lapply( 1:imp , function(i) fa.parallel( d2.imp$res.MI[[i]] , plot = F ) )

# look at the results of parallel tests
table( sapply( 1:imp , function(i) p.test[[i]]$nfact ) )

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

# calculate the percentages
for ( i in names(fat_perc) ) fat_perc[[i]] <- sapply(
  names(fat_perc[[i]]), function(j) 100 * sum( complete.cases( fat_perc[[i]][[j]] ) ) / imp
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
    labs( y = "Density", x = case_when( i == "TLI" ~ "TLI", i == "RMSEA_90_CI_upp" ~ "RMSEA (upper 90% CI)" ) ) +
    scale_x_continuous( limits = case_when( i == "TLI" ~ c(0.6,1.05), i == "RMSEA_90_CI_upp" ~ c(0.03, 0.12) ),
                        breaks = case_when( i == "TLI" ~ seq(.6,1,.1), i == "RMSEA_90_CI_upp" ~ seq(.04,.12,.02) ),
                        labels = case_when(
                          i == "TLI" ~ sprintf( "%.1f", round( seq(.6, 1, .1), 1) ),
                          i == "RMSEA_90_CI_upp" ~ sprintf( "%.2f", round( seq(.04,.12,.02), 2) )
                        ))
}

# arrange Fig S1 for printing
f.s1$TLI / f.s1$RMSEA_90_CI_upp + plot_layout( guides = "collect" ) + plot_annotation( tag_levels = "A" )

# save as Fig S1
ggsave(
  "figures/Fig S1 factor analysis performance indexes.png",
  height = 1.75 * 6.07, width = 1.75 * 11.5, dpi = "retina"
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
  # create an empty 25 ( 23 tests + 2 variance accounted) x 7 (factors) x 20 (imputations) array
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


# write a summary of loadings across all 20 imputed data sets
t3 <- matrix(
  data = NA, nrow = nrow(loads[ , , 1]), ncol = ncol(loads[ , , 1]),
  dimnames = list( rownames(loads[ , , 1 ]) , colnames(loads[ , , 1]) )
) %>% as.data.frame()

# fill-in averages and SDs across all 20 imputations 
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
sapply (
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
      cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ) # right censoring for DRS == 144
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
rm( list = ls()[ !( ls() %in% c( "cbPal", "d1", "d3", "doms", "imp", "s", "scl", "var_nms" ) ) ] )
gc() # garbage collection to ease RAM a bit if needed

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .95 # adapt_delta parameter


# ----------- baseline model with correlated random-effects  -----------

# set-up the linear model
f0.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time | id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf

# set-up priors
p0 <- c(
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
  prior( lkj(2), class = cor ),
  # other distributional parameters
  prior( exponential(1), class = sigma ),
  prior( gamma(2, 0.1), class = nu )
)

# fit the model with Student response function
m <- list(
  m0_base = brm_multiple(
    formula = f0.drs, family = student(), prior = p0, data = d3,
    sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
    control = list( adapt_delta = ad ), file = "models/m0_base.rds"
  )
)


# ----------- primary model with uncorrelated random-effects  -----------

# set-up the linear model
f1.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time || id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf

# set-up priors
p1 <- p0[ -which(p0$class == "cor"), ]

# fit the model with Student response function
m$m1_nocov <- brm_multiple(
  formula = f1.drs, family = student(), prior = p1, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m1_nocov.rds"
)


# ----------- sensitivity check with flat priors  -----------

# the same linear model as for m1
f2.drs <- f1.drs

# will be using brms default priors
p2 <- NULL

# fit the model with Student response function
m$m2_flat_priors <- brm_multiple(
  formula = f2.drs, family = student(), prior = p2, data = d3,
  sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu,
  control = list( adapt_delta = ad ), file = "models/m2_flat_priors.rds"
)


# ----------- model with covariates  -----------

# set contrast for sex as a factor
for( i in 1:imp ) contrasts(d3[[i]]$sex) <- -contr.sum(2)/2 # female = -0.5, male = 0.5

# set-up the linear model for outcome
f3.drs <- paste0(
  "drs | cens(cens_drs) ~ 1 + age + mi(bdi) + mi(led) + ", # outcome and intercept
  paste( "time", doms, sep = " * " , collapse = " + " ), # population-level effects/fixed-effects
  " + (1 + time || id)"  # varying-effects (patient-level)/random-effects
) %>% as.formula %>% bf + student()

# set-up linear models for covariates with missing values
f3.bdi <- bf( bdi | mi() ~ 1 + time + sex + age + mi(led) + (1 + time || id) ) + gaussian()
f3.led <- bf( led | mi() ~ t2(time) + (1 | id) ) + gaussian()

# set-up priors
p3 <- c(
  # DRS-2
  prior( normal(0.3, .1), class = Intercept, resp = drs ),
  prior( normal(-.2, .1), class = b, coef = time, resp = drs ),
  prior( normal(0, .1), class = b, coef = age, resp = drs ),
  prior( normal(0, .1), class = b, coef = mibdi, resp = drs ),
  prior( normal(0, .1), class = b, coef = miled, resp = drs ),
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
  prior( normal(0, .1), class = sd, coef = Intercept, group = id, resp = drs ),
  prior( normal(0, .1), class = sd, coef = time, group = id, resp = drs ),
  prior( exponential(1), class = sigma, resp = drs ),
  prior( gamma(2, 0.1), class = nu, resp = drs ),
  # BDI-II
  prior( normal(.6, .5), class = Intercept, resp = bdi ),
  prior( normal(0, .5), class = b, coef = time, resp = bdi ),
  prior( normal(0, .5), class = b, coef = sex1, resp = bdi ),
  prior( normal(0, .5), class = b, coef = miled, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
  prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
  prior( exponential(1), class = sigma, resp = bdi ),
  # LEDD
  prior( normal(0, 100), class = Intercept, resp = led ),
  prior( normal(0, 100), class = b, coef = t2time_1, resp = led ),
  prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = led ),
  prior( exponential(1), class = sigma, resp = led )
)

# fit the model as defined above
m$m3_wcov <- brm_multiple(
  formula = f3.drs + f3.bdi + f3.led, prior = p3, data = d3,
  sample_prior = F, seed = s, chains = ch, iter = it, warmup = wu, # not saving priors to spare some memory
  control = list( adapt_delta = .99 ), file = "models/m3_wcov.rds"
)

# ----------- post-processing: extract MCMC draws  -----------

# extract posteriors from all models
draws <- list( m0_base = posterior::as_draws_df(m$m0_base) )

# loop through the rest of the models to complete the draws list
for ( i in names(m)[-1] ) draws[[i]] <- posterior::as_draws_df( m[[i]] )

# save them as rds for sharing
saveRDS( draws , "models/dbs_longCOG_mcmc_draws.rds" )


# ----------- post-processing: compute PSIS-LOO  -----------

# compute PSIS-LOO for each imputed data set for each model
l <- list()

# loop through each model and data set
# doesn't work for the m4_wcov model in brms currently (likely because of the online imputation of covariates)
for ( i in names(m)[-4] ) {
  for ( j in 1:imp ) l[[i]][[j]] <- loo( m[[i]], newdata = d3[[j]] )
}

# save all PSIS-LOO as rds
saveRDS( l, "models/dbs_longCOG_psis-loo.rds" )


# ----------- post-processing: prepare posterior predictions  -----------

# prepare posterior predictions for Processing Speed and Episodic Memory for Fig. 4 in a multistep procedure
# first clear the environment of models we will not need anymore
m <- m$m1_nocov # using only the primary model to be reported in the main text
gc()

# write down a list of predictors effect of which we're going to visualize
prds <- c( "proc_spd" , "epis_mem" )

# prepare a raw data set for visualizations
d4 <- d1 %>%
  select( id, time_y , drs_tot ) %>%
  filter( complete.cases(drs_tot) ) %>%
  rename( "drs" = "drs_tot" ) %>%
  # add median of imputed Processing Speed and Episodic Memory
  mutate(
    proc_spd = sapply( 1:imp , function(i) d3[[i]]$proc_spd ) %>% apply( . , 1 , median ),
    epis_mem = sapply( 1:imp , function(i) d3[[i]]$epis_mem ) %>% apply( . , 1 , median )
  )

# get quantile groups for each patients according to Processing Speed and Episodic Memory
d5 <- d4[ d4$time_y < 0 , ] %>%
  mutate(
    # re-coding such that 1 means the lowermost quantile and 4 means the uppermost quantile
    proc_spd_quant = -ntile( proc_spd, 4)+5,
    epis_mem_quant = -ntile( epis_mem, 4)+5
  )

# add these groups to the longitudinal data set (d4)
d4 <- d4 %>% left_join( d5 %>% select( id, proc_spd_quant, epis_mem_quant) , by = "id" )

# prepare a prediction dummy data set
d_seq <- list()

# write down how many predictions (time-points) per predictor group/quantile to calculate
n_seq = 25

# loop through predictors
for ( i in prds ) d_seq[[i]] <- expand.grid(
  as.vector( by( d4[[i]], d4[[ paste0(i, "_quant") ]] , median ) ),
  seq( from = -2, to = 10, length.out = n_seq )
) %>% `colnames<-` ( c(i, "time_y") ) %>%
  # add all variables needed
  mutate(
    time = time_y + scl$Md$time,
    proc_spd = if( i == "proc_spd") proc_spd else 0,
    epis_mem = if( i == "epis_mem") epis_mem else 0,
    grp = rep( 1:4, n_seq ), # in the expand.grid above, need to write predictor first, time second, otherwise this is incorrect
    verb_wm = 0, visp_mem = 0, set_shift = 0, anxiety = 0, visp_wm = 0, id = NA
  )

# add predictions to d_seq
for ( i in prds ) d_seq[[i]] <- d_seq[[i]] %>%
  add_fitted_draws( m , re_formula = NA ) %>%
  mutate(.value = scl$M$drs + scl$SD$drs * .value) %>%
  median_hdi(.width = .95)

# prepare variable names for the plot
nms <- list( proc_spd = "Processing speed" , epis_mem = "Episodic memory" )


# ----------- post-processing: drawing the figures  -----------

# first prepare list with posteriors from each model
# make a list of variable names for all fixed-effects
pars <- list(
  intercept = names(draws$m1_nocov)[1], # global intercept
  base = names(draws$m1_nocov)[3:9], # baseline correlates
  time = names(draws$m1_nocov)[c(2,10:16)] # time-dependent effects
)

# prepare posterior lists
post <- list()

# loop through all models to fill-in posteriors in question
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
      grepl("m0_base", Model) ~ "Correlated RE",
      grepl("m1_nocov", Model) ~ "Uncorrelated RE",
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

# change character columns to ordered factors such that ggplot plots them in a correct order
# need to reverse the order by rev(.) because ggplot counts from bottom-up on the y-axis
post$Model <- factor( post$Model , levels = unique(post$Model)  , ordered = T )
post$Parameter <- factor( post$Parameter , levels = rev( unlist(pars, use.names = F) ) , ordered = T )
post$Group <- factor( post$Group , levels = rev( unique(post$Group) ) , ordered = T )

# plot it one by one per group
f.s2 <- list()

# loop through all three parameter groups
for ( i in rev( levels(post$Group) ) ) f.s2[[i]] <- post[ post$Group == i , ] %>%
  mutate(
    title = case_when( i == "intercept" ~ "Global intercept",
                       i == "base" ~ "Baseline correlates",
                       i == "time" ~ "Time-dependent effects")
  ) %>%
  ggplot( aes(y = Parameter, x = `DRS-2`, fill = Model, color = Model) ) +
  stat_halfeye( geom = "slab", slab_linetype = "solid" , slab_size = 1,  slab_alpha = .5 ) +
  geom_vline( xintercept = case_when( i == "intercept" ~ 139, i != "intercept" ~ 0 ), linetype = "dashed" ) +
  scale_color_manual( values = cbPal[2:5]) +
  scale_fill_manual( values = cbPal[2:5]) +
  scale_x_continuous( limits = case_when( i == "intercept" ~ c(138, 143), i != "intercept" ~ c(-3,2) ) ) +
  scale_y_discrete( name = NULL, labels = rev( var_nms[ pars[[i]] , ] ) ) +
  facet_wrap( . ~ title) +
  theme(
    legend.position = "bottom" , plot.title = element_text(hjust = 0.5) , axis.title.x = element_text(size = 18)
  )

# put them together
f.s2$intercept / f.s2$base / f.s2$time  +
  plot_layout( heights = c(1, 7, 8) , guides = "collect" ) & theme( legend.position = "bottom" )

# save as Fig S2B
ggsave(
  "figures/Fig S2B posteriors across-models.png",
  height = 3 * 8.53, width = 1.75 * 9.05, dpi = "retina"
)
