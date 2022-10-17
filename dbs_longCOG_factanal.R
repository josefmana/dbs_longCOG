# Ran in R version 4.2.0 (2022-04-22), on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.6.

# I used the following versions of packages employed: dplyr_1.0.9, tidyverse_1.3.1, missMDA_1.18, psych_2.2.5,
# ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", "tidyverse", # for data wrangling
  "missMDA", "psych", # for factor analyses 
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

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_classic(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# read the data set
d <- read.csv( "data/20220508_dbs_longCOG_data.csv", sep = "," ) %>% filter( included == 1 & ass_type == "pre" )

# read a file containing mapping of variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")


# ----------- pre-processing  -----------

# for EFA keep only id and cognitive tests in d2
d <- d[ , c( 2, which(names(d) == "tmt_a"):which(names(d) == "fp_dr"), which(names(d) %in% paste0("staix",1:2)) ) ]

# change names such that they get correct label in post-processing
# "cs_" stands for "centered & scaled" because the variables will be standardized before analyses
colnames(d)[-1] <- paste0( "cs_", colnames(d)[-1] )

# log-transform reaction times before the analysis
for ( i in c( paste0("cs_tmt_", c("a","b")), paste0("cs_pst_", c("d","w","c")) ) ) d[[i]] <- log( d[[i]] )

# standardize (center and scale) all test scores before analysis
for ( i in colnames(d)[-1] ) d[,i] <- as.vector( scale(d[,i], center = T, scale = T) )


# ----------- imputation  -----------

# find out the optimal number of components for multiple imputation
nb <- estim_ncpPCA( d[,-1] , ncp.min = 0, ncp.max = 10 , nbsim = imp )

# impute via PCA-based multiple imputation (n = 100 imputations)
set.seed(s) # set seed for reproducibility
d.imp <- MIPCA( d[ ,- 1] , ncp = nb$ncp , nboot = imp )


# ----------- model fitting & selection -----------

# fit EFAs with 3:8 factors to each imputed data set
efa <- lapply(
  # loop through 100 imputed data sets
  1:imp, function(i)
    # loop through three to eight latent factors for each imputation
    lapply( 3:8, function(j) fa( d.imp$res.MI[[i]], # one-by-one use each imputed data 
                                 nfactors = j, # fit 3-8 factor solutions
                                 rotate = "varimax", # rotate varimax to enforce orthogonality and for interpretation purposes
                                 scores = "regression" # compute regression scores for each patient
    )
  )
)

# one-by-one inspect all six-factor and seven-factor solutions' loading matrices
# create a convenience function so that I don't go crazy immediately
print_load <- function( i, c = .4, f = 7 ) print( efa[[i]][[f-2]]$loadings, cutoff = c, sort = T )

# choosing 7-factor solution due to good performance indexes,
# and theoretically sound loading patterns across imputed data sets
nf = 7


# ----------- post-processing  -----------

# prepare an array for labels of the seven-factor solution factors
# create an empty 2 (names/signs) x 7 (factors) x 100 (imputations) array
doms_sum <- array( data = NA, dim = c(2, nf, imp), dimnames = list( c("nms","sgn"), paste0("F", 1:nf), 1:imp ) )

# read the table with seven-factor labels
doms_sum["nms", , ] <- t( read.csv( "data/dbs_longCOG_efa_labels.csv" , sep = "," , row.names = 1, header = T) )

# fill-in signs of each factor in each imputation to know which scores should be reversed
doms_sum["sgn", , ] <- apply( doms_sum["nms", , ] , 2 , function(x) startsWith( x , "-") ) %>%
  t() %>% as.data.frame() %>% mutate( across( everything() , ~ ifelse( . == T , -1 , 1 ) ) ) %>% t()

# get rid of the minus sign in labels table
doms_sum["nms", , ] <- doms_sum["nms", , ] %>%
  t() %>% as.data.frame() %>% mutate( across( everything() , ~ gsub( "-" , "" , . ) ) ) %>% t()

# list all the domains
doms <- c("exec_fun", # loaded on primarily by PST, the first factor in 82% data sets
          "epis_mem", # loaded on primarily by RAVLT, the second factor in 79% data sets
          "verb_wm", # loaded on primarily by DS, the third factor in 62% data sets
          "visp_mem", # loaded on primarily by FP, the fourth factor in 45% data sets
          "set_shift", # loaded on primarily by TMT and RAVLT-B, the fifth factor in 28% data sets
          "anxiety", # loaded on primarily by STAI, the sixth factor in 60% data sets
          "visp_wm" # loaded on primarily by SS, the seventh factor in 49% data sets
)

# switch signs where appropriate in EFA loadings and scores, and rename and sort columns
for ( i in 1:imp ) {
  for ( j in c("loadings","scores","Vaccounted") ) {
    
    # multiply by a diagonal matrix of 1 and -1
    if ( j %in% c("loadings","scores") ) efa[[i]][[nf-2]][[j]] <- efa[[i]][[nf-2]][[j]] %*% diag(doms_sum["sgn", , i ] )
    
    # rename the columns
    colnames( efa[[i]][[5]][[j]] ) <- doms_sum["nms", , i ]
    
    # reorder the columns such that they are in the same order for each imputation
    efa[[i]][[nf-2]][[j]] <- efa[[i]][[nf-2]][[j]][, doms]
    
  }
}


# ----------- saving the outcomes  -----------

# save all EFA models
saveRDS( object = efa, file = "models/factanal.rds" )

# prepare a folder for imputed data sets
if ( !dir.exists("data/imputed") ) dir.create( "data/imputed" )

# loop through imputed data sets and save each to its own .csv
for ( i in 1:imp ) write.table( x = cbind.data.frame( id = d$id, d.imp$res.MI[[i]] ), # the data
                                file = paste0("data/imputed/imputed_df_",i,".csv"), # the file
                                sep = ",", row.names = F
)


# ----------- extracting performance indexes  -----------

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
      efa[[i]][[j]]$TLI, # Tucker-Lewis Index, TLI > 0.9 is considered good
      efa[[i]][[j]]$RMSEA[1], # Root-Mean Square Error of Approximation, RMSEA < 0.08 is considered good
      efa[[i]][[j]]$RMSEA[2], # 90% CI RMSEA, lower boundary
      efa[[i]][[j]]$RMSEA[3], # 90% CI RMSEA, upper boundary (ideally should be less than 0.08)
      max( efa[[i]][[j]]$Vaccounted["Cumulative Var",] ) # total variance "accounted for" by all included factors (i.e., the highest one from "Cumulative Var")
    )
  }
}


# ----------- tab s1 summary of performance index of factor analyses  -----------

# summarize the fat table
t.s1 <- data.frame( Model = paste0( 3:8, "-factor"), TLI = NA, RMSEA = NA, RMSEA_90_CI_upp = NA, var_account = NA )

# fill-in with models' summaries
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
    as.data.frame() %>%
    mutate( across( everything() , ~ replace( . , . > .08 , NA) ) ),
  # TLI higher than 0.9
  `TLI > 0.90 (%)` = t( fat[ , "TLI" , ] ) %>%
    as.data.frame() %>%
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


# ----------- fig s2 factor analyses performance indexes  -----------

# set-up a list to contain Fig S2 component figures
f.s2 <- list()

# loop through TLI and upper 90% CI RMSEA
for ( i in c("TLI","RMSEA_90_CI_upp") ) {
  f.s2[[i]] <- fat[ , i , ] %>%
    # re-format the table for plotting
    t() %>% as.data.frame() %>%
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
ggsave( "figures/Fig S2 factor analysis performance indexes.jpg", dpi = 600 )


# ----------- tab 2 factor loadings  -----------

# prepare an array for loading matrices of each imputed EFA
loads <- array(
  # create an empty 25 ( 23 tests + 2 variance accounted) x 7 (factors) x 100 (imputations) array
  data = NA, dim = c(25, nf, imp),
  dimnames = list( c( rownames(efa[[1]][[nf-2]]$loadings) , "Proportion Var", "Cumulative Var" ), doms, 1:imp )
)

# fill-in all loadings from efa objects prepared above
for( i in 1:imp ) loads[ , , i ] <- efa[[i]][[nf-2]]$loadings %>%
  as.data.frame() %>%
  bind_rows( efa[[i]][[nf-2]]$Vaccounted[ "Proportion Var", ] ) %>%
  bind_rows(
    apply( efa[[i]][[nf-2]]$Vaccounted[ "Proportion Var", ] %>% t(), 1 , cumsum ) %>% t() %>% as.data.frame()
  ) %>% as.matrix()

# write a summary of loadings across all imputed data sets
t2 <- matrix(
  data = NA, nrow = nrow(loads[ , , 1]), ncol = ncol(loads[ , , 1]),
  dimnames = list( rownames(loads[ , , 1 ]) , colnames(loads[ , , 1]) )
) %>% as.data.frame()

# fill-in averages and SDs across all imputations 
for ( i in rownames(t2) ) {
  t2[ i , ] <- paste0(
    sprintf( "%.2f" , round( loads[i, , ] %>% t() %>% colMeans() , 2 ) ), " (", # mean
    sprintf( "%.2f" , round( loads[i, , ] %>% t() %>% apply( . , 2 , sd ), 2) ), ")" # SD
  )
}

# make rownames to a column
t2 <- t2 %>% rownames_to_column( var = "Test" )

# rename tests (rows) and factors (columns) such that they are publication-ready
for ( i in 1:(nrow(t2)-2) ) t2[i,"Test"] <- var_nms[ t2[i,]$Test,  ]
for ( i in 2:ncol(t2) ) colnames(t2)[i] <- var_nms[ colnames(t2)[i],  ]

# save as csv
write.table( t2, file = "tables/Tab 2 factor loadings.csv", sep = ",", row.names = F, quote = F )


# ----------- session info -----------

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/factor_analysis.txt" )
