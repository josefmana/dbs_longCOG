# All analyses reported in the article (Mana et al., in review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 19043).

# I used the following versions of packages employed: dplyr_1.0.7, tidyverse_1.3.0,
# DiagrammeR_1.0.9, brms_2.16.3, psych_2.2.3, tidybayes_2.3.1, ggplot2_3.3.3 and patchwork_1.1.1

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
theme_set( bayesplot::theme_default(base_size = 25) )

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

# read the dataset
# In this file, my strategy is to label the sequential version
# of dataset as d# (where # goes for a number of the data version)
d0 <- read.csv( "data/20220423_dbs_longCOG_data.csv" , sep = "," )


# ----------- participant inclusion flowchart  -----------

# check that when selecting only rows containing pre-surgery assessment,
# there is no patient duplicated
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


# ----------- the total effect of time  -----------


# ----------- the direct effect of time  -----------


# ----------- the pre-surgery cognitive profile  -----------


# ----------- participant inclusion flowchart  -----------