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
t <- as.data.frame( table(d0[ d0$ass_type == "pre" , ]$why_excluded) )

# extract the numbers for each part of the flowchart
nums <- list(
  # total number of patients
  total = sum( t$Freq ),
  # patients that didn't receive STN-DBS treatment
  no_stn = list(
    total = sum(
      t$Freq[ t$Var1 %in% c("GPi","VIM","duodopa","rejected","suspended")]
    ),
    gpi = t$Freq[ t$Var1 == "GPi" ],
    vim = t$Freq[ t$Var1 == "VIM" ],
    duodopa = t$Freq[ t$Var1 == "duodopa" ],
    rejected = t$Freq[ t$Var1 == "rejected" ],
    suspended = t$Freq[ t$Var1 == "suspended" ]
  ),
  # patients treated with STN-DBS
  stn_dbs = list(
    total = sum(
      t$Freq[ t$Var1 %in% c("included","left_stn_only","right_stn_only","any_psycho_missing","only_drs_before_surgery","pre_psycho_missing",
                            "post_psycho_missing","language") ]
    ),
    included = t$Freq[ t$Var1 == "included" ],
    unilateral = sum( t$Freq[ t$Var1 %in% paste0( c("left","right") , "_stn_only" ) ] ),
    preop_miss = sum( t$Freq[ t$Var1 %in% c("any_psycho_missing","only_drs_before_surgery","pre_psycho_missing") ] ),
    follow_miss = t$Freq[ t$Var1 == "post_psycho_missing" ],
    language = t$Freq[ t$Var1 == "language" ]
  )
)

# create a flowchart in DiagrammeR
grViz(
  " digraph {
  
  # node definitions with substituted label text
  node [fontname = Calibri, shape = box, style = rounded , width = 2.5]
  
  1 [label =
      <
      <br/><b>200 consecutive PD patients </b><br/><br/>
      Local database 2000-2020<br/>
      General University Hospital<br/>
      in Prague<br/>
      >
    ]
  
  2 [label =
      <
      <br/><b>173 patients </b><br/><br/>
      implanted with STN-DBS<br/>
      >
    ]
  
  3 [label =
      <
      <br/><b>126 patients </b><br/><br/>
      followed-up longitudinally<br/>
      >
    ]
  
  # nodes for exluded patients
  node [fontname = Calibri, shape = box, style = rounded, width = 3.2]
  
  m1 [label = 
        <
        <br align = 'left'/>
        <b>27 patients excluded due to:</b><br align = 'left'/><br align = 'left'/>
            12 GPi-DBS<br align = 'left'/> 
            4 VIM-DBS<br align = 'left'/> 
            4 duodopa<br align = 'left'/> 
            5 rejected<br align = 'left'/> 
            2 suspended<br align = 'left'/>
        >
      ]

  m2 [label =
        <
        <br align = 'left'/>
        <b>47 patients excluded due to:</b><br align = 'left'/><br align = 'left'/>
            24 pre-surgery data missing<br align = 'left'/>
            18 follow-up data missing<br align = 'left'/>
            3 unilateral STN-DBS<br align = 'left'/>
            2 not speaking Czech<br align = 'left'/>
        >
     ]
  
  # dummy nodes
  node [ shape = rectangle, width = 0, height = 0, label = '' , fill = black]
  
  # edge definitions
  p1 -> 2; p2 -> 3;
  {rank = same; p1 -> m1}
  {rank = same; p2 -> m2}
  
  edge [dir = none]
  1 -> p1; 2 -> p2;

  # seperate nodes in the same rank (in inches)
  nodesep = 1
  
  }
  "
)

# ----------- directed acyclic graph -----------


# ----------- sample description  -----------


# ----------- the total effect of time  -----------


# ----------- the direct effect of time  -----------


# ----------- the pre-surgery cognitive profile  -----------


# ----------- participant inclusion flowchart  -----------