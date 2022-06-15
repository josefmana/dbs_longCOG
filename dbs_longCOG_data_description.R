# All analyses reported in the article (Mana et al., in review) ran in R version 4.2.0 (2022-04-22),
# on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.4.

# I used the following versions of packages employed: dplyr_1.0.9, DiagrammeR_1.0.9, DiagrammeRsvg_0.1,
# rsvg_2.3.1, missMDA_1.18, ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "tidyverse", # for rownames_to_column
  "DiagrammeR", # for flowcharts
  "DiagrammeRsvg", # for saving plots created in DiagrammeR
  "rsvg", # for saving plots created in DiagrammeR
  "ggplot2", # for plotting
  "patchwork" # for ggplots manipulations
)

# load required packages
# prints NULL if a package is already installed
sapply( c("models", "figures", "tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_classic(base_size = 25) )

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# read the data set and prepare subsets for individual analyses
d0 <- read.csv( "data/20220508_dbs_longCOG_data.csv" , sep = "," )
d1 <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data
d2 <- d1[ d1$ass_type == "pre" , ] # only pre-surgery assessments of included patients

# read a file with mapping variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")


# ----------- Fig 1 participant inclusion flowchart  -----------

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
  node [ fontname = Calibri, fontsize = 24, shape = box, style = rounded , width = 5 , margin= 0.5 , penwidth = 2 ];
  
  /// create a box for all patients, i.e., sum(t) = 200 patients
  all_pats [ label =<
  <b>200 consecutive<br/>PD patients </b><br/><br/>Local database 2000-2020<br/>General University Hospital<br/>in Prague
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
  node [ fixedsize = T, width = 6, height = 2, margin= 0.5 ];
  
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
grViz(f1) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_png("figures/Fig 1 inclusion-exclusion flowchart.png")


# ----------- Fig 2 distribution of assessments  -----------

# prepare a histogram of the distribution of assessments across time (Fig 2A)
f2 <- list()

# need to use the complete.cases command because three patients have duplicated rows
# due to more than one stimulation parameter
f2$hist <- d1[ complete.cases(d1$drs_tot) , ] %>% 
  ggplot( aes(x = time_y) ) +
  stat_bin( binwidth = .5, position = position_nudge(x = -.5*.5) ) + # creates bars
  stat_bin( binwidth = .5, geom = "text", aes(label = ..count..), vjust = -1.0,
            position = position_nudge(x = -.5*.5), size = 6 ) + # add numbers
  labs( x = "Time from STN-DBS surgery (years)", y = "Number of Assessments" ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0, 69), breaks = seq(0, 60, 10), labels = seq(0, 60, 10) ) +
  scale_x_continuous( limits = c(-2, 12), breaks = seq(-2, 12, 1), labels = seq(-2, 12, 1) )

# prepare a bin plot showing distribution of the number of assessments per patient (Fig 2B)
f2$bin <- table( d1[ complete.cases(d1$drs_tot) , ]$id ) %>%
  as.data.frame() %>%
  ggplot( aes(x = Freq) ) +
  geom_bar( width = .25 ) +
  geom_text( stat = "count", aes(label = ..count..), vjust = -1.0, size = 6 ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0, 65) ) +
  labs( x = "Number of Assessments per Patient", y = "Number of Patients" )

# arrange Fig 2A and Fig 2B for printing
f2$hist / f2$bin + plot_annotation( tag_levels = "A" )

# save as Fig 2
ggsave( "figures/Fig 2 distribution of assessments.png" , height = 2.5 * 6.12 , width = 1.5 * 11.6 , dpi = "retina" )


# ----------- Tab 1 stimulation parameters  -----------

# list all stimulation parameters
pars <- names(d1)[which(names(d1)=="current_r_mA"):which(names(d1)=="frequency_l_Hz")]

# prepare a data frame to be filled-in with stimulation parameters' summary
t1 <- data.frame( Md = rep( NA , length(pars) ), `Min-Max` = NA, M = NA, SD = NA, row.names = pars )

# fill-in statistics of all stimulation parameters
for ( i in pars ) {
  t1[ i , "Md" ] <- sprintf( "%.2f" , round( median( d1[[i]], na.rm = T ) , digits = 2 ) ) # median, rounded to hundredths
  t1[ i , "Min.Max" ] <- paste0( sprintf( "%.1f" , round( min( d1[[i]], na.rm = T ) , digits = 1 ) ), "-",
                                 sprintf( "%.1f" , round( max( d1[[i]], na.rm = T ) , digits = 1 ) ) ) # min-max, rounded to decimals
  t1[ i , "M" ] <- sprintf( "%.2f" , round( mean( d1[[i]], na.rm = T ) , digits = 2 ) ) # mean, rounded to hundredths
  t1[ i , "SD" ] <- sprintf( "%.2f" , round( sd( d1[[i]], na.rm = T ) , digits = 2 ) ) # SD, rounded to hundredths
}

# add empty rows and change names such that the table is publication-ready
t1 <- t1 %>%
  add_row( .before = 1 ) %>% # current
  add_row( .before = 4 ) %>% # voltage
  add_row( .before = 7 ) %>% # duration
  add_row( .before = 10 ) %>% # frequency
  # add a column with names
  mutate( Parameter = c("Current mode (N = 67, mA)", "right", "left",
                        "Voltage mode (N = 59, V)", "right", "left",
                        "Pulse duration (µs)", "right", "left",
                        "Frequency (Hz)", "right", "left"), .before = 1 )

# save as csv for import to word editor 
write.table( t1, file = "tables/Tab 1 stimulation parameters.csv", sep = ",", row.names = F, na = "" )


# ----------- Tab 2 baseline characteristic  -----------

# list all variables that will be included in Tab 2
vars <- names(d2)[which(names(d2)=="age_stim_y"):which(names(d2)=="fp_dr")]

# prepare a data frame to be filled-in with baseline characteristics
t2 <- data.frame( N = rep( NA, length(vars) ), Md = NA, `Min-Max` = NA, M = NA, SD = NA, row.names = vars )

# fill-in statistics of all variables but sex (which is nominal)
for ( i in vars[-3] ) {
  t2[ i , "N" ] <- sum( !is.na(d2[[i]]) ) # number of data points
  t2[ i , "Md" ] <- sprintf( "%.0f" , round( median(d2[[i]], na.rm = T ) , digits = 0 ) ) # median, rounded to integers
  t2[ i , "Min.Max" ] <- paste0( sprintf( "%.0f" , round( min(d2[[i]], na.rm = T ) , digits = 0 ) ), "-",
                                 sprintf( "%.0f" , round( max(d2[[i]], na.rm = T ) , digits = 0 ) ) ) # min-max, rounded to integers
  t2[ i , "M" ] <- sprintf( "%.2f" , round( mean(d2[[i]], na.rm = T ) , digits = 2 ) ) # mean, rounded to hundredths
  t2[ i , "SD" ] <- sprintf( "%.2f" , round( sd(d2[[i]], na.rm = T ) , digits = 2 ) ) # SD, rounded to hundredths
}

# add frequency of males to the table sex row
t2[ "sex" , c("N","Min.Max")] <- c(
  sum( !is.na(d2$sex) ), # number of entries
  # add frequency (percentage) to Min.Max
  paste0( table(d2$sex)["male"], " (", # frequency
          sprintf( "%.0f" , round( 100 * ( table(d2$sex)[ "male" ] / sum( table(d2$sex) ) ), 0 ) ), " %)" ) ) # percentage, rounded to integers

# prepare row names by making them a column and the changing names to a publication-ready format
t2 <- t2 %>% rownames_to_column( var = "Characteristic" )

# loop through all variables in Tab 2 and change their names accordingly
for ( i in 1:nrow(t2) ) t2[i, "Characteristic"] <- var_nms[ t2[i,]$Characteristic, ]

# save as csv for import to word editor 
write.table( x = t2, file = "tables/Tab 2 baseline characteristics.csv", sep = ",", row.names = F, na = "" )
