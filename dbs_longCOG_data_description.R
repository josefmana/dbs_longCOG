# Ran in R version 4.2.0 (2022-04-22), on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.6.

# I used the following versions of packages employed: dplyr_1.0.9, tidyverse_1.3.1, DiagrammeR_1.0.9,
# DiagrammeRsvg_0.1, rsvg_2.3.1, ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages
pkgs <- c(
  "dplyr", "tidyverse", # for data wrangling
  "DiagrammeR", "DiagrammeRsvg", "rsvg", # for flowcharts
  "ggplot2", "patchwork" # for general plotting
)

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_classic(base_size = 14) )

# read the data set and prepare subsets for individual analyses
d0 <- read.csv( "data/20220508_dbs_longCOG_data.csv" , sep = "," )
d1 <- d0[ d0$included == 1 , ] # only STN-DBS treated patients with pre- and post-surgery data
d2 <- d1[ d1$ass_type == "pre" , ] # only pre-surgery assessments of included patients

# read a file containing mapping of variables' names used in the script to variables' names for the manuscript
var_nms <- read.csv( "data/var_nms.csv" , sep = ";" , row.names = 1 , encoding = "UTF-8")


# ----------- Fig 1 participant inclusion flowchart  -----------

# check that when selecting only rows containing pre-surgery assessment there ain't no patient duplicated
isTRUE( all.equal( d0[ d0$ass_type == "pre" , ]$id, unique( d0[ d0$ass_type == "pre" , ]$id ) ) ) # TRUE

# print a table summarizing reasons for excluding patients
table(d0[ d0$ass_type == "pre" , ]$why_excluded) %>% print()

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
grViz(f1) %>% export_svg %>% charToRaw %>% rsvg_png("figures/Fig 1 inclusion-exclusion flowchart.png")


# ----------- Fig S1 distribution of assessments  -----------

# prepare a histogram of the distribution of assessments across time (Fig S1a)
f.s1 <- list()

# need to use the complete.cases command because three patients have duplicated rows
# due to more than one stimulation parameter
f.s1$hist <- d1[ complete.cases(d1$drs_tot) , ] %>% 
  ggplot( aes(x = time_y) ) +
  stat_bin( binwidth = .5, position = position_nudge(x = -.5*.5) ) + # creates bars
  stat_bin( binwidth = .5, geom = "text", aes(label = ..count..), vjust = -1.0,
            position = position_nudge(x = -.5*.5), size = 4 ) + # add numbers
  labs( x = "Time from STN-DBS surgery (years)", y = "Number of Assessments" ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0, 74), breaks = seq(0, 70, 10), labels = seq(0, 70, 10) ) +
  scale_x_continuous( limits = c(-2, 12), breaks = seq(-2, 12, 1), labels = seq(-2, 12, 1) )

# prepare a bin plot showing distribution of the number of assessments per patient (Fig S1b)
f.s1$bin <- table( d1[ complete.cases(d1$drs_tot) , ]$id ) %>%
  as.data.frame() %>%
  ggplot( aes(x = Freq) ) +
  geom_bar( width = .25 ) +
  geom_text( stat = "count", aes(label = ..count..), vjust = -1.0, size = 4 ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0, 65) ) +
  labs( x = "Number of Assessments per Patient", y = "Number of Patients" )

# arrange for printing
f.s1$hist / f.s1$bin + plot_annotation( tag_levels = "a" ) & theme( plot.tag = element_text(face = "bold") )

# save as Fig S1
ggsave( "figures/Fig S1 distribution of assessments.jpg", dpi = 600 )


# ----------- Tab 1 & Tab 2 sample characteristics -----------

# list all variables for Tab 1 (clinics and stimulation parameters) and Tab 2 (neuropsychology)
nms <- list(
  dems = names(d2)[which(names(d2)=="age_stim_y"):which(names(d2)=="mds_updrs_iii_med_off")],
  pars = names(d1)[which(names(d1)=="current_r_mA"):which(names(d1)=="frequency_l_Hz")],
  tests = names(d2)[which(names(d2)=="drs_tot"):which(names(d2)=="fp_dr")]
)

# prepare a data frame to be filled-in with baseline characteristics and stimulation parameters
t <- lapply(
  names(nms), function(i)
    data.frame( N = rep(NA,length(nms[[i]])), Md = NA, `Min-Max` = NA, M = NA, SD = NA, row.names = nms[[i]] )
) %>% `names<-`( names(nms) )

# fill-in statistics of all baseline variables but sex (which is nominal)
for ( i in names(nms)[c(1,3)] ) {
  for ( j in nms[[i]] ) {
    
    # skip sex
    if (j == "sex" ) next
    
    # fill-in the rest
    t[[i]][ j , ] <- c(
      sum( !is.na(d2[[j]]) ), # number of data points
      sprintf( "%.0f" , round( median(d2[[j]], na.rm = T ) , digits = 0 ) ), # median, rounded to integers
      paste0( sprintf( "%.0f" , round( min(d2[[j]], na.rm = T ) , digits = 0 ) ), "-",
              sprintf( "%.0f" , round( max(d2[[j]], na.rm = T ) , digits = 0 ) ) ), # min-max, rounded to integers
      sprintf( "%.2f" , round( mean(d2[[j]], na.rm = T ) , digits = 2 ) ), # mean, rounded to hundredths
      sprintf( "%.2f" , round( sd(d2[[j]], na.rm = T ) , digits = 2 ) ) # SD, rounded to hundredths
    )
    
  }
}

# add frequency of males to the table sex row
t$dems[ "sex" , c("N","Min.Max")] <- c(
  sum( !is.na(d2$sex) ), # number of entries
  paste0(
    table(d2$sex)["male"], " (", # frequency
    sprintf( "%.0f" , round( 100 * ( table(d2$sex)[ "male" ] / sum( table(d2$sex) ) ), 0 ) ), " %)" # percentage, rounded to integers
  )
)

# fill-in statistics of stimulation parameters (all columns but N which I will fill-in by hand because
# the data are longitudinal and the N column is supposed to represent number of patients, not events)
for ( i in nms$pars ) t$pars[ i , 2:ncol(t$pars) ] <- c(
  sprintf( "%.1f" , round( median( d1[[i]], na.rm = T ) , digits = 1 ) ), # median, rounded to decimals
  paste0( sprintf( "%.1f" , round( min( d1[[i]], na.rm = T ) , digits = 1 ) ), "-",
          sprintf( "%.1f" , round( max( d1[[i]], na.rm = T ) , digits = 1 ) ) ), # min-max, rounded to decimals
  sprintf( "%.2f" , round( mean( d1[[i]], na.rm = T ) , digits = 2 ) ), # mean, rounded to hundredths
  sprintf( "%.2f" , round( sd( d1[[i]], na.rm = T ) , digits = 2 ) ) # SD, rounded to hundredths
)

# fill-in number of patients in each row of stimulation parameters table
t$pars$N <- c( rep(67,2) , rep(59,2) , rep(nrow(d2),4) )

# prepare row names by making them into a column
for ( i in names(t) ) t[[i]] <- t[[i]] %>% rownames_to_column( var = "Characteristic" )

# loop through all variables the pre-tables and change their names accordingly
for ( i in names(t) ) {
  for ( j in 1:nrow(t[[i]]) ) {
    # add indentions in Tab 1 (clinics and parameters) but not Tab 2 (neuropsychology)
    if (i != "tests" ) t[[i]][j, "Characteristic"] <- paste0( "\t", var_nms[ t[[i]][j,]$Characteristic, ] )
    else t[[i]][j, "Characteristic"] <- var_nms[ t[[i]][j,]$Characteristic, ]
  }
}

# prepare Tab 1 and Tab 2 proper
t1 <- do.call( rbind.data.frame, t[1:2] )
t2 <- t$tests %>% rename( "Test" = "Characteristic" )

# add header before each subsection of Tab 1
t1 <- t1 %>%
  add_row( Characteristic = "Baseline characteristics" , .before = 1 ) %>%
  add_row( Characteristic = "Stimulation parameters" , .before = length(nms$dems)+2 )

# save as csv for import to word editor
write.table( t1, file = "tables/Tab 1 clinical characteristics.csv", sep = ",", row.names = F, na = "", quote = F )
write.table( t2, file = "tables/Tab 2 baseline neuropsychology.csv", sep = ",", row.names = F, na = "", quote = F )


# ----------- descriptive stats for in-text reporting -----------

# duration of follow-up
mean( d1[ complete.cases(d1$drs_tot) & d1$ass_type != "pre", ]$time_y ) # 3.54
sd( d1[ complete.cases(d1$drs_tot) & d1$ass_type != "pre", ]$time_y ) # 2.32
median( d1[ complete.cases(d1$drs_tot) & d1$ass_type != "pre", ]$time_y ) # 3.07
min( d1[ complete.cases(d1$drs_tot) & d1$ass_type != "pre", ]$time_y ) # 0.72
max( d1[ complete.cases(d1$drs_tot) & d1$ass_type != "pre", ]$time_y ) # 11.38

# number of assessments per patient
median( table(d1[ complete.cases(d1$drs_tot) , ]$id) ) # 3
min( table(d1[ complete.cases(d1$drs_tot) , ]$id) ) # 2
max( table(d1[ complete.cases(d1$drs_tot) , ]$id) ) # 6


# ----------- session info -----------

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/data_description.txt" )
