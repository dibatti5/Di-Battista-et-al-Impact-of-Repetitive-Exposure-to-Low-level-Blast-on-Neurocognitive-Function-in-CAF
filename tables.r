###Code to create tables used in manuscript

#libraries
library(gtsummary)
library(gt)

#functions
HDILow<- function(x, HDI=0.90) {
  sortedPts = sort( x)
  ciIdxInc = ceiling( HDI * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc 
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ] }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ] 
  return( HDImin)
}
HDIHigh<- function(x, HDI=0.90) {
  sortedPts = sort( x)
  ciIdxInc = ceiling( HDI * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc 
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ] }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ] 
  return( HDImax)
}

#load data (Not available in public repository version for privacy reasons)
#df <- read.csv('data_df.csv', stringsAsFactors = FALSE)

### Table Creation ###
#create a basic demos table, and then a table of the variables of interest

#Demos
colnames(df)
demo_df <- df[c(46, 1:8)]

#relable demo_df according to the variable list sent with the original excel file.
colnames(demo_df)
demo_df$Gender <- ifelse(demo_df$Gender==1,"Male",'Female')
demo_df$Status <- ifelse(demo_df$Status==1,'Regular Force','Primary Reserve')

# Renaming rank for table
demo_df$Rank[demo_df$Rank==1] <- 'Junior NCM'
demo_df$Rank[demo_df$Rank=='2']  <- 'Senior'
demo_df$Rank[demo_df$Rank=='3'] <- 'Subordinate Officer'
demo_df$Rank[demo_df$Rank=='4'] <- 'Junior Officer'
demo_df$Rank[demo_df$Rank=='5'] <- 'Senior Officer'
demo_df$Rank[demo_df$Rank=='6'] <- 'General Officer'
colnames(demo_df)[7] <- 'Years of Service'
demo_df <- demo_df[c(1, 3, 4, 7, 6, 5)]

# Reorder ranks for table
demo_df$Rank <- factor(demo_df$Rank, levels=c('Junior NCM','Senior',
                                              'Subordinate Officer','Junior Officer',
                                              'Senior Officer','General Officer'))

# Ordering group for table
demo_df$group <- factor(demo_df$group, levels = c("Breacher", "Sniper",
                                                  "Military Control"))


#create table with gt summary
t1 <- tbl_summary(demo_df, by = group,missing = 'no',
                  digits = all_continuous() ~ 1) %>%
  bold_labels() %>%
  modify_caption("**Table 1. Participant Demographics**") %>%
  as_gt() %>%
  gt::gtsave(filename = "t1.html") 

#relevant variables table
#removing some subscale items for brevity
colnames(df)
t2_df <- df[c(46, 23, 24, 29, 25, 26, 27, 7, 8, 12, 13, 17:22)]

#renaming binary variables
colnames(t2_df)
t2_df[c(2:7)][t2_df[c(2:7)]==1] <- 'No'
t2_df[c(2:7)][t2_df[c(2:7)]=='2'] <- 'Yes'

## If putting groups together
t2_df$group <- demo_df$group

#renaming variables
colnames(t2_df) <- c("Group", "War Zone Deployment", "Concussion", "Blast", 
                    "Impact", "MVA", "Fell", "RPQ 3", "RPQ 13", "PCL-5", "BSI", 
                    "4-choice RT task (mean RT correct)", 
                    "Delayed matching-to-sample task (% accurate)", "1-back",
                    "2-back", "3-back", "Stroop (RT difference)")


t2 <- tbl_summary(t2_df,by = Group, missing = 'no',
                  digits = all_continuous() ~ 2) %>%
  modify_caption("**Table 2. Psychological, Neurological and Brain Injury Measures**") %>%
  as_gt() %>%
  tab_row_group(
    label = "Neurocognitive",
    rows = 11:16) %>%
  tab_row_group(
    label = "Mental Health",
    rows = 9:10) %>%
  tab_row_group(
    label = "Concussion Symptoms",
    rows = 7:8) %>%
  tab_row_group(
    label = "History of Head Injury",
    rows = 1:6) %>%
  tab_options(row_group.font.weight = 'bold')%>%
  tab_source_note(
    source_note = md('MVA, motor vehicle accident; RPQ, rivermead post concussion
                      symptoms questionnaire; PCL, Posttraumatic stress disorder 
                      checklist; BSI, brief symptom inventory; RT, reaction time'))%>%
  gt::gtsave(filename = "t2.html") 
