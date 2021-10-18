# Program Information  ----------------------------------------------------

# Program:      run_sccs
# Author:       Ema Alsina; Anna Schultze 
# Description:  run sccs on already formatted data, format and export output
# Requirements: 
#               input:  data in csv with 1 row per event [confirm folder]
#                       assumes doses in separate columns 
#               output: table2, table3, log [confirm folder]
              
# Housekeeping  -----------------------------------------------------------

# install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, data.table, here, SCCS, gtsummary, ggplot2, survivalAnalysis)

# save output and error message as a log *sink to be entered 

# Import Data -------------------------------------------------------------
# need to confirm folder structure and input data format 
# input <- fread(here::here("datafolder","inputdata"), data.table = FALSE, na.strings = "")
# in the interim use data from SCCS package and rename variables as applicable 
sample_data <- SCCS::dtpdat

sample_data <- sample_data %>% 
  rename(study_start = sta, 
         censor = end, 
         myocarditis = conv,
         vaccd1 = dtp, 
         vaccd3 = dtpd2) %>% 
  select(-c(dtpd3,sex)) %>% 
  # all times anticipated to be relative to study start 
  mutate(across(c(study_start, censor, myocarditis, vaccd1, vaccd3), ~ .x - study_start)) %>% 
  # negative days assumed to be due to random data generation and filtered out (confirm, seems odd)
  filter(vaccd1 > 0) %>% 
  # make vaccine dates slightly greater to fit pre-exposure window in our analyses 
  mutate(across(c(vaccd1, vaccd3), ~.x + 90)) 

# Data Management ---------------------------------------------------------
# some data management required before fitting model 
# should add controls for variable type here when dummy data received 
# need to consider if risk windows will be split further 
# Levels of interest: 
# 0 - Control Window 
# 1 - Pre-exposure window 
# 2 - Day 0 Dose 1 
# 3 - Dose 1 Risk Window 
# 4 - Time between Dose 1 and 2 
# 5 - Day 0 Dose 2 
# 6 - Dose 2 Risk Window 

# note; assume convention is negative days are exclusive, positive inclusive 
tidy_data <- sample_data %>% 
  # start of control window 
  mutate(start_control = vaccd1 - 91) %>% 
  # start of pre-exposure window 
  mutate(vaccd = vaccd1 - 31) %>% 
  # start of washout
  mutate(vaccd2 = vaccd1 + 28) %>% 
  # constrain end to be end of second risk window 
  mutate(end_risk = case_when(!is.na(vaccd3) ~ vaccd3 + 28, 
                                TRUE ~ vaccd1 + 28)) %>% 
  rowwise %>% 
  mutate(end_study = min(end_risk, censor)) %>% 
  ungroup() %>% 
  # keep those with events 
  filter(myocarditis >= start_control & myocarditis <= end_study)


  # calculate number of 'control' (during control period) and 'events' (during risk period)
  # this is relevant for tables and figures later
tidy_data2 <- tidy_data %>%
  mutate(
    # create one variable for start second risk window (should not be empty), which is either date of 2nd vacc or censor date
    vaccd3_censor = ifelse(is.na(vaccd3), end_study, vaccd3),
    # create logical whether individual has event as control or 'exposed'
    ncontrol = ifelse(between(myocarditis, start_control, vaccd1-1), 1, 0),
    nevent1 = ifelse(between(myocarditis, vaccd1, vaccd2) , 1, 0),
    nevent2 = ifelse(between(myocarditis, vaccd3_censor, end_risk), 1, 0))

# create a vector which is every 30 days since study start until end (last week not complete)
# note 60 day in dummy data to ensure convergence, need to investigate 
min_start <- min(tidy_data$start_control)
max_end <- max(tidy_data$end_study)
calendar_time <- seq(from = min_start+60, to = max_end, by = 60)

# Descriptives ------------------------------------------------------------
# format data to check length of intervals and nr of events per interval 
# need to double check interval specification, have assumed adrug specification is inclusive, aedrug exclusive 
# unclear why final risk windows ends up 29 days; may be clearer to specify all these in variables... 
sccs_data <- formatdata(indiv = case, 
                        astart = start_control, 
                        aend = end_study, 
                        aevent = myocarditis, 
                        adrug = cbind(vaccd, vaccd1-1, vaccd1, vaccd2, vaccd3-1, vaccd3), 
                        aedrug = cbind(vaccd1-1, vaccd1, vaccd1+28, vaccd3-1, vaccd3, vaccd3+28), 
                        dataformat = "multi",
                        sameexpopar = F, 
                        data = tidy_data) 

# check interval distribution and event nr per interval 
interval_distribution <- sccs_data %>%
  select(vaccd, interval, event) %>% 
  tbl_summary(by = vaccd, 
              type = all_continuous() ~ "continuous2", 
              statistic = list(all_continuous() ~ c("{median} ({p25}, {p75})", 
                                                    "{min}, {max}"),
                               all_categorical() ~ "{n}")) %>% 
  modify_header(stat_by = "**{level}**") 

# need to either print to move code to markdown, print to console or export 
interval_distribution

# Unadjusted Analyses -----------------------------------------------------
# precedence given to most recent risk interval - default when data is in wide ("multi") format 

myocard.mod1 <- standardsccs(event ~ vaccd, 
                             indiv = case, 
                             astart = start_control, 
                             aend = end_study, 
                             aevent = myocarditis, 
                             adrug = cbind(vaccd, vaccd1-1, vaccd1, vaccd2, vaccd3-1, vaccd3), 
                             aedrug = cbind(vaccd1-1, vaccd1, vaccd1+28, vaccd3-1, vaccd3, vaccd3+28), 
                             dataformat = "multi",
                             sameexpopar = F, 
                             data = tidy_data) 

myocard.mod1 

# extract relevant results and components for meta-analysis from object 
output.mod1 <- survivalAnalysis::cox_as_data_frame(myocard.mod1)
output.meta.mod1 <- as.data.frame(myocard.mod1$coefficients)

# Adjusted Analyses -------------------------------------------------------

# Note - calendar time adjustment is added through a variable called age in this function, as per Farrington et al 2019. 
myocard.mod2 <- standardsccs(event ~ vaccd + age, 
                             indiv = case, 
                             astart = start_control, 
                             aend = end_study, 
                             aevent = myocarditis, 
                             adrug = cbind(vaccd, vaccd1-1, vaccd1, vaccd2, vaccd3-1, vaccd3), 
                             aedrug = cbind(vaccd1-1, vaccd1, vaccd1+28, vaccd3-1, vaccd3, vaccd3+28), 
                             agegrp = calendar_time, 
                             dataformat = "multi",
                             sameexpopar = F, 
                             data = tidy_data) 
myocard.mod2 

# extract relevant results and components for meta-analysis from object 
output.mod2 <- survivalAnalysis::cox_as_data_frame(myocard.mod2)
output.meta.mod2 <- as.data.frame(myocard.mod2$coefficients) 

# Format and Output Results  ----------------------------------------------

# placeholder - need to format outputs, merge in event counts (don't think this is in model output)
# placeholder - need to apply redaction to small N 

output.meta.mod1 <- output.meta.mod1 %>% 
  rename(yi = `exp(coef)`,
         sei = `se(coef)`) %>%
  select(yi, sei) %>%
  mutate(
    ncase = myocard.mod1$nevent,
    nevent1 = sum(tidy_data2$nevent1),
    nevent2 = sum(tidy_data2$nevent2),
    ncontrol = sum(tidy_data2$ncontrol),
    # I think the ones below won't need to be hardcoded in the final version as this information is
    # included in the analytical datasets as variables we can just copy from an earlier dataframe in this script
    # similar to what I did with the follow-up time and number of events above
    dlab = "CPRD",
    vaccine_dose = factor(x = c(1,2,3,4,5,6), levels = c(1,2,3,4,5,6), 
                          labels = c("pre-exposure", "dose1_day0", "dose1risk", "inbetween", "dose2_day0", "dose2risk")),
    vaccine_brand = "moderna",
    subset ="all"
  )

write.csv(output.meta.mod1, "output.meta.mod1.csv", row.names = FALSE)

# Sensitivity Analyses ----------------------------------------------------
## PLACEHOLDER - WHICH ONES TO PRIORITISE 

# send output back to screen
sink()
sink(type="message")
