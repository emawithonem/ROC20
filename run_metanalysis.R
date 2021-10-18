# Program Information  ----------------------------------------------------

# Program:      run_sccs
# Author:       Ema Alsina; Anna Schultze 
# Description:  run sccs on already formatted data, format and export output
# Requirements: 
#               input:  data in csv with 1 row per event [confirm folder]
#                       assumes doses in separate columns 
#               output: table2, table3, log [confirm folder]

# Housekeeping  -----------------------------------------------------------

## install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, data.table, here, SCCS, gtsummary, ggplot2)

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

tidy_data <- sample_data %>% 
  # start of control window 
  mutate(start_control = vaccd1 - 89) %>% 
  # start of pre-exposure window NOTE: think this might be simplest way of incorporating? treated as level of the exposure var
  mutate(vaccd = vaccd1 - 29) %>% 
  # start of washout
  mutate(vaccd2 = vaccd1 + 29) %>% 
  # constrain end to be end of second risk window 
  mutate(end_risk = case_when(!is.na(vaccd3) ~ vaccd3 + 28, 
                              TRUE ~ vaccd1 + 28)) %>% 
  rowwise %>% 
  mutate(end_study = min(end_risk, censor)) %>% 
  ungroup() %>% 
  # keep those with events 
  filter(myocarditis >= start_control & myocarditis <= end_study)

# Descriptives ------------------------------------------------------------
# format data to check length of intervals and nr of events per interval 
sccs_data <- formatdata(indiv = case, 
                        astart = start_control, 
                        aend = end_study, 
                        aevent = myocarditis, 
                        adrug = cbind(vaccd, vaccd1+1, vaccd2, vaccd3+1), 
                        aedrug = cbind(vaccd1, vaccd1 + 28, vaccd3, vaccd3 + 28), 
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

interval_distribution
# need to either print to move code to markdown, print to console or export 

# Unadjusted Analyses -----------------------------------------------------
# precedence given to most recent risk interval - default when data is in wide ("multi") format 

myocard.mod1 <- standardsccs(event ~ vaccd, 
                             indiv = case, 
                             astart = start_control, 
                             aend = end_study, 
                             aevent = myocarditis, 
                             adrug = cbind(vaccd, vaccd1+1, vaccd2, vaccd3+1), 
                             aedrug = cbind(vaccd1, vaccd1 + 28, vaccd3, vaccd3 + 28), 
                             dataformat = "multi",
                             sameexpopar = F, 
                             data = tidy_data) 
myocard.mod1 

# Adjusted Analyses -------------------------------------------------------
# create a vector which is every 2 weeks since study start until end (last week not complete)
# note, error message if I did not change these to not match astart and aend exactly 
# requirement for "unique" cut-points, may run into issues here - AS to confirm internally at LSHTM
min_start <- min(tidy_data$start_control)
max_end <- max(tidy_data$end_study)
calendar_time <- seq(from = min_start+1, to = max_end-1, by = 7)

myocard.mod2 <- standardsccs(event ~ vaccd, 
                             indiv = case, 
                             astart = start_control, 
                             aend = end_study, 
                             aevent = myocarditis, 
                             adrug = cbind(vaccd, vaccd1+1, vaccd2, vaccd3+1), 
                             aedrug = cbind(vaccd1, vaccd1 + 28, vaccd3, vaccd3 + 28), 
                             agegrp = calendar_time, 
                             dataformat = "multi",
                             sameexpopar = F, 
                             data = tidy_data) 
myocard.mod2 
# Note - need to figure out how to add adjustment into model formula (modeling guide recommends as age, data has no age atm)

# Output Results  ---------------------------------------------------------

# creating dataset
# for now I have just randomly assumed some stuff about vaccine brand and subsets
# I also realise my table in the SAP as it currently is might not be detailed enough.
# And this might be a very laborious and roundabout way to make the dataset, because I hardcoded stuff
# and we'd need to include a variety of results from different vaccines and subgroups
results <- as.data.frame(myocard.mod1$coefficients) %>%
  rename(yi = `exp(coef)`,
         sei = `se(coef)`) %>%
  select(yi, sei) %>%
  mutate(
    ncase = nrow(tidy_data),
    #nevent = sum(ifelse(between(tidy_data$myocarditis, tidy_data$))),
    dlab = c("BIFAP", "PHARMO", "ARS","CPRD"),
    vaccine_dose = c(1,2,3,4),
    vaccine_brand = "moderna",
    subset ="all"
    )

# running the meta-analysis
# for now meta-analysing on vaccine dose but should of course be done across datasets per vaccine and dose
#install.packages("metafor")
library(metafor)

# turn data into right format for meta-analysis function (rma.uni)
meta_data <- escalc(yi = yi, sei = sei, slab = dlab, data = results,  measure = "IRR")

# meta-analyse and create forest plot
# I don't know why the 95%CIs in the forest plot are different from the ones from the sccs function. Have to look into that
meta_analysis <- rma.uni(yi = yi, sei = sei, measure = "IRR", data = meta_data,
                         slab = dlab)

forest(meta_analysis, header = c("Dataset", "IRR (95% CI)"), xlab = "Incidence Rate Ratio",
       mlab = "", ilab = meta_data$ncase,
       xlim = c(-10, 10), ilab.xpos = -5)

op <- par(cex=0.75, font=2)
text(-5, "Cases")

# Sensitivity Analyses ----------------------------------------------------
## PLACEHOLDER - WHICH ONES TO PRIORITISE 

