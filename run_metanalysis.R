# Program Information  ----------------------------------------------------

# Program:      run_metaanalysis
# Author:       Sophie Bots
# Description:  meta-analyse aggregated scri output data
# Requirements: 
#               input:  data in csv with 1 row per analysis stratum
#               output: table 1, table 2, figures

# Housekeeping  -----------------------------------------------------------

## install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, meta)

# Import Data -------------------------------------------------------------
# need to confirm folder structure on the anDREa platform
# I also imagine each DAP uploads their CSV separately, so I will need to mix and match a bit
# I can write a suggestion, but I have no clue if this is how it'll work so probably
# will need to adapt when the actual data are there

## NEEDS FIXING ##
# because CSV output acts weirdly on my local machine, just copied the all_results that was in my R environment
# but this really needs to be fixed after we decided on updating output code
# pretending it's from different locations so I can write template code for that
bifap_data <- all_results
pharmo_data <- all_results
ars_data <- all_results
cprd_data <- all_results

# Format Data -------------------------------------------------------------

# adding dlab variable
bifap_data <- bifap_data %>%
  mutate(dlab = "BIFAP")

pharmo_data <- pharmo_data %>%
  mutate(dlab = "PHARMO")

ars_data <- ars_data %>%
  mutate(dlab = "ARS")

cprd_data <- cprd_data %>%
  mutate(dlab = "CPRD")

# add all together
scri_data <- rbind(bifap_data, pharmo_data, ars_data, cprd_data)

# Running the meta-analysis -----------------------------------------------
# meta-analysis should be done:
# (a) for each main analysis subtype
# (b) for each sensitivity analysis subtype
# for both: meta-analyse across DAPs separately for each
# - vaccine type
# - risk window
# - adjustment type (unadjusted / adjusted for seasonality)

## TABLE 1 ------------------------------------------------

# write function so it's easy to change settings without having to copy-paste code
create_tab1 <- function(adjustment, riskwindow, outcome) {
  
  meta_analysis <- metagen(data = scri_data %>% filter(analysis == adjustment & label == riskwindow &
                                                    eventtype == outcome),
                           TE = yi, seTE = sei, studlab = dlab,
                           sm = "IRR", lower = lci, upper = uci,
                           random = T, fixed = F,
                           subgroup = vacctype)
  return(meta_analysis)
}

# loop over the different outcomes
# this can also easily be adjusted to accommodate more outcomes
outcomes <- c("Myocarditis", "Pericarditis")
tab1 <- vector(mode = "list", length = length(outcomes))
names(tab1) <- outcomes


for (i in 1:length(outcomes)) {
  
  # Run analysis
  u_dose1 <- create_tab1(adjustment = "unadjusted", riskwindow = "dose 1 risk window", outcome = outcomes[i])
  u_dose2 <- create_tab1(adjustment = "unadjusted", riskwindow = "dose 2 risk window", outcome = outcomes[i])
  a_dose1 <- create_tab1(adjustment = "adjusted", riskwindow = "dose 1 risk window", outcome = outcomes[i])
  a_dose2 <- create_tab1(adjustment = "adjusted", riskwindow = "dose 2 risk window", outcome = outcomes[i])
  
  # make the table
  tab <- data.frame(
    vacc = u_dose1$bylevs,
    irr_1 = exp(u_dose1$TE.random.w),
    lci_1 = exp(u_dose1$lower.random.w),
    uci_1 = exp(u_dose1$upper.random.w),
    irr_2 = exp(u_dose2$TE.random.w),
    lci_2 = exp(u_dose2$lower.random.w),
    uci_2 = exp(u_dose2$upper.random.w))
    
  # save table
  tab1[[i]] <- tab
  
}

# output tables as csv or something else, depends on how anDREa platform works

## FOREST PLOTS --------------------------------------------------------------------
# placeholder now -- needs to be updated to actual data and functions


# first turn data into required format for plotting function
meta_analysis <- metagen(data = results, TE = yi, seTE = sei, studlab = dlab, sm = "IRR",
                          lower = lci, upper = uci,
                          random = T, fixed = F, n.c = ncase, n.e = nevent,
                          label.e = "", label.c = "",
                          label.left = "lower risk", label.right = "higher risk")

# then plot
forest.meta(meta_analysis2,
            sortvar = TE,
            smlab = "Incidence Rate Ratio",
            print.tau2 = F,
            leftcols = c("dlab", "ncase", "nevent"),
            leftlabs = c("Data source", "Cases", "Events"),
            just.addcols = "left")

# Sensitivity Analyses ----------------------------------------------------


