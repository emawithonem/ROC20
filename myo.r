#Author: Ema Alsina, MSc.
#email: e.m.alsina-2@umcutrecht.nl
#Organisation: UMC Utrecht, Utrecht, The Netherlands
#Date: 12/10/2021

#demo myocarditis SCRI model
#SCRI is because the pre-exposure (first vaccine) is set in the data

myo.mod1 <- standardsccs(indiv=case, astart=start_date, aend=end_date,
                          aevent=myocarditis, adrug=cbind(dose1,dose2),
                          aedrug=cbind(rv+21,rvd2+21), data=myodat)