#Author: Ema Alsina, MSc.
#email: e.m.alsina-2@umcutrecht.nl
#Organisation: UMC Utrecht, Utrecht, The Netherlands
#Date: 12/10/2021

#demo myocarditis SCRI model
#SCRI is because the pre-exposure (first vaccine) is set in the data

myo.mod1 <- standardsccs(indiv=case, astart=start_date, aend=end_date,
                          aevent=myocarditis, adrug=cbind(dose1,dose2),
                          aedrug=cbind(rv+21,rvd2+21), data=myodat)

#are any of the SCCS assumptions violated?
  #Event must not influence likelihood of exposure. Example: event (stroke) is a contraindication for the exposure (medication)
    #probably- is a myocarditis/pericarditis diagnosis a contraindication for vaccination?
    #use eventdepenexp model ?
  #Event must not influence observation time. Example: death
    #OK
  #Exposure must not influence observation of event. Example: misinformation making hypervigilant parents look for autism symptoms following vaccination
    #OK
  # Recurrent events must be independent. Example: first stroke increases likelihood of having another
    #OK