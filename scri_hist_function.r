#creates histogram of events over time relative to exposure
#one vaccine dose: tidy_data$vaccd then multiple risk window cut offs (vaaccd1, 2,3)

# mock data
# myocarditis<-sample(1:200, 100)
# vaccd<-sample(1:200,100)
# tidy_data<-as.data.frame(cbind(myocarditis, vaccd))

scri_hist<-function(data=tidy_data, exposure="vaccd", event="myocarditis",cutoffs=c(-31, 1, 0, 30, 60, 100), 
         cutnames=c("control", "buffer","day 0", "risk 1", "risk 2", "risk 3"))
{
  time_from_exp<-data[[event]]-data[[exposure]]
  b <- min(time_from_exp)  # Set the minimum for the breakpoints
  e <- max(time_from_exp) # Set the maximum for the breakpoints
  ax <- pretty(b:e, n = (round(((e-b)/1),digits=0))) # Make a neat vector for the breakpoints
  
  hg_all<-hist(time_from_exp,breaks=ax, plot=F)
  ymax<-(max(hg_all$counts)+3)
  
  hg_unex1 <- hist((time_from_exp[time_from_exp<=cutoffs[1]]), breaks = ax, plot = FALSE) 
  
  hg_burnin <- hist((time_from_exp[between(time_from_exp, (cutoffs[1]+1), cutoffs[2])]), breaks = ax, plot = FALSE)
  
  hg_day0 <- hist((time_from_exp[between(time_from_exp, (cutoffs[2]+1), cutoffs[3])]), breaks = ax, plot = FALSE)
  
  hg_ex1 <- hist((time_from_exp[between(time_from_exp, (cutoffs[3]+1), cutoffs[4])]), breaks = ax, plot = FALSE) 
  
  hg_ex2 <- hist((time_from_exp[between(time_from_exp, (cutoffs[4]+1), cutoffs[5])]), breaks = ax, plot = FALSE)
  
  hg_ex3 <- hist((time_from_exp[between(time_from_exp, (cutoffs[5]+1), cutoffs[6])]), breaks = ax, plot = FALSE)
  
  hg_unex2 <- hist((time_from_exp[time_from_exp>=(cutoffs[6]+1)]), breaks = ax, plot = FALSE) 
  
  # xlim = c(min(time_from_exp), max(time_from_exp)), ylim(c(0, (max(hg_all$counts)+2))),                                   
  
  plot(hg_unex1, col = "grey", main="SCRI event frequency", cex.main=2, ylim=c(0,ymax), xlab="Time from 1st Vaccination", cex.lab=1.5) # Plot 1st histogram using a transparent color
  plot(hg_burnin, col="blue", add=T)
  plot(hg_day0, col=1, add=T)
  plot(hg_ex1, col="red", add=T)
  plot(hg_ex2, col="orange", add=T)
  plot(hg_ex3, col="yellow", add=T)
  plot(hg_unex2, col="grey", add=T)
  abline(v = 0, lwd=3, lty=3)
  legend("topright", pch=16, cex=1.5, pt.cex=2, col=c("grey", "lightblue", "black","red", "orange", "yellow"),
         cutnames, bty="n")
  
  
}

scri_hist(tidy_data)

