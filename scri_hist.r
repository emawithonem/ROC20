#creates histogram of events over time relative to exposure
#one vaccine dose: tidy_data$vaccd then multiple risk window cut offs (vaaccd1, 2,3)

#1st: center on vaccination exposure T0


#center time on 1st vaccine exposure
time_from_exp<- tidy_data$myocarditis-tidy_data$vaccd

#set cut offs for time categories
control1=-31
burnin= -1
day0<-0
exposure1=30
exposure2=60
exposure3=100
control2= max(time_from_exp)

# set plot dimensions and bin widths

b <- min(time_from_exp)  # Set the minimum for the breakpoints
e <- max(time_from_exp) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = (round(((e-b)/1),digits=0))) # Make a neat vector for the breakpoints

hg_all<-hist(time_from_exp,breaks=ax, plot=F)

ymax<-(max(hg_all$counts)+3)

hg_unex1 <- hist((time_from_exp[time_from_exp<=control1]), breaks = ax, plot = FALSE) 
                  
hg_burnin <- hist((time_from_exp[between(time_from_exp, (control1+1), burnin)]), breaks = ax, plot = FALSE)

hg_day0 <- hist((time_from_exp[(time_from_exp==day0)]), breaks = ax, plot = FALSE)

hg_ex1 <- hist((time_from_exp[between(time_from_exp, (burnin+1), exposure1)]), breaks = ax, plot = FALSE) 

hg_ex2 <- hist((time_from_exp[between(time_from_exp, (exposure1+1), exposure2)]), breaks = ax, plot = FALSE)

hg_ex3 <- hist((time_from_exp[between(time_from_exp, (exposure2+1), exposure3)]), breaks = ax, plot = FALSE)

hg_unex2 <- hist((time_from_exp[time_from_exp>=(exposure3+1)]), breaks = ax, plot = FALSE) 
 
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
       c("control", "buffer","day 0", "risk 1", "risk 2", "risk 3"), bty="n")
