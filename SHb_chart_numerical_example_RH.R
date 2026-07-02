#library(rattle)
# data1<-as.data.frame(weatherAUS)
#weatherAUS
#data1[data1$Location=="Sydney",]
########################################
# required packages
library(dplyr)
library(readxl)
library(gamlss.dist)
library(gamlss)
####################
RHsydney <- read_excel("RHsydney.xls") # the file is in the default directory
# View(RHsydney) # uncomment to view the data in R
RHsydney1<-as.data.frame(RHsydney) # transform it into a data frame
attach(RHsydney1)
min_humidity <- RHsydney1 %>%
  group_by(Year, Month) %>%
  summarise(
    min_relative_humidity = min(RH, na.rm = TRUE),
    .groups = "drop"
  )
# Print the result
print(min_humidity,n=61)
x<-min_humidity$min_relative_humidity # X is the monthly relative humidity
###############
# create Figure 1
plot(x[1:61],type='l',axes=F,ylim=c(0,60),ylab='Relative Humidity (%)',xlab='time')
axis(1,c(1,10,20,30,40,50,60),
     labels=c("10-2010",'08-2011','06-2012','06-2013','04-2014','02-2015','12-2015'))
axis(2)
abline(v=50,lty=2,col='grey',lwd=2)
points(x[1:61],pch=16)
###############
summary(x) # basic descriptive measures
###############
# use the 50 first values as the Phase I dataset
x0<-x[1:50]/100 # Phase I sample
# estimate mu and sigma
A<-gamlss(x0~1,family=BE,trace=F)
v1est<-1-1/(1+exp(A$mu.coefficients[[1]][1])) # estimation of mu
v2est<-1-1/(1+exp(A$sigma.coefficients[[1]][1])) # estimation of sigma
mu.hat<-round(v1est,digits=5)
sigma.hat<-round(v2est,digits=5)
paste('the estimated mu is:',mu.hat,' and the estimated sigma is',sigma.hat)
###############################
# plug-in limits
alpha<-0.0027 # nominal FAR
UCLest<-qBE(1-alpha/2,mu=mu.hat,sigma=sigma.hat)
LCLest<-qBE(alpha/2,mu=mu.hat,sigma=sigma.hat)
# Plot the Phase I data against the plug in limits
plot(x[1:61],type='l',axes=F,ylim=c(0,90),
     ylab='Relative Humidity (%)',xlab='time',
     main= expression(SH[B] * "-chart for Sydney RH data "))
axis(1,c(1,10,20,30,40,50,60),
     labels=c("01-2010",'08-2011','06-2012','06-2013','04-2014','02-2015','12-2015'))
axis(2)
abline(v=50,lty=1,lwd=1)
points(x[1:61],pch=16,cex=1.5)
abline(h=(100*UCLest),lty=1,lwd=1)
abline(h=(100*LCLest),lty=1,lwd=1)
# Adjustment method-A, find limits & plot the data
aa1<-0.00340 # use the a' from the respective table in the paper
UCL1<-qBE(1-aa1/2,mu=mu.hat,sigma=sigma.hat)
LCL1<-qBE(aa1/2,mu=mu.hat,sigma=sigma.hat)
abline(h=(100*UCL1),lty=2,lwd=1)
abline(h=(100*LCL1),lty=2,lwd=1)
# Adjustment method-B, epsilon=0, find limits & plot the data
aa2<-0.00022 # use the a'' from the respective table in the paper
UCL2<-qBE(1-aa2/2,mu=mu.hat,sigma=sigma.hat)
LCL2<-qBE(aa2/2,mu=mu.hat,sigma=sigma.hat)
abline(h=(100*UCL2),lty=3,lwd=1)
abline(h=(100*LCL2),lty=3,lwd=1)
# Adjustment method-B, epsilon=0.2, find limits & plot the data
aa3<-0.00052 # use the a'' from the respective table in the paper
UCL3<-qBE(1-aa3/2,mu=mu.hat,sigma=sigma.hat)
LCL3<-qBE(aa3/2,mu=mu.hat,sigma=sigma.hat)
abline(h=(100*UCL3),lty=4,lwd=1)
abline(h=(100*LCL3),lty=4,lwd=1)
# Bootstrap limits
# below the code to obtain them
simsB<-25000 # B=25000 bootstrap samples from the Phase I dataset (resampling)
B<-matrix(0,ncol=4,nrow=simsB) # empty matrix
# we simulate 25000 bootstrap samples
# and from each sample we estimate mu and sigma
# for each pair of estimate, a pair of control limits is calculated
# in columns 3 and 4, LCL and UCL, respectively, of matrix B
for(j in 1:simsB){
  xB<-sample(x0,size=length(x0),replace=T)
  AB<-gamlss(xB~1,family=BE,trace=F)
  v1estB<-1-1/(1+exp(AB$mu.coefficients[[1]][1]))
  v2estB<-1-1/(1+exp(AB$sigma.coefficients[[1]][1]))
  B[j,1]<-v1estB
  B[j,2]<-v2estB
  UCLestB<-qBE(1-alpha/2,mu=B[j,1],sigma=B[j,2])
  LCLestB<-qBE(alpha/2,mu=B[j,1],sigma=B[j,2])
  B[j,3]<-LCLestB
  B[j,4]<-UCLestB
}
# use the 0.025-percentile of the 25000 LCL values as the Bootstrap-based LCL
LCLstarB<-quantile(B[,3],c(0.025))
# use the 0.975-percentile of the 25000 UCL values as the Bootstrap-based UCL
UCLstarB<-quantile(B[,4],c(0.975))
# add lines on the chart
abline(h=(100*UCLstarB),lty=5,lwd=1)
abline(h=(100*LCLstarB),lty=5,lwd=1)
##########################
# add a legend
legend('topright',c('Plug-in','AdjA','AdjB,eps=0','AdjB,eps=0.2','Bootstrap'),
lty=1:5,lwd=c(1,1,1,1,1),cex=0.75)
text(42,70,"Phase I")
text(53,70,"Phase II")
#######################
# below we provide the OOC scenario when there is a shift only in mu
# OOC data
y<-c(0.30609619, 0.25720683, 0.29109522,
     0.37668141, 0.50084142, 0.35030771,
     0.38610248, 0.31420693, 0.52388937,
     0.21405942, 0.17209139, 0.14343708,
     0.31683782, 0.11597719, 0.15894399,
     0.31830692, 0.33706144, 0.16367507,
     0.27481577, 0.15476869, 0.53893121,
     0.31816739, 0.18837781, 0.40100539,
     0.38775144, 0.36651343, 0.19352635,
     0.20936857, 0.33370933, 0.22000786,
     0.18772538, 0.56640794, 0.16146448,
     0.33519866, 0.44298871, 0.12165616,
     0.47767826, 0.35177966, 0.22916423,
     0.06598294, 0.10970764, 0.42896833,
     0.28774871, 0.34526811, 0.09874269,
     0.17692411, 0.35719114, 0.30169661,
     0.24562682, 0.57983994, 0.32418274,
     0.11011865, 0.37507631, 0.44071340,
     0.16974362, 0.21863353, 0.15890170,
     0.16431547, 0.32858313, 0.30379717,
     0.22211175, 0.59030055, 0.24904023,
     0.35161064, 0.34471037, 0.17930245,
     0.39895906, 0.17159297, 0.18355666,
     0.21480203)
###########################################
# create the chart
plot(y,type='l',axes=F,ylim=c(0,0.9),
     xlim=c(0,70),ylab='',xlab='time',main= expression(SH[B] * "-chart for OOC data, shift in " * mu))
axis(1,c(0,10,20,30,40,50,60,70))
axis(2,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
points(y,pch=16,cex=1.5)
# plug-in limits
abline(h=(UCLest),lty=1,lwd=1)
abline(h=(LCLest),lty=1,lwd=1)
# Adjustement-A method limits
abline(h=(UCL1),lty=2,lwd=1)
abline(h=(LCL1),lty=2,lwd=1)
# Adjustment-B method limits, eps=0
abline(h=(UCL2),lty=3,lwd=1)
abline(h=(LCL2),lty=3,lwd=1)
# Adjustment-B method limits, eps=0.2
abline(h=(UCL3),lty=4,lwd=1)
abline(h=(LCL3),lty=4,lwd=1)
# Bootstrap-based limits
abline(h=(UCLstarB),lty=5,lwd=1)
abline(h=(LCLstarB),lty=5,lwd=1)
# Add a Legend
legend('topright',c('Plug-in','AdjA','AdjB,eps=0','AdjB,eps=0.2','Bootstrap'),
       lty=1:5,lwd=c(1,1,1,1,1),cex=0.75)
abline(v=20,lty=2,col='darkgrey')
text(12,0.75,"IC period")
text(28,0.75,"OOC period")
####################################
# below we provide the OOC scenario when there is a shift only in sigma
# OOC data
y<-c(0.40890923, 0.20711377, 0.48565290,
     0.41342518, 0.14514276, 0.28553799,
     0.36154274, 0.17148259, 0.25833854,
     0.44677121, 0.20605114, 0.20945256,
     0.33630733, 0.33168441, 0.35256853,
     0.47304573, 0.17692752, 0.19831857,
     0.23240015, 0.23455703, 0.13494582,
     0.21935793, 0.25749392, 0.26100636,
     0.09398065, 0.22532295, 0.35740527,
     0.37906132, 0.17381027, 0.10006924,
     0.18711972, 0.24825863, 0.31336222,
     0.16107138, 0.24012080, 0.40356519,
     0.16573125, 0.48085574, 0.13661168,
     0.17103081, 0.22010640, 0.50253456,
     0.16494011, 0.60033655, 0.30104390,
     0.22931254, 0.59511069, 0.18845682,
     0.34120781, 0.56904273, 0.20563767,
     0.64621361, 0.24440051, 0.20950493,
     0.12702669, 0.18782709, 0.40409009,
     0.15381170, 0.33829964, 0.26553242,
     0.11436794, 0.32380358, 0.23113463,
     0.20722643, 0.47252334, 0.26407311,
     0.19531615, 0.43436949, 0.04664434,
     0.16630466)
###########################################
# create the chart
plot(y,type='l',axes=F,ylim=c(0,0.9),
     xlim=c(0,70),ylab='',xlab='time',main= expression(SH[B] * "-chart for OOC data, shift in " * sigma))
axis(1,c(0,10,20,30,40,50,60,70))
axis(2,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
points(y,pch=16,cex=1.5)
# plug-in limits
abline(h=(UCLest),lty=1,lwd=1)
abline(h=(LCLest),lty=1,lwd=1)
# Adjustement-A method limits
abline(h=(UCL1),lty=2,lwd=1)
abline(h=(LCL1),lty=2,lwd=1)
# Adjustment-B method limits, eps=0
abline(h=(UCL2),lty=3,lwd=1)
abline(h=(LCL2),lty=3,lwd=1)
# Adjustment-B method limits, eps=0.2
abline(h=(UCL3),lty=4,lwd=1)
abline(h=(LCL3),lty=4,lwd=1)
# Bootstrap-based limits
abline(h=(UCLstarB),lty=5,lwd=1)
abline(h=(LCLstarB),lty=5,lwd=1)
# Add a Legend
legend('topright',c('Plug-in','AdjA','AdjB,eps=0','AdjB,eps=0.2','Bootstrap'),
       lty=1:5,lwd=c(1,1,1,1,1),cex=0.75)
abline(v=20,lty=2,col='darkgrey')
text(12,0.75,"IC period")
text(28,0.75,"OOC period")