# required packages
library(MASS)
library(gamlss.dist)
library(gamlss)
### 
sims<-25000 # number of simulation runs
alpha0<-0.0027 # nominal FAR
# parameters of the re-parameterized Beta distribution
# here (mu,sigma) = (0.15,0.20)
v01<-0.15 # mu
v02<-0.20 # sigma
#######################
m<-30 # size of the Phase I sample
B<-matrix(0,ncol=2,nrow=sims) # empty matrix
# generate 25000 random samples of size m
# estimate mu and sigma
# store the estimates in matrix B
for(j in 1:sims){
  x<-rBE(m,mu=v01,sigma=v02)
  A<-gamlss(x~1,family=BE,trace=F)
  v1est<-1-1/(1+exp(A$mu.coefficients[[1]][1]))
  v2est<-1-1/(1+exp(A$sigma.coefficients[[1]][1]))
  B[j,1]<-v1est # estimate of mu
  B[j,2]<-v2est # estimate of sigma
}
#########################################
listCARL0<-c() # empty vector to store the 25000 CARL0 values
for(j0 in 1:sims){
  # plug-in limits
  UCLest<-qBE(1-alpha0/2,mu=B[j0,1],sigma=B[j0,2])
  LCLest<-qBE(alpha0/2,mu=B[j0,1],sigma=B[j0,2])
  # conditional false alarm rate
  alpha0est<-pBE(LCLest,mu=v01,sigma=v02)+1-pBE(UCLest,mu=v01,sigma=v02)
  # IC CARL
  CARL0<-1/alpha0est
  # append the CARL0 in the vector
  listCARL0<-c(listCARL0,CARL0)
}
#############################
AARL<-mean(listCARL0) # estimate the E(CARL0) as the sample mean of the CARL0 values
SDARL<-sd(listCARL0) # estimate the sqrt(V(CARL0)) as the standard deviation of the CARL0 values
perc<-sum(listCARL0<1/alpha0)/sims # relative frequency of the cases CARL0<1/alpha0 = 370.4
#############################
# print the results
paste('The mu is',v01,' and the sigma is',v02)
cat("m:",m,"FAR0:",alpha0,'AARL0:',AARL,'SDARL0:',SDARL,'perc:',perc*100,'%','\n')
print('percentiles of the distribution of CARL0')
print(quantile(listCARL0,c(0.05,0.25,0.50,0.75,0.95)))


