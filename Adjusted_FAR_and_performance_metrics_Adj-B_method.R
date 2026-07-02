# required packages
library(MASS)
library(gamlss.dist)
library(gamlss)
### 
sims<-25000 # number of simulation runs
alpha0<-0.0027 # nominal FAR
ARL0<-1/alpha0 # nominal ARL0 value
eps<-0.2 # the epsilon value
pstar<-0.1 # threshold p', here p'=0.05
# parameters of the re-parameterized Beta distribution
# here (mu,sigma) = (0.75,0.40)
v01<-0.75 # mu
v02<-0.40 # sigma
#######################
m<-100 # size of the Phase I sample
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
# perform a numerical search for the alpha'' value (adjusted FAR)
# we try several values (denoted as a1)
# and the one that satisfies the requirement RF = #{CARL0<((1+eps)^(-1))*370.4}/sims < p' is the alpha''
for(a1 in seq(0.003,0.00002,by=-0.00001)){
  listCARL0<-c() # empty vector to store the 25000 CARL0 values for the given
  for(j0 in 1:sims){
    # plug-in limits for the given a1 value
    UCLest<-qBE(1-a1/2,mu=B[j0,1],sigma=B[j0,2])
    LCLest<-qBE(a1/2,mu=B[j0,1],sigma=B[j0,2])
    # conditional false alarm rate
    alpha0est<-pBE(LCLest,mu=v01,sigma=v02)+1-pBE(UCLest,mu=v01,sigma=v02)
    # IC CARL
    CARL0<-1/alpha0est
    # append the CARL0 in the vector
    listCARL0<-c(listCARL0,CARL0)
  }
  AARL<-mean(listCARL0) # estimate the E(CARL0) as the sample mean of the CARL0 values
  RF<-sum(listCARL0<(370.4/(1+eps)))/sims # calculate the Relative Frequency
  if(RF<pstar){break} # check if the requirement is satisfied
}
# since the requirement is satisfied
# calculate the performance metrics
# for this a1 value
SDARL<-sd(listCARL0) # estimate the sqrt(V(CARL0)) as the standard deviation of the CARL0 values
perc<-sum(listCARL0<1/alpha0)/sims # relative frequency of the cases CARL0<1/alpha0 = 370.4
#############################
# print the results
paste('The mu is',v01,' and the sigma is',v02)
cat("m:",m,"adjusted FAR:",a1,'AARL:',AARL,'SDARL:',SDARL,'perc:',perc*100,'%',"RF:",RF,'\n')
print('percentiles of the distribution of CARL0')
print(quantile(listCARL0,c(0.05,0.25,0.50,0.75,0.95)))
##########################################################