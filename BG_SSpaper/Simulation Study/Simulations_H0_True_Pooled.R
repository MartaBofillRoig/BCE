
################################################################
# SIMULATIONS SAMPLE SIZE FOR COMPOSITE BINARY ENDPOINTS
# Simulating H0 true (pooled variance)
# oct-2018  
# Marta Bofill and Guadalupe Gómez
################################################################

# Reading the database
######################
setwd("C:/Users/marta.bofill/Simulations")
data<-read.table(file="./DATABASE_Scenarios_H0_True.csv",sep = ";", header = T, dec=",")


# SIMULATION FUNCTION
# without sampling
#######################

# PREAMBLE
###########
alpha=0.025; beta=0.2
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1)


library(snowfall)
library(snow)

# Init Snowfall with explicit settings: Initialisation of cluster usage
sfInit( parallel=TRUE, cpus=5 )

# FUNCTION
###########
# Simulating and calculating the statistic test 
set.seed(321)
# nsim: number of simulations
nsim=100000

################################################################
# RISK DIFFERENCE
f <- function(samplesize,s1_group0,s2_group0,s3_group0,s1_group1,s2_group1,s3_group1){
  
  # estimation of the probability of observing the composite event in control group
  R = runif(samplesize)
  phat_group0 = 1 - sum(R>=(s1_group0+s2_group0+s3_group0) & R<1)/samplesize
  
  # estimation of the probability of observing the composite event in test group
  R = runif(samplesize)
  phat_group1 = 1 - sum(R>=(s1_group1+s2_group1+s3_group1) & R<1)/samplesize
  
  # estimation of the probability under the null hypotheses
  phat_pooled = (phat_group0 + phat_group1)/2
  
  # test risk difference with pooled variance
  TestD_pooled = (phat_group1-phat_group0)/sqrt(( 2*phat_pooled*(1-phat_pooled))/samplesize) 
  
  return(TestD_pooled)
}

################################################################
# RISK RATIO
fRR <- function(samplesize,s1_group0,s2_group0,s3_group0,s1_group1,s2_group1,s3_group1){
  
  # estimation of the probability of observing the composite event in control group
  R = runif(samplesize)
  phat_group0 = 1 - sum(R>=(s1_group0+s2_group0+s3_group0) & R<1)/samplesize
  
  # estimation of the probability of observing the composite event in test group
  R = runif(samplesize)
  phat_group1 = 1 - sum(R>=(s1_group1+s2_group1+s3_group1) & R<1)/samplesize
  
  # estimation of the probability under the null hypotheses
  phat_pooled = (phat_group0 + phat_group1)/2
  
  # test risk difference with pooled variance
  TestR_pooled = log(phat_group1/phat_group0)*( 2*(1-phat_pooled)/(samplesize*phat_pooled) )^(-1/2) 
  
  return(TestR_pooled)
}

################################################################
# ODDS RATIO
fOR <- function(samplesize,s1_group0,s2_group0,s3_group0,s1_group1,s2_group1,s3_group1){
  
  # estimation of the probability of observing the composite event in control group
  R = runif(samplesize)
  phat_group0 = 1 - sum(R>=(s1_group0+s2_group0+s3_group0) & R<1)/samplesize
  
  # estimation of the probability of observing the composite event in test group
  R = runif(samplesize)
  phat_group1 = 1 - sum(R>=(s1_group1+s2_group1+s3_group1) & R<1)/samplesize
  
  # estimation of the probability under the null hypotheses
  phat_pooled = (phat_group0 + phat_group1)/2
  
  # test risk difference with pooled variance
  TestOR_pooled = log((phat_group1/(1-phat_group1))/(phat_group0/(1-phat_group0)))*(2/(samplesize*phat_pooled*(1-phat_pooled)))^(-1/2)
  
  return(TestOR_pooled)
}


################################################################

data$Test_Reject=0 
data$Test_Reject_RR=0 
data$Test_Reject_OR=0 

t0=Sys.time()
for(i in 1:dim(data)[1]){
  # corr: correlation 
  corr= data$rho_true[i]
  
  # Group 0 (Control Group)
  ############################ 
  
  # p1_group0, p2_group0: marginal probabilities for the single endpoints in the control group
  p1_group0 = data$p1.0[i]
  p2_group0 = data$p2.0[i]
  
  # phi_group0: Joint problability two correlated binary endpoints P(X1 =1, X2 =1)
  phi_group0 = corr*sqrt(p1_group0*(1-p1_group0)*p2_group0*(1-p2_group0)) + p1_group0*p2_group0
  
  # 2x2 Table of probabilities
  s1_group0 = phi_group0
  s2_group0 = p1_group0 - phi_group0
  s3_group0 = p2_group0 - phi_group0
  s4_group0 = 1 - s1_group0 - s2_group0 - s3_group0 
  
  # Group 1 (Treatment group)
  ############################
  
  # p1_group1, p2_group1: marginal probabilities for the single endpoints in the test group
  p1_group1 = p1_group0
  p2_group1 = p2_group0
  
  # phi_group1: Joint problability two correlated binary endpoints P(X1 =1, X2 =1)
  phi_group1 = phi_group0
  
  # 2x2 Table of probabilities
  s1_group1 = phi_group1
  s2_group1 = p1_group1 - phi_group1
  s3_group1 = p2_group1 - phi_group1
  s4_group1 = 1 - s1_group1 - s2_group1 - s3_group1 
  
  # Simulating and calculating the statistic test RISK DIFFERENCE
  # Test<-data.frame()
  data$Test_Reject[i] <- sum(replicate(nsim,f(data$samplesizeP[i],s1_group0,s2_group0,s3_group0,s1_group1,s2_group1,s3_group1))< - z.alpha)/nsim
  
  # Simulating and calculating the statistic test RISK RATIO 
  data$Test_Reject_RR[i] <- sum(replicate(nsim,fRR(data$samplesizeP_RR[i],s1_group0,s2_group0,s3_group0,s1_group1,s2_group1,s3_group1))< - z.alpha)/nsim
  
  # Simulating and calculating the statistic test ODDS RATIO
  data$Test_Reject_OR[i] <- sum(replicate(nsim,fOR(data$samplesizeP_OR[i],s1_group0,s2_group0,s3_group0,s1_group1,s2_group1,s3_group1))< - z.alpha)/nsim
  
  
  
  # cat(i, "\n", file="LOG.txt", append=TRUE)
  if(i%%100==0){cat(i, "\n", file="LOG.txt", append=TRUE)}
}


t1=Sys.time()-t0
cat(t1, "\n", file="LOG.txt", append=TRUE)

#####################################################################################

rm(alpha,beta,corr,i,nsim,z.alpha,z.beta)
rm(p1_group0,p1_group1,p2_group0,p2_group1,phi_group0,phi_group1,s1_group0,s1_group1,s2_group0,s2_group1,s3_group0,s3_group1,s4_group0,s4_group1)

write.csv2(data,'RESULTS_H0True_CompletePooled.csv2')
(t1)

sfStop()



