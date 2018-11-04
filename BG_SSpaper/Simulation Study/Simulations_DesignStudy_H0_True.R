
################################################################
# SIMULATIONS SAMPLE SIZE FOR COMPOSITE BINARY ENDPOINTS
# Defining scenarios H0 true  
# oct-2018  
# Marta Bofill and Guadalupe Gómez
################################################################

setwd("C:/Users/marta.bofill/Simulations")

#####################################################################################
#####################################################################################
# FUNCTIONS
#####################################################################################


# FUNCTION: Bahadur.composite
#####################################################################################
# Description: 
# this function returns the probability of a composite binary endpoint given 
# the probabilities of its components and the correlation between them.
# 
# Arguments:
# pk: probability of the k-th endpoint (k=1,2)
# rho: correlation between the endpoints

Bahadur.composite<- function(p1, p2, rho){
  
  p.composite<- 1- (1-p1)*(1-p2)*( 1+ rho*sqrt(p1*p2/((1-p1)*(1-p2)) ))
  return(p.composite)
}

# FUNCTION: Correlation bounds
#####################################################################################
# Description: 
# these functions return the lower and upper bound of the correlation given 
# the probabilities of the single events.
# 
# Arguments:
# pk: probability of the k-th endpoint (k=1,2)

correlation.min.function<- function(p1, p2){
  rho.min <- max(  -sqrt(p1*p2/((1-p1)*(1-p2))), -sqrt(((1-p1)*(1-p2))/(p1*p2)  ) )  
  return(rho.min)
}

correlation.max.function <- function(p1, p2){
  rho.max <- min(  sqrt( ( p1/ (1-p1) )/( p2/ (1-p2) ) ), sqrt((p2/(1-p2))/(p1/(1-p1)) ) )  
  return(rho.max)
}

# FUNCTION: Sample Size Binary Diff
#####################################################################################
# Sample size calculations in terms of risk difference (assuming unpooled variance)

SampleSize.Diff <- function(p0, d, alpha=0.05, beta=0.2){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  # sample size per group
  n = ((z.alpha+z.beta)/d)^2*( p0*(1-p0) + (d+p0)*(1-p0-d))
  
  
  return(n)
}

# FUNCTION: Sample Size Binary Diff 
#####################################################################################
# Sample size calculations in terms of risk difference (assuming pooled variance)

SampleSizePooled.Diff <- function(p0, d, alpha=0.05, beta=0.2){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  p1 = p0 + d
  p = (p1 + p0)/2
  
  # sample size per group
  n = ((z.alpha* sqrt(2*p*(1-p)) +  z.beta* sqrt( p0*(1-p0) + (d+p0)*(1-p0-d)))/d)^2
  
  
  return(n)
}


# FUNCTION: Sample Size Binary R
#####################################################################################
# Sample size calculations in terms of risk ratio (assuming unpooled variance)

SampleSize.RR <- function(p0, R, alpha=0.05, beta=0.2){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  # sample size per group
  n = ((z.alpha+z.beta)/(log(R)))^2*( (1-R*p0)/(R*p0) + (1-p0)/p0 )    
  
  return(n)
}


# FUNCTION: Sample Size Binary R
#####################################################################################
# Sample size calculations in terms of risk ratio (assuming pooled variance)

SampleSizePooled.RR <- function(p0, R, alpha=0.05, beta=0.2){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  p1= R*p0
  p = (p1 + p0)/2
  
  # sample size per group
  n = ( (z.alpha* sqrt(2*(1-p)/p) + z.beta* sqrt( (1-R*p0)/(R*p0) + (1-p0)/p0) )/(log(R)) )^2    
  
  return(n)
}

# FUNCTION: Sample Size Binary OR
#####################################################################################
# Sample size calculations in terms of odds ratio (assuming unpooled variance)

SampleSize.OR <- function(p0, OR, alpha=0.05, beta=0.2){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  p1 = (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  
  # sample size per group
  n = ((z.alpha+z.beta)/(log(OR)))^2*( 1/(p0*(1-p0)) + 1/(p1*(1-p1)) )  
  
  return(n)
}


# FUNCTION: Sample Size Binary OR
#####################################################################################
# Sample size calculations in terms of odds ratio (assuming pooled variance)

SampleSizePooled.OR <- function(p0, OR, alpha=0.05, beta=0.2){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  p1 = (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  p = (p1 + p0)/2
  
  # sample size per group
  n = ((z.alpha* sqrt(2/(p*(1-p))) + z.beta* sqrt(1/(p0*(1-p0)) + 1/(p1*(1-p1))))/(log(OR)))^2 
  
  return(n)
}



#####################################################################################
#####################################################################################
# DESIGN OF SIMULATION STUDY
# Design
#####################################################################################

# Scenarios considered

p1.0 = c(0.01, 0.05, 0.10) 
p2.0 = c(0.01, 0.05, 0.10, 0.15, 0.20)

R1 = c(0.6, 0.7, 0.8)
R2 = c(0.6, 0.7, 0.8)

rho_true = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) 



# Defining the database 
# X0 =  num of the simulated scenario
# X0.1 = 'p1.0': probability of observing the relevant endpoint in the control group
# X0.2 = 'R1': risk ratio for the relevant endpoint
# X0.3 = 'p2.0': probability of observing the additional endpoint in the control group
# X0.4 = 'R2': risk ratio for the additional endpoint
# X0.5 = 'p1.1': probability of observing the relevant endpoint in the treatment group
# X0.6 = 'p2.1': probability of observing the additional endpoint in the treatment group
# X0.7 = 'rho.min': lower bound for the correlation given p1.0, p1.1, p2.0 and p2.1
# X0.8 = 'rho.max': upper bound for the correlation given p1.0, p1.1, p2.0 and p2.1
# X0.9 = 'rho_true': value of the correlation 
# X0.10 = 'rho_assumed': assumed correlation for the sample size calculations
# X0.11 = 'v': Sample size approach (0: true value correlation, 1: weak, 2:moderate, 3:strong).
# X0.12 = 'rho.rank': rho.max-rho.min.


library(plyr)
data<- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0)

data <- rename(data,c('c(0)'='Scenario',
                      X0 = 'num_scenario', 
                      X0.1 = 'num_scenario_theta',
                      X0.2 = 'p1.0',
                      X0.3 = 'R1',
                      X0.4 = 'p2.0',
                      X0.5 = 'R2',
                      X0.6 = 'p1.1',
                      X0.7 = 'p2.1',
                      X0.8 = 'rho.min',
                      X0.9 = 'rho.max', 
                      X0.10 = 'rho_true',
                      X0.11 = 'rho_assumed',
                      X0.12 = 'v' ))



# Creating the scenarios
i=1;j=1;k=1;l=1;m=1;it=1;
n=0;
epsilon=0.1;
num_scenario_theta=0;

for(i in 1:length(p1.0)){
  for(j in 1:length(R1)){
    
    p1.1= R1[j]*p1.0[i]
    
    for(k in 1:length(p2.0)){
      if(p1.0[i] < p2.0[k]){
        for(l in 1:length(R2)){
          
          p2.1= R2[l]*p2.0[k]
          
          # Rank(rho)
          rho0.min <- correlation.min.function(p1.0[i], p2.0[k])
          rho0.max <- correlation.max.function(p1.0[i], p2.0[k])
          
          rho1.min <- correlation.min.function(p1.1, p2.1)
          rho1.max <- correlation.max.function(p1.1, p2.1)
          
          rho.min <- max(rho0.min, rho1.min)
          rho.max <- min(rho0.max, rho1.max) 
          
          num_scenario_theta = num_scenario_theta+1;
          
          for(m in 1:length(rho_true)){
            if(rho_true[m] >=rho.min && rho_true[m] <=rho.max){
              n=n+1; 
              
              v=0
              rho_assumed = rho_true[m]
              data[it,]<- c(n, num_scenario_theta, p1.0[i], R1[j], p2.0[k], R2[l], p1.1, p2.1, rho.min, rho.max, rho_true[m], rho_assumed,  v)
              it=it+1;
              
              for(v in 1:3){
                rho_assumed = v*rho.max/3
                
                data[it,]<- c(n, num_scenario_theta, p1.0[i], R1[j], p2.0[k], R2[l], p1.1, p2.1, rho.min, rho.max, rho_true[m], rho_assumed, v)
                it=it+1;
                
              }
            }
          }
        }
      }
    } 
  }
}

rm(i,j,k,l,m,v,it,n,epsilon,num_scenario_theta)
rm(p1.0,p2.0,p1.1,p2.1,R1,R2,rho_true,rho0.min,rho0.max,rho1.min,rho1.max,rho.min,rho.max,rho_assumed)


# total number of scenarios 421
# total number of scenarios considering different rho_assumed 1684 (421*4)

data$rho.rank = data$rho.max - data$rho.min
data$error = data$rho_assumed - data$rho_true
data$indicator = (data$error>0)  


data$R1_true = 1
data$R2_true = 1


#####################################################################################
#####################################################################################
# DESIGN OF SIMULATION STUDY
# Sample size and parameters
#####################################################################################

alpha=0.025; beta=0.2
z.alpha <- qnorm(1-alpha,0,1)  
z.beta <-  qnorm(1-beta,0,1)



for(i in 1:dim(data)[1]){  
  data$pcomp.0[i] = Bahadur.composite(data$p1.0[i], data$p2.0[i], data$rho_assumed[i]);
  data$pcomp.1[i] = Bahadur.composite(data$p1.1[i], data$p2.1[i], data$rho_assumed[i]); 
  
  data$Rcomp[i] = data$pcomp.1[i]/data$pcomp.0[i]
  
  data$diffcomp[i] = data$pcomp.1[i]-data$pcomp.0[i]
  
  data$ORcomp[i] = (data$pcomp.1[i]*(1-data$pcomp.0[i]))/(data$pcomp.0[i]*(1-data$pcomp.1[i]))
  
  # sample size for risk difference using unpooled variance estimator
  data$samplesize[i] = round(SampleSize.Diff(data$pcomp.0[i], data$diffcomp[i], alpha=0.025))
  
  # sample size for risk ratio using unpooled variance estimator
  data$samplesize_RR[i] = round(SampleSize.RR(data$pcomp.0[i], data$Rcomp[i], alpha=0.025))
  
  # sample size for odds ratio using unpooled variance estimator
  data$samplesize_OR[i] = round(SampleSize.OR(data$pcomp.0[i], data$ORcomp[i], alpha=0.025))
  
  # sample size for risk difference using pooled variance estimator
  data$samplesizeP[i] = round(SampleSizePooled.Diff(data$pcomp.0[i], data$diffcomp[i], alpha=0.025))
  
  # sample size for risk ratio using pooled variance estimator
  data$samplesizeP_RR[i] = round(SampleSizePooled.RR(data$pcomp.0[i], data$Rcomp[i], alpha=0.025))
  
  # sample size for odds ratio using pooled variance estimator
  data$samplesizeP_OR[i] = round(SampleSizePooled.OR(data$pcomp.0[i], data$ORcomp[i], alpha=0.025))
  
}  

rm(i,alpha,beta,z.alpha,z.beta)


#####################################################################################

write.csv2(data,'DATABASE_Scenarios_H0_True.csv2')
# db<-read.table()

