################################################################
# CASE STUDY
# CompARE paper - 2019
# Marta Bofill 
################################################################


################################################################
# PREAMBLE
################################################################ 
library(ggplot2)


# FUNCTION: Odds Ratio
################################################################ 
# pk: probability of the k-th endpoint (k=1,2) 

OR.function<- function(p0, p1){
  
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}


# FUNCTION: Bahadur.composite
################################################################ 
# pk: probability of the k-th endpoint (k=1,2)
# rho: correlation between the endpoints

Bahadur.composite<- function(p1, p2, rho){
  
  p.composite<- 1- (1-p1)*(1-p2)*( 1+ rho*sqrt(p1*p2/((1-p1)*(1-p2)) ))
  return(p.composite)
}


# FUNCTION: Correlation bounds
################################################################ 
# pk: probability of the k-th endpoint (k=1,2)

correlation.min.function = function(p1, p2){
  rho.min <- max(  -sqrt(p1*p2/((1-p1)*(1-p2))), -sqrt(((1-p1)*(1-p2))/(p1*p2)  ) )  
  return(rho.min)
}

correlation.max.function = function(p1, p2){
  rho.max <- min(  sqrt( ( p1/ (1-p1) )/( p2/ (1-p2) ) ), sqrt((p2/(1-p2))/(p1/(1-p1)) ) )  
  return(rho.max)
}

# FUNCTION: OR.composite
################################################################ 
# pk.0: probability of the k-th endpoint under control group (k=1,2)
# ORk: Odds ratio for the k-th endpoint
# rho: correlation between the endpoints

OR.composite.function = function(p1.0, p2.0, OR1, OR2, rho){
  
  O10= p1.0/(1-p1.0)
  O20= p2.0/(1-p2.0)
  
  OR= ((O10*OR1+1)*(O20*OR2+1)-1-rho*sqrt(OR1*OR2*O10*O20))*(1+rho*sqrt(O10*O20))/
    (((1+O10)*(1+O20)-1-rho*sqrt(O10*O20))*(1+rho*sqrt(OR1*OR2*O10*O20))) 
  return(OR)
}

# FUNCTION: ARE.betaOR.binary.endpoints (beta parametrization) 
################################################################ 
# pk.0: probability of the k-th endpoint under control group (k=1,2)
# ORk: Odds ratio for the k-th endpoint
# rho: correlation between the endpoints

ARE.betaOR.binary.endpoints = function(p1.0, p2.0, OR1, OR2, rho){
  
  OR.composite <- OR.composite.function(p1.0, p2.0, OR1, OR2, rho)
  
  p.composite.0 <- Bahadur.composite(p1.0, p2.0, rho)
  
  are <- log(OR.composite)^2*p.composite.0*(1- p.composite.0)/(log(OR1)^2*p1.0*(1-p1.0))
  return(are)
}

# FUNCTION: ARE.p.binary.endpoints (p parametrization)
################################################################ 

ARE.p.binary.endpoints <- function(p1.0, p1.1, p2.0, p2.1, rho){
  
  p.composite.0 <- Bahadur.composite(p1.0, p2.0, rho)
  p.composite.1 <- Bahadur.composite(p1.1, p2.1, rho)
  
  are <- ( (p.composite.0 - p.composite.1)^2*(p1.0*(1-p1.0)) )/( (p1.0-p1.1)^2*(p.composite.0*(1-p.composite.0))  )
  return(are)
}

# FUNCTION: Sample Size Composite Binary Endpoint Diff
################################################################ 
SampleSize.CBE.Diff <- function(p1.0, p2.0, d1, d2, rho, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){ 
  
  p1.1= d1+p1.0
  p2.1 = d2+p2.0
  
  p0.CBE = 1- (1-p1.0)*(1-p2.0)*( 1+ rho*sqrt(p1.0*p2.0/((1-p1.0)*(1-p2.0)) )) 
  d.CBE <- diff.pcomp(p1.0, p2.0, p1.1, p2.1, rho) 
  
  n = SampleSize.Diff(p0.CBE,d.CBE,alpha,beta,Unpooled) 
  
  return (n)
}  

# FUNCTION: Sample Size Binary Diff
#####################################################################################
# Sample size calculations in terms of risk difference 

SampleSize.Diff <- function(p0, d, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  if(Unpooled=="Unpooled Variance"){
    # sample size per group
    n1 = ((z.alpha+z.beta)/d)^2*( p0*(1-p0) + (d+p0)*(1-p0-d))
  }else{
    p1 = p0 + d
    p = (p1 + p0)/2
    # sample size per group
    n1 = ((z.alpha* sqrt(2*p*(1-p)) +  z.beta* sqrt( p0*(1-p0) + (d+p0)*(1-p0-d)))/d)^2
  }
  n = 2*n1      
  
  return(n)
}


# FUNCTION: DIFFERENCES IN PROPORTIONS COMPOSITE
#####################################################################################
# pk: probability of the k-th endpoint (k=1,2)
# rho: correlation between the endpoints
diff.pcomp <- function(p0.RE, p0.AE, p1.RE, p1.AE, rho){ 
  diff = Bahadur.composite(p1.RE, p1.AE, rho)-Bahadur.composite(p0.RE, p0.AE, rho) 
  
  return(diff)
}


#####################################################################################
#####################################################################################
# E1: Target-vessel revascularization
p10 = 70/1195
p11 = 94/2416

d1 = p11-p10
round(d1,3)
OR1 = OR.function(p10,p11)

c(p10,p11,d1,OR1)

# E2: Cardiac death or target-vessel myocardial infarction
p20 = 38/1195
p21 = 53/2416

d2 = p21-p20
round(d2,3)
OR2 = OR.function(p20,p21)

c(p20,p21,d2,OR2)

# correlation 
rho.min= max(correlation.min.function(p10,p20),correlation.min.function(p11,p21)) 
rho.max= min(correlation.max.function(p10,p20),correlation.max.function(p11,p21))

#####################################################################################

# COMPOSITE UNDER CONTROL GROUP
rho0.min= correlation.min.function(p10,p20)
rho0.max= correlation.max.function(p10,p20)

c(rho0.min,rho0.max)

Bahadur.composite(p10,p20,rho0.min)
Bahadur.composite(p10,p20,rho0.max)

# ARE in OR

f <- function(x,p10, p20, OR1, OR2) ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho=x)
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))

ggplot(data.frame(x=c(0, 1)), aes(x)) + stat_function(fun=f,args=list(p10=p10, p20=p20, OR1=OR1, OR2=OR2),size=1.3)+  
  geom_hline(yintercept = 1,col='black',linetype=2) + ylim(ymin, ymax)  + xlab('Correlation') + ylab('ARE')


# ARE in p
fp <- function(x,p10, p11, p20, p21) ARE.p.binary.endpoints(p10, p11, p20, p21, rho=x)
# ARE.p.binary.endpoints(p1.0, p1.1, p2.0, p2.1, rho)
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))

ggplot(data.frame(x=c(rho.min, rho.max)), aes(x)) + stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=p21),size=1.3)+  
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=0.020),size=1.3,col=2)+ 
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=0.025),size=1.3,col=3)+  
  geom_hline(yintercept = 1,col='black',linetype=2)  + ylim(ymin, ymax)  + xlab('Correlation') + ylab('ARE') 
# +guides(col=guide_legend(title="p21 Endpoint 2")) 

#####################################################################################

# COMPOSITE UNDER CONTROL GROUP


round(c(rho.min,(rho.max-rho.min)/2,rho.max),2 )

round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=rho.min),2)
round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=(rho.max-rho.min)/2),2)
round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=rho.max),2)

round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, rho.min))
round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, (rho.max-rho.min)/2))
round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, rho.max))

round(Bahadur.composite(p10,p20,rho.min),2)
round(Bahadur.composite(p10,p20,(rho.max-rho.min)/2),2)
round(Bahadur.composite(p10,p20,rho.max),2)

