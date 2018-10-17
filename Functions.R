

# FUNCTION: Odds Ratio
#####################################################################################
# pk: probability of the k-th endpoint (k=1,2) 

OR.function<- function(p0, p1){
  
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
}
  

# FUNCTION: Bahadur.composite
#####################################################################################
# pk: probability of the k-th endpoint (k=1,2)
# rho: correlation between the endpoints

Bahadur.composite<- function(p1, p2, rho){
  
  p.composite<- 1- (1-p1)*(1-p2)*( 1+ rho*sqrt(p1*p2/((1-p1)*(1-p2)) ))
  return(p.composite)
}

# FUNCTION: DIFFERENCES IN PROPORTIONS COMPOSITE
#####################################################################################
# pk: probability of the k-th endpoint (k=1,2)
# rho: correlation between the endpoints
diff.pcomp <- function(p0.RE, p0.AE, p1.RE, p1.AE, rho){ 
  diff = Bahadur.composite(p1.RE, p1.AE, rho)-Bahadur.composite(p0.RE, p0.AE, rho) 
  
  return(diff)
}


# FUNCTION: RISK RATIO COMPOSITE
#####################################################################################
# pk: probability of the k-th endpoint (k=1,2)
# rho: correlation between the endpoints
RR.pcomp <- function(p0.RE, p0.AE, p1.RE, p1.AE, rho){ 
  RR = Bahadur.composite(p1.RE, p1.AE, rho)/Bahadur.composite(p0.RE, p0.AE, rho)
  
  return(RR)
}

# FUNCTION: OR.composite
#####################################################################################
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
############################################################
# pk.0: probability of the k-th endpoint under control group (k=1,2)
# ORk: Odds ratio for the k-th endpoint
# rho: correlation between the endpoints

ARE.betaOR.binary.endpoints = function(p1.0, p2.0, OR1, OR2, rho){
  
  OR.composite <- OR.composite.function(p1.0, p2.0, OR1, OR2, rho)
  
  p.composite.0 <- Bahadur.composite(p1.0, p2.0, rho)
  
  are <- log(OR.composite)^2*p.composite.0*(1- p.composite.0)/(log(OR1)^2*p1.0*(1-p1.0))
  return(are)
}


# FUNCTION: Correlation bounds
############################################################
# pk: probability of the k-th endpoint (k=1,2)

correlation.min.function = function(p1, p2){
  rho.min <- max(  -sqrt(p1*p2/((1-p1)*(1-p2))), -sqrt(((1-p1)*(1-p2))/(p1*p2)  ) )  
  return(rho.min)
}

correlation.max.function = function(p1, p2){
  rho.max <- min(  sqrt( ( p1/ (1-p1) )/( p2/ (1-p2) ) ), sqrt((p2/(1-p2))/(p1/(1-p1)) ) )  
  return(rho.max)
}


# FUNCTION: Standardized Treatment effect CE OR
#####################################################################################

Standardized_TreatEff_OR <- function(p1.0, p2.0, OR1, OR2, rho){
  p0 = Bahadur.composite(p1.0, p2.0, rho)
  OR = OR.composite.function(p1.0, p2.0, OR1, OR2, rho)
  p1 = (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  
  Stad_TreatOR = (log(OR))/(sqrt( 1/(p0*(1-p0)) + 1/(p1*(1-p1)) ))
  
  return(Stad_TreatOR)
}


# FUNCTION: Sample Size Binary OR
#####################################################################################
# Sample size calculations in terms of odds ratio  

SampleSize.OR <- function(p0, OR, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  p1 = (OR*p0/(1-p0))/(1+(OR*p0/(1-p0)))
  
  if(Unpooled=="Unpooled Variance"){
    # sample size per group
    n1 = ((z.alpha+z.beta)/(log(OR)))^2*( 1/(p0*(1-p0)) + 1/(p1*(1-p1)) )  
  }else{
    p = (p1 + p0)/2
    # sample size per group
    n1 = ((z.alpha* sqrt(2/(p*(1-p))) + z.beta* sqrt(1/(p0*(1-p0)) + 1/(p1*(1-p1))))/(log(OR)))^2 
  }
   
  n = 2*n1  
  
  return(n)
}

# FUNCTION: Sample Size Composite Binary Endpoint OR
##################################################################################### 
SampleSize.CBE.OR <- function(p1.0, p2.0, OR1, OR2, rho, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){ 
  
  p1.1= (OR1*p1.0/(1-p1.0))/(1+(OR1*p1.0/(1-p1.0)))
  p2.1 = (OR2*p2.0/(1-p2.0))/(1+(OR2*p2.0/(1-p2.0)))
  
  p0.CBE = 1- (1-p1.0)*(1-p2.0)*( 1+ rho*sqrt(p1.0*p2.0/((1-p1.0)*(1-p2.0)) ))
  p1.CBE = 1- (1-p1.1)*(1-p2.1)*( 1+ rho*sqrt(p1.1*p2.1/((1-p1.1)*(1-p2.1)) ))
  
  OR.CBE = OR.composite.function(p1.0, p2.0, OR1, OR2, rho)
  
  n =  SampleSize.OR(p0.CBE, OR.CBE, alpha, beta, Unpooled)
  
  return (n)
}


# FUNCTION: Sample Size Binary R
#####################################################################################
# Sample size calculations in terms of risk ratio (assuming unpooled variance)

SampleSize.RR <- function(p0, R, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){
  z.alpha <- qnorm(1-alpha,0,1)  
  z.beta <-  qnorm(1-beta,0,1)
  
  
  if(Unpooled=="Unpooled Variance"){
    # sample size per group
    n1 = ((z.alpha+z.beta)/(log(R)))^2*( (1-R*p0)/(R*p0) + (1-p0)/p0 )
  }else{
    p1= R*p0
    p = (p1 + p0)/2
    
    # sample size per group
    n1 = ( (z.alpha* sqrt(2*(1-p)/p) + z.beta* sqrt( (1-R*p0)/(R*p0) + (1-p0)/p0) )/(log(R)) )^2  
  }
  
  n = 2*n1      
  
  return(n)
}

# FUNCTION: Sample Size Composite Binary Endpoint R
##################################################################################### 
SampleSize.CBE.RR <- function(p1.0, p2.0, R1, R2, rho, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){ 
  
  p1.1=  R1*p1.0 
  p2.1 = R2*p2.0 
  
  p0.CBE = 1- (1-p1.0)*(1-p2.0)*( 1+ rho*sqrt(p1.0*p2.0/((1-p1.0)*(1-p2.0)) )) 
  RR.CBE <- RR.pcomp(p1.0, p2.0, p1.1, p2.1, rho)
  
  n = SampleSize.RR(p0.CBE, RR.CBE, alpha, beta, Unpooled)  
  
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



# FUNCTION: Sample Size Composite Binary Endpoint Diff
##################################################################################### 
SampleSize.CBE.Diff <- function(p1.0, p2.0, d1, d2, rho, alpha=0.05, beta=0.2, Unpooled="Unpooled Variance"){ 
  
  p1.1= d1+p1.0
  p2.1 = d2+p2.0
  
  p0.CBE = 1- (1-p1.0)*(1-p2.0)*( 1+ rho*sqrt(p1.0*p2.0/((1-p1.0)*(1-p2.0)) )) 
  d.CBE <- diff.pcomp(p1.0, p2.0, p1.1, p2.1, rho) 
  
  n = SampleSize.Diff(p0.CBE,d.CBE,alpha,beta,Unpooled) 
  
  return (n)
}  




# FUNCTION: Loss in power Diff
#####################################################################################   

Power.CBE.Diff <- function(p1.0, p2.0, d1, d2, rho, n, alpha=0.05, Unpooled="Unpooled Variance"){
  z.alpha <- qnorm(1-alpha,0,1)  
  
  p1.1= d1+p1.0
  p2.1 = d2+p2.0
  
  p0.CBE = 1- (1-p1.0)*(1-p2.0)*( 1+ rho*sqrt(p1.0*p2.0/((1-p1.0)*(1-p2.0)) )) 
  d.CBE <- diff.pcomp(p1.0, p2.0, p1.1, p2.1, rho) 
  
  if(Unpooled=="Unpooled Variance"){ 
    z = sqrt(n/2)*abs(d.CBE)/sqrt(p0.CBE*(1-p0.CBE)+(p0.CBE+d.CBE)*(1-p0.CBE-d.CBE)) - z.alpha 
  }else{
    p1 = p0.CBE + d.CBE
    p = (p1 + p0.CBE)/2
    # sample size per group 
    z = ( abs(d.CBE)*sqrt(n/2) - z.alpha* sqrt(2*p*(1-p)) )/sqrt( p0.CBE*(1-p0.CBE) + (d.CBE+p0.CBE)*(1-p0.CBE-d.CBE))
  }
  
  
  power = pnorm(z)
  
  return (power)
}  

