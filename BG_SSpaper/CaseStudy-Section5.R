################################################################
# CASE STUDY
# Statistics in Medicine - 2018 Oct
# Marta Bofill and Guadalupe Gómez
################################################################


################################################################
# PREAMBLE
################################################################ 


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



# FUNCTION: Loss in power Diff CBE
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


# FUNCTION: Loss in power Diff
#####################################################################################   

Power.Diff <- function(p0, d, n, alpha=0.05){
  z.alpha <- qnorm(1-alpha,0,1)    
  
  z = sqrt(n/2)*abs(d)/sqrt(p0*(1-p0)+(p0+d)*(1-p0-d)) - z.alpha
  
  power = pnorm(z)
  
  return (power)
}  

#####################################################################################



################################################################
# PARAMETERS OF THE PROBLEM
################################################################

beta=0.2

# ENDPOINT 1
# Death or nonfatal myocardial infarction
p1.0 = 0.095
p1.1 = 0.073
d1=p1.1-p1.0 

# CI
alpha=0.05/2
z.alpha <- qnorm(1-alpha,0,1)  
p = 0.095
n=1106
v1 = z.alpha*sqrt(p*(1-p)/n)

# ENDPOINT 2
# Rehospitalization for acute coronary syndrome  
p2.0 = 0.137 
p2.1 = 0.110
d2=p2.1-p2.0 

# CI
alpha=0.05/2
z.alpha <- qnorm(1-alpha,0,1) 
p = 0.137
n=1106
v2 = z.alpha*sqrt(p*(1-p)/n)

# INTERVALS
round(c(p1.0-v1,p1.0+v1),3) 
round(c(p2.0-v2,p2.0+v2),3)
#####################################################################################   


################################################################
# MAIN DOCUMENT
################################################################

# Table 5 MANUSCRIPT

################################################################
# CORRELATION BOUNDS  based on p1,p2 and d1,d2
################################################################
rho0.min <- correlation.min.function(p1.0, p2.0)
rho0.max <- correlation.max.function(p1.0, p2.0) 

rho1.min <- correlation.min.function(p1.1, p2.1)
rho1.max <- correlation.max.function(p1.1, p2.1)

rho.min <- max(rho0.min, rho1.min)
rho.max <- min(rho0.max, rho1.max)

round(c(rho.min, rho.max),2)


################################################################
# SAMPLE SIZE based on p1,p2 and d1,d2
################################################################

r.weak = rho.min + (rho.max-rho.min)/3
r.mod = rho.min + 2*(rho.max-rho.min)/3
r.str = rho.max

n.weak = SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, r.weak, alpha = 0.025,Unpooled="Pooled Variance")  
n.mod = SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, r.mod, alpha = 0.025,Unpooled="Pooled Variance") 
n.str = SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, r.str, alpha = 0.025,Unpooled="Pooled Variance") 

c(round(r.weak,2), round(n.weak))
c(round(r.mod,2), round(n.mod))
c(round(r.str, 2), round(n.str))


################################################################
# POWER based on p1,p2 and d1,d2
################################################################

# range of corr in the approach 1: rho.min to r.weak
round(Power.CBE.Diff(p1.0, p2.0, d1, d2, rho.min, n.weak, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0, p2.0, d1, d2, r.weak, n.weak, alpha = 0.025,Unpooled="Pooled Variance"),2)

# range of corr in the approach 1: r.weak to r.mod
round(Power.CBE.Diff(p1.0, p2.0, d1, d2, r.weak, n.mod, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0, p2.0, d1, d2, r.mod, n.mod, alpha = 0.025,Unpooled="Pooled Variance"),2) 

# range of corr in the approach 1: r.mod to r.str
round(Power.CBE.Diff(p1.0, p2.0, d1, d2, r.mod, n.str, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0, p2.0, d1, d2, r.str, n.str, alpha = 0.025,Unpooled="Pooled Variance"),2) 



################################################################
# CORRELATION BOUNDS example
################################################################

round(
  c(max(correlation.min.function(p1.0, p2.0),correlation.min.function(p1.1, p2.1)),
    min(correlation.max.function(p1.0, p2.0), correlation.max.function(p1.1, p2.1))),2)

round(
  c(max(correlation.min.function(p1.0-v1, p2.0-v2),correlation.min.function(d1+p1.0-v1, d2+p2.0-v2)),
    min(correlation.max.function(p1.0-v1, p2.0-v2),correlation.max.function(d1+p1.0-v1, d2+p2.0-v2) )),2 )

round(
  c(max(correlation.min.function(p1.0+v1, p2.0+v2),correlation.min.function(d1+p1.0+v1, d2+p2.0+v2)),
    min(correlation.max.function(p1.0+v1, p2.0+v2),correlation.max.function(d1+p1.0+v1, d2+p2.0+v2))) ,2)




################################################################
# CORRELATION BOUNDS based on I1,I2 and d1,d2
################################################################
rho0.min <- correlation.min.function(p1.0, p2.0)
rho0.max <- correlation.max.function(p1.0, p2.0) 

rho0.minv1 <- correlation.min.function(p1.0+v1, p2.0+v2)
rho0.minv2 <- correlation.min.function(p1.0-v1, p2.0-v2) 

rho0.maxv1 <- correlation.max.function(p1.0+v1, p2.0+v2)
rho0.maxv2 <- correlation.max.function(p1.0-v1, p2.0-v2) 

rho1.min <- correlation.min.function(d1+p1.0, d2+p2.0)
rho1.max <- correlation.max.function(d1+p1.0, d2+p2.0) 

rho1.minv1 <- correlation.min.function(d1+p1.0+v1, d2+p2.0+v2)
rho1.minv2 <- correlation.min.function(d1+p1.0-v1, d2+p2.0-v2) 

rho1.maxv1 <- correlation.max.function(d1+p1.0+v1, d2+p2.0+v2)
rho1.maxv2 <- correlation.max.function(d1+p1.0-v1, d2+p2.0-v2) 

rho.minI <- max(rho0.min, rho0.minv1, rho0.minv2, rho1.min, rho1.minv1, rho1.minv2)
rho.maxI <- min(rho0.max, rho0.maxv1, rho0.maxv2, rho1.max, rho1.maxv1, rho1.maxv2)

round(c(rho.minI, rho.maxI),2)

################################################################
# SAMPLE SIZE based on I1,I2 and d1,d2
################################################################

r.weakI = rho.minI + (rho.maxI-rho.minI)/3
r.modI = rho.minI + 2*(rho.maxI-rho.minI)/3
r.strI = rho.maxI

n.weakI = SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.weakI, alpha = 0.025,Unpooled="Pooled Variance")  
n.modI = SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.modI, alpha = 0.025,Unpooled="Pooled Variance") 
n.strI = SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.strI, alpha = 0.025,Unpooled="Pooled Variance") 

c(round(r.weakI,2), round(n.weakI))
c(round(r.modI,2), round(n.modI))
c(round(r.strI, 2), round(n.strI))


################################################################
# EXAMPLE text
round(SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, 0.3, alpha = 0.025,Unpooled="Pooled Variance"))
round(SampleSize.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, 0.3, alpha = 0.025,Unpooled="Pooled Variance"))
round(SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, 0.3, alpha = 0.025,Unpooled="Pooled Variance"))

 

################################################################
# POWER based on I1,I2 and d1,d2
################################################################

# range of corr in the approach 1: rho.min to r.weak
round(Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.weakI, n.weakI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, rho.minI, n.weakI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, r.weakI, n.weakI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, rho.minI, n.weakI, alpha = 0.025,Unpooled="Pooled Variance"),2)

# range of corr in the approach 1: r.weak to r.mod
round(Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.modI, n.modI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, r.weakI, n.modI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.weakI, n.modI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, r.modI, n.modI, alpha = 0.025,Unpooled="Pooled Variance"),2)

# range of corr in the approach 1: r.mod to r.str
round(Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.strI, n.strI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.modI, n.strI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, r.strI, n.strI, alpha = 0.025,Unpooled="Pooled Variance"),2)
round(Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, r.modI, n.strI, alpha = 0.025,Unpooled="Pooled Variance"),2)




##################################################################################### 
# SAMPLE SIZE  AND POWER PLOTS
##################################################################################### 
# using polygon

windows(width=10, height=5)
par(mfrow=c(1,2))

# SAMPLE SIZE
plot(NA,
     xlab = expression(paste("Correlation ", "(", rho, ")")),
     # xlab = expression(paste("Correlation ", "(", rho["true"], ")")),
     ylab = paste("Sample Size"),
     xlim=c(rho.min,rho.max),  
     ylim=c(800, 5000), 
     n = 1001, col = 1, lty=1, lwd=2, yaxt="n") 
x=c(1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)
axis(2, at=x,labels=x,   las=2, cex.axis = 0.8)

###################################
# AREAS
x <- seq(r.modI,r.strI,0.001)
x2 <- seq(r.strI,r.modI,-0.001) 
y <- SampleSize.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, x, alpha = 0.025,Unpooled="Pooled Variance")
y2 <- SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, x2, alpha = 0.025,Unpooled="Pooled Variance")

X <- c(x,x2)
Y <- c(y,y2)
polygon(X,Y,col= "Pink 1 ", border = NA)


x <- seq(r.weakI,r.modI,0.001)
x2 <- seq(r.modI,r.weakI,-0.001) 
y <- SampleSize.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, x, alpha = 0.025,Unpooled="Pooled Variance")
y2 <- SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, x2, alpha = 0.025,Unpooled="Pooled Variance")

X <- c(x,x2)
Y <- c(y,y2)
polygon(X,Y,col= "Dark Sea Green 1", border = NA)


x <- seq(rho.minI,r.weakI,0.001)
x2 <- seq(r.weakI,rho.minI,-0.001)
y <- SampleSize.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, x, alpha = 0.025,Unpooled="Pooled Variance")
y2 <- SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, x2, alpha = 0.025,Unpooled="Pooled Variance")

X <- c(x,x2)
Y <- c(y,y2)
polygon(X,Y,col= "Light Steel Blue 1", border = NA)

###################################
curve(SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, x,alpha = 0.025,Unpooled="Pooled Variance"),
      from =rho.min, to = rho.max, lty=1, lwd=1, add=T)  

curve(SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, r.weak, alpha = 0.025,Unpooled="Pooled Variance")+0*x, from =rho.min,  to = r.weak, 
      ylab = paste("Sample Size"),
      xlim=c(rho.min,rho.max),  yaxt="n", 
      n = 1001, col = 4, lwd=2, add=T)
curve(SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, r.mod, alpha = 0.025,Unpooled="Pooled Variance")+0*x, from =r.weak, to =r.mod, 
      xlim=c(rho.min,rho.max),  yaxt="n", 
      n = 1001, col = 3, lwd=2, add=T)
curve(SampleSize.CBE.Diff(p1.0, p2.0, d1, d2, rho.max, alpha = 0.025,Unpooled="Pooled Variance")+0*x, from =r.mod, to = rho.max, 
      xlim=c(rho.min,rho.max), yaxt="n", 
      n = 1001, col = 2, lwd=2, add=T)

curve(SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.weak, alpha = 0.025,Unpooled="Pooled Variance")+0*x, from =rho.minI,  to = r.weakI, 
      ylab = paste("Sample Size"),
      xlim=c(rho.min,rho.max),  yaxt="n", 
      n = 1001, lty=2, col = 4, lwd=2, add=T)
curve(SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, r.mod, alpha = 0.025,Unpooled="Pooled Variance")+0*x, from =r.weakI, to =r.modI, 
      xlim=c(rho.min,rho.max),  yaxt="n", 
      n = 1001, lty=2, col = 3, lwd=2, add=T)
curve(SampleSize.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, rho.max, alpha = 0.025,Unpooled="Pooled Variance")+0*x, from =r.modI, to = rho.maxI, 
      xlim=c(rho.min,rho.max), yaxt="n", 
      n = 1001, lty=2, col = 2, lwd=2, add=T)

legend("topleft", legend = c("Weak", "Moderate", "Strong"), 
       col=c(4,3,2), ncol=2 , lwd=2, cex=0.7, pt.cex=0.9) 


###################################
# POWER

plot(NA,
     xlab = expression(paste("Correlation ", "(", rho, ")")),
     # xlab = expression(paste("Correlation ", "(", rho["true"], ")")),
     ylab = paste("Power"),
     xlim=c(rho.min,rho.max), ylim=c(0.6,1), 
     n = 1001, col = 2, lty=1, lwd=2, yaxt="n") 

x=c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
axis(2, at=x,labels=x,   las=2, cex.axis = 0.8) 

###################################
# AREAS
x <- seq(r.modI,rho.maxI,0.001)
x2 <- seq(rho.maxI,r.modI,-0.001) 
y <- Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, x, n.strI,alpha = 0.025,Unpooled="Pooled Variance")
y2 <- Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, x2, n.strI,alpha = 0.025,Unpooled="Pooled Variance")

X <- c(x,x2)
Y <- c(y,y2)
polygon(X,Y,col= "Pink 1 ", border = NA)

x <- seq(r.weakI,r.modI,0.001)
x2 <- seq(r.modI,r.weakI,-0.001) 
y <- Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, x, n.modI,alpha = 0.025,Unpooled="Pooled Variance")
y2 <- Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, x2, n.modI,alpha = 0.025,Unpooled="Pooled Variance")
X <- c(x,x2)
Y <- c(y,y2)
polygon(X,Y,col= "Dark Sea Green 1", border = NA)

x <- seq(rho.minI,r.weakI,0.001)
x2 <- seq(r.weakI,rho.minI,-0.001) 
y <- Power.CBE.Diff(p1.0-v1, p2.0-v2, d1, d2, x, n.weakI,alpha = 0.025,Unpooled="Pooled Variance")
y2 <- Power.CBE.Diff(p1.0+v1, p2.0+v2, d1, d2, x2, n.weakI,alpha = 0.025,Unpooled="Pooled Variance")
X <- c(x,x2)
Y <- c(y,y2)
polygon(X,Y,col= "Light Steel Blue 1", border = NA)

###################################
curve(Power.CBE.Diff(p1.0, p2.0, d1, d2, x, n.str,alpha = 0.025,Unpooled="Pooled Variance" ), from =r.mod, to =rho.max, add=T, lwd=2, col=2)
curve(Power.CBE.Diff(p1.0, p2.0, d1, d2, x, n.mod,alpha = 0.025,Unpooled="Pooled Variance" ), from =r.weak, to =r.mod, add=T, lwd=2, col=3)
curve(Power.CBE.Diff(p1.0, p2.0, d1, d2, x, n.weak,alpha = 0.025,Unpooled="Pooled Variance" ), from =rho.min,  to = r.weak,  add=T, lwd=2, col=4)


abline(h=0.8,   lty=2)  
legend("topright", legend = c("Weak", "Moderate", "Strong"), 
       col=c(4,3,2), ncol=2 , lwd=2, cex=0.7, pt.cex=0.9)
###################################



