################################################################
# CASE STUDY
# CompARE paper - 2019
# Marta Bofill 
################################################################

################################################################
# PREAMBLE
################################################################ 
library(ggplot2)
source('RFunctions.R')

################################################################ 
# VALUES FROM SPIRIT IV TRIAL
################################################################ 
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

################################################################ 
# Correlation Bounds
################################################################ 
# correlation under the control group
rho0.min= correlation.min.function(p10,p20)
rho0.max= correlation.max.function(p10,p20)

c(rho0.min,rho0.max)

# correlation 
rho.min= max(correlation.min.function(p10,p20),correlation.min.function(p11,p21)) 
rho.max= min(correlation.max.function(p10,p20),correlation.max.function(p11,p21))

c(rho.min,rho.max)

################################################################ 
# ARE
################################################################ 

# ARE in OR
f <- function(x,p10, p20, OR1, OR2) ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho=x)
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))

ggplot(data.frame(x=c(0, 1)), aes(x)) + stat_function(fun=f,args=list(p10=p10, p20=p20, OR1=OR1, OR2=OR2),size=1.3)+  
  geom_hline(yintercept = 1,col='black',linetype=2) + ylim(ymin, ymax)  + xlab('Correlation') + ylab('ARE')


# ARE in p
fp <- function(x,p10, p11, p20, p21) ARE.p.binary.endpoints(p10, p11, p20, p21, rho=x) 
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))

ggplot(data.frame(x=c(rho.min, rho.max)), aes(x)) + stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=p21),size=1.3)+  
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=0.020),size=1.3,col=2)+ 
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=0.025),size=1.3,col=3)+  
  geom_hline(yintercept = 1,col='black',linetype=2)  + ylim(ymin, ymax)  + xlab('Correlation') + ylab('ARE') 
# +guides(col=guide_legend(title="p21 Endpoint 2")) 


################################################################ 
# TABLE PROBABILITIES UNDER CONTROL GROUP / ARE / SAMPLE SIZE 
################################################################ 
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

