################################################################
# CASE STUDY
# CompARE paper - 2019
# Marta Bofill 
################################################################

################################################################
# PREAMBLE
################################################################ 
setwd("C:/Users/marta.bofill/Dropbox/C5/Scripts/2017 Shiny/CompARE/BG_CompAREpaper")
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

# ggplot(data.frame(x=c(0, 1)), aes(x)) + stat_function(fun=f,args=list(p10=p10, p20=p20, OR1=OR1, OR2=OR2),size=1.3)+  
#   geom_hline(yintercept = 1,col='black',linetype=2) + ylim(ymin, ymax)  + xlab('Correlation') + ylab('ARE')


# ARE in p
fp <- function(x,p10, p11, p20, p21) ARE.p.binary.endpoints(p10, p11, p20, p21, rho=x) 
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))


windows(width = 9, height = 6)
ggplot(data.frame(x=c(rho.min, rho.max)), aes(x)) + 
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=p21),size=1.3)+  
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=0.020),size=1.3,col=2)+ 
  stat_function(fun=fp,args=list(p10=p10, p20=p20, p11=p11, p21=0.025),size=1.3,col=3)+  
  geom_hline(yintercept = 1,col='black',linetype=2) +  
  scale_y_log10(breaks=c(0.1,0.2,0.5,0.65,0.8,1,1.25,1.5,2,5,10),limits=c(ymin,ymax)) +scale_x_continuous(limits=c(0,1))+
  guides(col=guide_legend(title="Correlation")) + xlab('Correlation') + ylab('ARE') +
  theme(legend.position="bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold"))   

round(c(p21-p20,0.020-p20,0.025-p20),3)

################################################################ 

fp <- function(x,p10, p11, p20, p21) ARE.p.binary.endpoints(p10, p11, p20, p21, rho=x) 

RHO <- seq(0,1,0.01)
d <- data.frame(rho=rep(RHO,3),
                p10=p10,p20=p20,p11=p11,p21=rep(c(p21,0.020,0.025),each=length(RHO)))
d$f <- with(d,fp(rho,p10, p11, p20, p21))
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))

windows(width = 9, height = 6)
ggplot(d, aes(x=rho,y=f,col=as.factor(round(p21,3)))) + 
  geom_line(size=1.3)+  
  scale_y_log10(breaks=c(0.1,0.2,0.5,0.65,0.8,1,1.25,1.5,2,5,10),limits=c(ymin,ymax)) +scale_x_continuous(limits=c(0,1))+
  geom_hline(yintercept = 1,col='black',linetype=2) +
  guides(col=guide_legend(title=expression(bold(Probability~E[2]~"treated group")))) + xlab('Correlation') + ylab('ARE') +
  theme(legend.position="bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold"))

################################################################ 

fp <- function(x,p10, p11, p20, p21) ARE.p.binary.endpoints(p10, p11, p20, p21, rho=x) 

RHO <- seq(0,1,0.01)
d <- data.frame(rho=rep(RHO,3),
                p10=p10,p20=p20,p11=p11,p21=rep(c(p21,0.020,0.025),each=length(RHO)))
d$f <- with(d,fp(rho,p10, p11, p20, p21))
d$RD <- with(d,p20-p21)
ymax <- min(max(2,ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.min)+0.5), 100)
ymin <- min(0.5, max(0.005,min(ARE.betaOR.binary.endpoints(p10, p20, OR1, OR2, rho0.max) -0.15)))

windows(width = 9, height = 6)
ggplot(d, aes(x=rho,y=f,col=factor(formatC(RD,digits=3,format='f'),levels=c('0.012','0.010','0.007')))) + 
  geom_line(size=1.3)+  
  scale_y_log10(breaks=c(0.1,0.2,0.5,0.65,0.8,1,1.25,1.5,2,5,10),limits=c(ymin,ymax)) +scale_x_continuous(limits=c(0,1))+
  geom_hline(yintercept = 1,col='black',linetype=2) +
  guides(col=guide_legend(title=expression(bold("Absolute Reduction"~E[2])))) + xlab('Correlation') + ylab('ARE') +
  theme(legend.position="bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold"))
################################################################ 
# TABLE PROBABILITIES UNDER CONTROL GROUP / ARE / SAMPLE SIZE 
################################################################ 
round(c(rho.min,(rho.max-rho.min)/2,rho.max),2 )


################################################################ 
round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=0),2)
round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=0.1),2)
round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=0.4),2)
round(ARE.p.binary.endpoints(p10, p11, p20, p21, rho=0.7),2)

round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, 0))
round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, 0.1))
round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, 0.4))
round(SampleSize.CBE.Diff(p10, p20, p11-p10, p21-p20, 0.7))

round(Bahadur.composite(p10,p20,0),2)
round(Bahadur.composite(p10,p20,0.1),2)
round(Bahadur.composite(p10,p20,0.4),2)
round(Bahadur.composite(p10,p20,0.7),2)

round(cond.prob(p10,p20,0),2)
round(cond.prob(p10,p20,0.1),2)
round(cond.prob(p10,p20,0.4),2)
round(cond.prob(p10,p20,0.7),2) 
# round(cond.prob(p10,p20,rho.max),2) 

round(cond.prob(p20,p10,0),2)
round(cond.prob(p20,p10,0.1),2)
round(cond.prob(p20,p10,0.4),2)
round(cond.prob(p20,p10,0.7),2) 

round(diff.pcomp(p10, p20, p11, p21, 0),3)
round(diff.pcomp(p10, p20, p11, p21, 0.1),3)
round(diff.pcomp(p10, p20, p11, p21, 0.4),3)
round(diff.pcomp(p10, p20, p11, p21, 0.7),3)

corr.composite(p10, p20, 0.084) 
 

