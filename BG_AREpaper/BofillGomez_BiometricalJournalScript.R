

#####################################################################################
# Title: Selection of composite binary endpoints in clinical trials 
# Authors: Marta Bofill Roig, Guadalupe Gómez Melis
# Journal: Biometrical Journal
# Last updated: 2017-08-18
# 
# Software version: R version 3.3.3 (2017-03-06)
# package 'plyr' version 1.8.4
#####################################################################################



#####################################################################################
#####################################################################################
# PART I: FUNCTIONS
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


# FUNCTION: OR.composite
#####################################################################################
# Description: 
# this function returns the odds ratio for a composite binary endpoint given 
# the probabilities of its components under control group, the odds ratio for the
# components of the composite endpoint and the correlation between them.
# 
# Arguments:
# pk.0: probability of the k-th endpoint under control group (k=1,2)
# ORk: Odds ratio for the k-th endpoint
# rho: correlation between the endpoints

OR.composite.function <- function(p1.0, p2.0, OR1, OR2, rho){
  
  O10= p1.0/(1-p1.0)
  O20= p2.0/(1-p2.0)
  
  OR= ((O10*OR1+1)*(O20*OR2+1)-1-rho*sqrt(OR1*OR2*O10*O20))*(1+rho*sqrt(O10*O20))/
    (((1+O10)*(1+O20)-1-rho*sqrt(O10*O20))*(1+rho*sqrt(OR1*OR2*O10*O20))) 
  return(OR)
}


# FUNCTION: OR
#####################################################################################
# Description: 
# this function returns the odds ratio given the probability of observing the event
# in the control and treatment group.
# 
# Arguments:
# p0: probability under control group
# p1: probability under test group

OR.function<- function(p0, p1){
  
  OR<- (p1/(1-p1))/(p0/(1-p0))
  return(OR)
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


# FUNCTION: ARE.betaOR.binary.endpoints 
#####################################################################################
# Description: 
# this function calculates the ARE method given the probabilities of its components 
# under control group, the odds ratio for the components of the composite endpoint 
# and the correlation between them.
# 
# Arguments:
# pk.0: probability of the k-th endpoint under control group (k=1,2)
# ORk: Odds ratio for the k-th endpoint
# rho: correlation between the endpoints

ARE.betaOR.binary.endpoints <- function(p1.0, p2.0, OR1, OR2, rho){
  
  OR.composite <- OR.composite.function(p1.0, p2.0, OR1, OR2, rho)
  
  p.composite.0 <- Bahadur.composite(p1.0, p2.0, rho)
  
  are <- log(OR.composite)^2*p.composite.0*(1- p.composite.0)/(log(OR1)^2*p1.0*(1-p1.0))
  return(are)
}


#####################################################################################
#####################################################################################
# PART II: ASYMPTOTIC RELATIVE EFFICIENCY
# Code to reproduce Figure 1 in Section 4
#####################################################################################

win.graph(width=8,height=6)
curve(dnorm(x,mean=-7.5,sd=1),from=-10,to=3,ylab="",xlab="",col=2,axes=FALSE)
curve(dnorm(x,mean=0,sd=1),from=-10,to=10,ylab="",xlab="",col=1,add=TRUE)
curve(dnorm(x,mean=-4.5,sd=1),from=-10,to=10,ylab="",xlab="",col=4,add=TRUE)
axis(1, xaxp=c(-12,3,10), fon=3,pos=-0.05)

segments(x0=0, y0=0, x1 = 0, y1 = 0.39,lty=2)
segments(x0=-4.5, y0=0, x1 = -4.5, y1 = 0.39,lty=2,col=4)
segments(x0=-7.5, y0=0, x1 = -7.5, y1 = 0.39,lty=2,col=2)

text(x=0,y=0.407,label="T under Ho")
text(x=-4.5,y=0.407,label="Under H1",col=4)
text(x=-7.5,y=0.407,label="Under H1",col=2)

text(x=-4.5,y=-0.009,label=expression(paste(delta,"*")),col=4)
text(x=-7.5,y=-0.009,label=expression(paste(delta,1)),col=2)

# savePlot("AN_behavior",type='pdf')
graphics.off()




#####################################################################################
#####################################################################################
# PART III: CASE STUDY
# Code to reproduce Figure 2 in Section 5
#####################################################################################

# Specifying the x axis values to plot 
fun <- function(x) 0*x+1

curve(fun(x),  log="y",  from= -0.2 , to =2, xlim=c(-0.2,0.7),  
      ylim=c(10^(-0.3), 10^(0.35)), 
      yaxt="n",
      col = 1,  ylab = "ARE(correlation)", xlab="Correlation")
 
x=c(0.5, 0.6, 0.7, 0.8, 1, 1.25, 1.5, 1.75, 2)
axis(2, at=x,labels=x,   las=2)

# Function ARE in terms of rho
ARE.rho <- function(rho){ 
  ARE.betaOR.binary.endpoints(p1.0, p2.0, OR1, OR2, rho)
}
# Function to plot the ARE curves in terms of rho
Plot.ARE.betaOR.rho <- function(p1.0, p1.1, p2.0, p2.1, i=2){
  # Rank(rho)
  rho0.min <- correlation.min.function(p1.0, p2.0)
  rho0.max <- correlation.max.function(p1.0, p2.0)  
  
  rho1.min <- correlation.min.function(p1.1, p2.1)
  rho1.max <- correlation.max.function(p1.1, p2.1)
  
  rho.min <- max(rho0.min, rho1.min)
  rho.max <- min(rho0.max, rho1.max) 
  
  
  OR1<- OR.function(p1.0, p1.1)
  OR2<- OR.function(p2.0, p2.1)
  
  ARE.rho <- function(rho){ 
    ARE.betaOR.binary.endpoints(p1.0, p2.0, OR1, OR2, rho)
  }
  # Plot ARE.rho vs rho
  curve(ARE.rho,  
        log="y",  
        from =rho.min, to = rho.max,  
        n = 1001, col = i, add = TRUE)
  
}



# Values in TAXUS V:
p1.0= 0.173;  p2.0= 0.055 ; 
p1.1=0.121; p2.1=0.057;


# Treatment effect in terms of Odds Ratio:
a= c(0.057, 0.050, 0.045, 0.040, 0.035)
OR1 = OR.function(p1.0, p1.1)
OR2= OR.function(p2.0, a)
OR2 = round(OR2,2)


Plot.ARE.betaOR.rho(p1.0, p1.1, p2.0, 0.057, 2)
Plot.ARE.betaOR.rho(p1.0, p1.1, p2.0, 0.050, 3)
Plot.ARE.betaOR.rho(p1.0, p1.1, p2.0, 0.045, 4)
Plot.ARE.betaOR.rho(p1.0, p1.1, p2.0, 0.040, 5)
Plot.ARE.betaOR.rho(p1.0, p1.1, p2.0, 0.035, 6)
par(xpd=F)
legend("topright", legend = OR2, col=c(2,3,4,5,6), ncol=2 , lwd=2, cex=0.6, inset=0.05, x.intersp=0.2, y.intersp=0.5)




#####################################################################################
#####################################################################################
# PART IV: STATISTICAL EFFICIENCY GUIDELINES 
#####################################################################################


# DESIGN
#####################################################################################

# Scenarios considered
p1.0 = c(0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100)
p2.0 = c(0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100)

OR1 = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)
OR2 = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

rho = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


# Defining the database
# 'c(0)'='Scenario': num scenario
# X0 = 'p1.0': probability of observing the relevant endpoint in the control group
# X0.1 = 'OR1': odds ratio for the relevant endpoint
# X0.2 = 'p2.0': probability of observing the additional endpoint in the control group
# X0.3 = 'OR2': odds ratio for the additional endpoint
# X0.4 = 'p1.1': probability of observing the relevant endpoint in the treatment group
# X0.5 = 'p2.1': probability of observing the additional endpoint in the treatment group
# X0.6 = 'rho.min': lower bound for the correlation given p1.0, p1.1, p2.0 and p2.1
# X0.7 = 'rho.max': upper bound for the correlation given p1.0, p1.1, p2.0 and p2.1
# X0.8 = 'rho': value of the correlation
# X0.9 = 'indicator': if indicator takes value 1, then rho is between rho.min and rho.max. Otherwise takes value 0.
# X0.11 = 'ARE.beta': computation of the ARE
# X0.12 = 'ORcomp': odds ratio for the composite endpoint given p1.0, p1.1, p2.0, p2.1 and rho
# X0.13 = 'pcomp.0': probability of observing the composite endpoint in the control group  given the values of p1.0, p1.1, p2.0, p2.1 and rho
# X0.14 = 'pcomp.1': probability of observing the composite endpoint in the treatment group given the values of p1.0, p1.1, p2.0, p2.1 and rho


library(plyr)
data<- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0)

data <- rename(data,c('c(0)'='Scenario',
                      X0 = 'p1.0',
                      X0.1 = 'OR1',
                      X0.2 = 'p2.0',
                      X0.3 = 'OR2',
                      X0.4 = 'p1.1',
                      X0.5 = 'p2.1',
                      X0.6 = 'rho.min',
                      X0.7 = 'rho.max',
                      X0.8 = 'rho',
                      X0.9 = 'indicator', 
                      X0.10 = 'ARE.beta',
                      X0.11 = 'ORcomp',
                      X0.12 = 'pcomp.0',
                      X0.13 = 'pcomp.1'))



# Creating the scenarios
i=1;j=1;k=1;l=1;m=1;it=1; n=0

# Start the clock
ptm <- proc.time()


for(i in 1:length(p1.0)){
  for(j in 1:length(OR1)){
    
    O10= p1.0[i]/(1-p1.0[i])
    p1.1= OR1[j]*O10/(1+OR1[j]*O10)
    
    for(k in 1:length(p2.0)){
      for(l in 1:length(OR2)){
        O20= p2.0[k]/(1-p2.0[k])
        p2.1= OR2[l]*O20/(1+OR2[l]*O20)  
        
        # Rank(rho)
        rho0.min <- correlation.min.function(p1.0[i], p2.0[k])
        rho0.max <- correlation.max.function(p1.0[i], p2.0[k])
        
        rho1.min <- correlation.min.function(p1.1, p2.1)
        rho1.max <- correlation.max.function(p1.1, p2.1)
        
        rho.min <- max(rho0.min, rho1.min)
        rho.max <- min(rho0.max, rho1.max)
        
        for(m in 1:length(rho)){
          if(rho[m] >=rho.min && rho[m] <=rho.max){
            indicator=1;
            n=n+1;
            
            pcomp.0 = Bahadur.composite(p1.0[i], p2.0[k], rho[m]);
            pcomp.1 = Bahadur.composite(p1.1, p2.1, rho[m]); 
            
            ORcomp = OR.composite.function(p1.0[i], p2.0[k], OR1[j], OR2[l], rho[m])  
            
            ARE.beta = ARE.betaOR.binary.endpoints(p1.0[i], p2.0[k], OR1[j], OR2[l], rho[m])
          }else{
            indicator=0;
            pcomp.0 = NA;
            pcomp.1 = NA;
            
            ORcomp = NA;
            ARE.p = NA;
          }
          
          data[it,]<- c(p1.0[i], OR1[j], p2.0[k], OR2[l], p1.1, p2.1, rho.min, rho.max, rho[m], indicator, 
                        ARE.beta, ORcomp, pcomp.0, pcomp.1)
          it=it+1
        }
      }
    } 
  }
}
(n)

# Stop the clock
elapsed<- proc.time() - ptm
time<-elapsed[3]

# Saving data
write.csv2(data,'DATABASE_Scenarios.csv2')
db<-read.table()



#####################################################################################
#####################################################################################
# GENERAL PATTERN OF THE PERCENTAGE OF CASES IN WHICH are>1 
# Section 6.2
#####################################################################################

# make sure the current working directory is the folder that contains DATABASE_Scenarios.csv 
# Reading the previous database
db<-read.table(file="./DATABASE_Scenarios.csv",sep = ";", header = T, dec=",")



# PATTERN %are>1 as a function of  RHO
# Data obtained correspond to Figure 3
#####################################################################################

data = subset(db, db$rho<1)
aux.rho = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ) 

guidelines.rho<- data.frame(0,0,0,0 )
guidelines.rho <- rename(guidelines.rho, c('c(0)'='Scenario',
                                           X0 = 'rho',
                                           X0.1 = 'total cases', 
                                           X0.2 = 'cases ARE.beta>1',
                                           X0.3 = '% ARE.beta>1'
)
)

for(i in 1:length(aux.rho)){
  # Number of cases rho= 0 and indicator=1
  ntotal = length(which(data$rho==aux.rho[i] & data$indicator==1)) 
  
  # Number cases ARE.beta>1
  n.arebeta = length(which(data$rho==aux.rho[i] & data$indicator==1 & data$ARE.beta>1)) 
  
  # % cases ARE.beta>1
  percent.arebeta = n.arebeta/ntotal
  
  guidelines.rho[i,]<- c(aux.rho[i], ntotal, n.arebeta, percent.arebeta)
}



# PATTERN %are>1 as a function of  OR1
# Data obtained correspond to Figure 5
#####################################################################################

auxOR = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)  

dataguidelines.OR1 <- data.frame(0,0,0,0)
dataguidelines.OR1 <- rename(dataguidelines.OR1, c('c(0)'='Scenario',
                                                   X0 = 'OR2',
                                                   X0.1 = 'total cases', 
                                                   X0.2 = 'cases ARE.beta>1',
                                                   X0.3 = '% ARE.beta>1'
)
)


for(i in 1:length(auxOR)){
  # Number of cases rho= 0 and indicator=1
  ntotal = length(which(data$OR1==auxOR[i] & data$indicator==1))
  
  # Number cases ARE.beta>1
  n.arebeta = length(which(data$OR1==auxOR[i] & data$indicator==1 & data$ARE.beta>1)) 
  
  # % cases ARE.beta>1
  percent.arebeta = n.arebeta/ntotal
  
  dataguidelines.OR1[i,]<- c(auxOR[i], ntotal, n.arebeta, percent.arebeta)
}




# PATTERN %are>1 as a function of OR2
# Data obtained correspond to Figure 6
#####################################################################################

guidelines.OR2 <- data.frame(0,0,0,0)
guidelines.OR2 <- rename(guidelines.OR2, c('c(0)'='Scenario',
                                           X0 = 'OR2',
                                           X0.1 = 'total cases', 
                                           X0.2 = 'cases ARE.beta>1',
                                           X0.3 = '% ARE.beta>1'
)
)


for(i in 1:length(auxOR)){
  # Number of cases rho= 0 and indicator=1
  ntotal = length(which(data$OR2==auxOR[i] & data$indicator==1)) 
  
  # Number cases ARE.beta>1
  n.arebeta = length(which(data$OR2==auxOR[i] & data$indicator==1 & data$ARE.beta>1)) 
  
  # % cases ARE.beta>1
  percent.arebeta = n.arebeta/ntotal
  
  guidelines.OR2[i,]<- c(auxOR[i], ntotal, n.arebeta, percent.arebeta)
}



# PATTERN %are>1 as a function of OR2, when OR1= 0.6
# Data obtained correspond to Figure 4
#####################################################################################

data1 <- subset(data, data$OR1 ==0.6)

dataguidelines.OR2.1 <- data.frame(0,0,0,0)
dataguidelines.OR2.1 <- rename(dataguidelines.OR2.1, c('c(0)'='Scenario',
                                                       X0 = 'OR2',
                                                       X0.1 = 'total cases',
                                                       X0.2 = 'cases ARE.beta>1',
                                                       X0.3 = '% ARE.beta>1'
)
)

auxOR = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95) 
for(i in 1:length(auxOR)){
  # Number of cases rho= 0 and indicator=1
  ntotal = length(which(data1$OR2 ==auxOR[i] & data1$indicator==1))
  
  # Number cases ARE.beta>1
  n.arebeta = length(which(data1$OR2 ==auxOR[i] & data1$indicator==1 & data1$ARE.beta>1)) 
  
  # % cases ARE.beta>1
  percent.arebeta = n.arebeta/ntotal
  
  dataguidelines.OR2.1[i,]<- c(auxOR[i], ntotal, n.arebeta, percent.arebeta)
}


# PATTERN %are>1 as a function of OR2, when OR1= 0.7
# Data obtained correspond to Figure 7
#####################################################################################

data2 <- subset(data, data$OR1 ==0.7)

dataguidelines.OR2.2 <- data.frame(0,0,0,0)
dataguidelines.OR2.2 <- rename(dataguidelines.OR2.2, c('c(0)'='Scenario',
                                                       X0 = 'OR2',
                                                       X0.1 = 'total cases', 
                                                       X0.2 = 'cases ARE.beta>1',
                                                       X0.3 = '% ARE.beta>1'
)
)

auxOR = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95)  
for(i in 1:length(auxOR)){
  # Number of cases rho= 0 and indicator=1
  ntotal = length(which(data2$OR2 ==auxOR[i] & data2$indicator==1))
  
  # Number cases ARE.beta>1
  n.arebeta = length(which(data2$OR2 ==auxOR[i] & data2$indicator==1 & data2$ARE.beta>1)) 
  
  # % cases ARE.beta>1
  percent.arebeta = n.arebeta/ntotal
  
  dataguidelines.OR2.2[i,]<- c(auxOR[i], ntotal, n.arebeta, percent.arebeta)
}


# PATTERN %are>1 as a function of OR2, when OR1= 0.8
# Data obtained correspond to Figure 8
#####################################################################################


data3 <- subset(data, data$OR1 ==0.8)

dataguidelines.OR2.3 <- data.frame(0,0,0,0)
dataguidelines.OR2.3 <- rename(dataguidelines.OR2.3, c('c(0)'='Scenario',
                                                       X0 = 'OR2',
                                                       X0.1 = 'total cases', 
                                                       X0.2 = 'cases ARE.beta>1',
                                                       X0.3 = '% ARE.beta>1'
)
)

auxOR = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95) 
for(i in 1:length(auxOR)){
  # Number of cases rho= 0 and indicator=1
  ntotal = length(which(data3$OR2 ==auxOR[i] & data3$indicator==1))
  
  # Number cases ARE.beta>1
  n.arebeta = length(which(data3$OR2 ==auxOR[i] & data3$indicator==1 & data3$ARE.beta>1)) 
  
  # % cases ARE.beta>1
  percent.arebeta = n.arebeta/ntotal
  
  dataguidelines.OR2.3[i,]<- c(auxOR[i], ntotal, n.arebeta, percent.arebeta)
}



#####################################################################################
#####################################################################################
# RECOMMENDATIONS FOR THE CHOICE OF THE PRIMARY ENDPOINT
# Section 6.3
#####################################################################################

# make sure the current working directory is the folder that contains DATABASE_Scenarios.csv
# load the data
db<-read.table(file="./DATABASE_Scenarios.csv",sep = ";", header = T, dec=",")

# the relevant endpoint (RE) is choosen when (are<=1) or composite endpoint (CE) (are > 1) 
# CE =1 indicates that the composite endpoint is choosen as the primary endpoint, otherwise CE =0
db$CE[db$ARE.beta>1] = 1 
db$CE[db$ARE.beta<=1] = 0 

# Define the categories of association, of strengths of the relative effect and of levels of frequency of the events:
# Correlation categories:
#   0= weak
#   1= medium-weak
#   2= medium-strong
#   3= strong
db$ClassRHO[db$rho<=0.2] = 0 
db$ClassRHO[db$rho>0.2 & db$rho<=0.5] = 1
db$ClassRHO[db$rho>0.5 & db$rho<=0.7] = 2
db$ClassRHO[db$rho>0.7] = 3

# OR categories:
#   0= weak
#   1= medium
#   2= strong
db$ClassOR2[db$OR2>=0.5 & db$OR2<0.7] = 0 
db$ClassOR2[db$OR2>=0.7 & db$OR2<0.9] = 1
db$ClassOR2[db$OR2>=0.9 & db$OR2<1] = 2

db$ClassOR1[db$OR1>=0.5 & db$OR1<0.7] = 0 
db$ClassOR1[db$OR1>=0.7 & db$OR1<0.9] = 1
db$ClassOR1[db$OR1>=0.9 & db$OR1<1] = 2


# event rate categories:
#   0= low
#   1= medium-low
#   2= medium-large
#   3= large
db$Classp1[db$p1.0<=0.025] = 0 
db$Classp1[db$p1.0>0.025 & db$p1.0<=0.05] = 1
db$Classp1[db$p1.0>0.05 & db$p1.0<=0.075] = 2
db$Classp1[db$p1.0>0.075] = 3

db$Classp2[db$p2.0<=0.025] = 0 
db$Classp2[db$p2.0>0.025 & db$p2.0<=0.05] = 1
db$Classp2[db$p2.0>0.05 & db$p2.0<=0.075] = 2
db$Classp2[db$p2.0>0.075] = 3


db$Classp1 <- as.factor(db$Classp1)
db$Classp2 <- as.factor(db$Classp2)
db$ClassOR1 <- as.factor(db$ClassOR1)
db$ClassOR2 <- as.factor(db$ClassOR2)
db$ClassRHO <- as.factor(db$ClassRHO)





# CASE I: when the correlation takes values between 0 and 1
#####################################################################################

dbExc = subset(db, db$rho> 0 & db$rho<1)


# Recommendations in terms of treatment effects of the relevant and the additional endpoint
# Data obtained correspond to Table 3 in the paper
#####################################################################################
TABLE3 = with(dbExc,tapply(CE,list(dbExc$ClassOR1,dbExc$ClassOR2),mean, na.rm=TRUE)) 
# write.csv2(TABLE3,'Table3dbExc.csv2')



# Recommendations in terms of degree of association between endpoints, 
# treatment effects of the relevant and the additional endpoint,
# event rates in control group for the relevant and additional endpoints
# Data obtained correspond to Table 4 in the paper
#####################################################################################
TABLE1 = with(dbExc,tapply(CE,list(dbExc$ClassOR2,dbExc$ClassRHO),mean, na.rm=TRUE)) 
# write.csv2(TABLE1,'Table1dbExc.csv2')

TABLE2 = with(dbExc,tapply(CE,list(dbExc$ClassOR1,dbExc$ClassRHO),mean, na.rm=TRUE)) 
# write.csv2(TABLE2,'Table2dbExc.csv2')

TABLE4 = with(dbExc,tapply(CE,list(dbExc$Classp1,dbExc$ClassRHO),mean, na.rm=TRUE))
# write.csv2(TABLE4,'Table4dbExc.csv2')

TABLE5 = with(dbExc,tapply(CE,list(dbExc$Classp2,dbExc$ClassRHO),mean, na.rm=TRUE))
# write.csv2(TABLE5,'Table5dbExc.csv2')



# CASE II: when the relevant and the additional endpoint are independent, i.e., rho=0
#####################################################################################
db0 = subset(db, db$rho==0)

# Recommendations in case of independence between the relevant and the additional endpoint
# (rho= 0) in terms of treatment effects of the relevant and the additional endpoint
# Data obtained correspond to Table 5 in the paper
#####################################################################################

with(db0,tapply(CE,list(db0$ClassOR1,db0$ClassOR2),mean, na.rm=TRUE))




#####################################################################################
#####################################################################################
# APPENDIX RECOMMENDATIONS using the threshold 1.1
#####################################################################################
# CE =1 if composite is choosen
db$CE_GL[db$ARE.beta>1.1] = 1 
db$CE_GL[db$ARE.beta<=1.1] = 0 


# CASE I: when the correlation takes values between 0 and 1
# using the threshold 1.1
#####################################################################################

# Recommendations in terms of treatment effects of the relevant and the additional endpoint
# the relevant endpoint (RE) is choosen when (ARE<1.1) or composite endpoint (CE) (ARE > 1.1) 
# Data obtained correspond to Table 6 in the paper
#####################################################################################

TABLE3 = with(dbExc,tapply(CE_GL,list(dbExc$ClassOR1,dbExc$ClassOR2),mean, na.rm=TRUE)) 
# write.csv2(TABLE3,'Table3dbExcGL.csv2')




# Recommendations in terms of degree of association between endpoints, 
# treatment effects of the relevant and the additional endpoint,
# event rates in control group for the relevant and additional endpoints
# Data obtained correspond to Table 7 in the paper
#####################################################################################

TABLE1 = with(dbExc,tapply(CE_GL,list(dbExc$ClassOR2,dbExc$ClassRHO),mean, na.rm=TRUE))
# write.csv2(TABLE1,'Table1dbExcGL.csv2')

TABLE2 = with(dbExc,tapply(CE_GL,list(dbExc$ClassOR1,dbExc$ClassRHO),mean, na.rm=TRUE)) 
# write.csv2(TABLE2,'Table2dbExcGL.csv2')

TABLE4 = with(dbExc,tapply(CE_GL,list(dbExc$Classp1,dbExc$ClassRHO),mean, na.rm=TRUE))
# write.csv2(TABLE4,'Table4dbExcGL.csv2')

TABLE5 = with(dbExc,tapply(CE_GL,list(dbExc$Classp2,dbExc$ClassRHO),mean, na.rm=TRUE))
# write.csv2(TABLE5,'Table5dbExcGL.csv2')




# CASE II: when the relevant and the additional endpoint are independent, i.e., rho=0
# using the threshold 1.1
#####################################################################################

db0 = subset(db, db$rho==0)

# Recommendations in case of independence between the relevant and the additional endpoint
# (rho= 0) in terms of treatment effects of the relevant and the additional endpoint
# using the threshold 1.1 (results not shown in the paper) 
#####################################################################################

with(db0,tapply(CE_GL,list(db0$ClassOR1,db0$ClassOR2),mean, na.rm=TRUE))





