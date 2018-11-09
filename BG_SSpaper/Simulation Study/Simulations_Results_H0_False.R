
################################################################
# SIMULATIONS SAMPLE SIZE FOR COMPOSITE BINARY ENDPOINTS
# Results H0 false
# oct-2018  
# Marta Bofill and Guadalupe Gómez
################################################################

setwd("C:/Users/marta.bofill/Simulations")
data<-read.table(file="./RESULTS_H0False_CompletePooled.csv",sep = ";", header = T, dec=",")


summary(data) 
data$indicator = as.factor(data$indicator)
data$v = as.factor(data$v) 

################################################################

db0 = subset(data, data$v==0)
db1 = subset(data, data$v==1)
db2 = subset(data, data$v==2)
db3 = subset(data, data$v==3)



case_plot = c("Correct correlation", "Weak correlation", "Moderate correlation", "Strong correlation")

windows(height = 7, width = 21)
par(mfrow = c(1, 3))
plot(db0$error, db0$Test_Reject, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])),
     # xlab = expression(paste(rho["assumed"])-paste(rho["true"])),
     ylab="Empirical Power Risk Difference  (pooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0.65,1)
)  
points(db1$error, db1$Test_Reject, cex = 1, pch=20,
     # xlab="error", ylab="Test_Reject",
     xlim=c(-0.5,0.8), ylim=c(0.65,1), col=2
)  
points(db2$error, db2$Test_Reject, cex = 1, pch=20,
     # xlab="error", ylab="Test_Reject",
     xlim=c(-0.5,0.8), ylim=c(0.65,1), col=3
)  
points(db3$error, db3$Test_Reject, cex = 1, pch=20,
     # xlab="error", ylab="Test_Reject",
     xlim=c(-0.5,0.8), ylim=c(0.65,1), col=4
) 
abline(h=0.8, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


# RISK RATIO
plot(db0$error, db0$Test_Reject_RR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])),
     # xlab = expression(paste(rho["assumed"])-paste(rho["true"])),
     ylab="Empirical Power Relative Risk (pooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0.65,1)
)  
points(db1$error, db1$Test_Reject_RR, cex = 1, pch=20,
       # xlab="error", ylab="Test_Reject",
       xlim=c(-0.5,0.8), ylim=c(0.65,1), col=2
)  
points(db2$error, db2$Test_Reject_RR, cex = 1, pch=20,
       # xlab="error", ylab="Test_Reject",
       xlim=c(-0.5,0.8), ylim=c(0.65,1), col=3
)  
points(db3$error, db3$Test_Reject_RR, cex = 1, pch=20,
       # xlab="error", ylab="Test_Reject",
       xlim=c(-0.5,0.8), ylim=c(0.65,1), col=4
) 
abline(h=0.8, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 

# ODDS RATIO
plot(db0$error, db0$Test_Reject_OR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])),
     # xlab = expression(paste(rho["assumed"])-paste(rho["true"])),
     ylab="Empirical Power Odds Ratio  (pooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0.65,1)
)  
points(db1$error, db1$Test_Reject_OR, cex = 1, pch=20,
       # xlab="error", ylab="Test_Reject",
       xlim=c(-0.5,0.8), ylim=c(0.65,1), col=2
)  
points(db2$error, db2$Test_Reject_OR, cex = 1, pch=20,
       # xlab="error", ylab="Test_Reject",
       xlim=c(-0.5,0.8), ylim=c(0.65,1), col=3
)  
points(db3$error, db3$Test_Reject_OR, cex = 1, pch=20,
       # xlab="error", ylab="Test_Reject",
       xlim=c(-0.5,0.8), ylim=c(0.65,1), col=4
) 
abline(h=0.8, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


################################################################ 
################################################################
# MISSPECIFICATION ERROR: Correlation within the category AND Correlation out of the category

# within category
data$within=0
data$within[data$rho_assumed>data$rho_true & data$v==1]=1 
data$within[data$rho_assumed>data$rho_true & (data$rho_assumed-data$rho.rank/3)<data$rho_true & data$v==2]=2 
data$within[data$rho_assumed>data$rho_true & (data$rho_assumed-data$rho.rank/3)<data$rho_true & data$v==3]=3 
data$within[data$v==0]=0 


data$within = as.factor(data$within)
summary(data)


# Correlation within the category
################################################################
# under weak correlation assumption
summary(subset(data, data$within==1))
# under moderate correlation assumption
summary(subset(data, data$within==2))
# under strong correlation assumption
summary(subset(data, data$within==3))

# Correlation out the category
################################################################ 
# under weak correlation assumption
summary(subset(data, data$v==1 & data$within!=1))
# under moderate correlation assumption 
summary(subset(data, data$v==2 & data$within!=2))
# under strong correlation assumption 
summary(subset(data, data$v==3 & data$within!=3))


################################################################
################################################################

#############################################################
setwd("C:/Users/marta.bofill/Simulations")
data<-read.table(file="./RESULTS_H0False_CompleteUnpooled.csv",sep = ";", header = T, dec=",")


summary(data) 
data$indicator = as.factor(data$indicator)
data$v = as.factor(data$v) 

################################################################

db0 = subset(data, data$v==0)
db1 = subset(data, data$v==1)
db2 = subset(data, data$v==2)
db3 = subset(data, data$v==3)



case_plot = c("Correct correlation", "Weak correlation", "Moderate correlation", "Strong correlation")

# PRINT BOTH APPROACHES
# windows(height = 14, width = 21)
# par(mfrow = c(2, 3))

windows(height = 7, width = 21)
par(mfrow = c(1, 3))
plot(db0$error, db0$Test_Reject, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Power Risk Difference  (unpooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0.65,1)
)  
points(db1$error, db1$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=2)  
points(db2$error, db2$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=3)  
points(db3$error, db3$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=4) 
abline(h=0.8, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


# RISK RATIO
plot(db0$error, db0$Test_Reject_RR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), ylab="Empirical Power Relative Risk  (unpooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0.65,1)
)  
points(db1$error, db1$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=2)  
points(db2$error, db2$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=3)  
points(db3$error, db3$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=4) 
abline(h=0.8, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 

# ODDS RATIO
plot(db0$error, db0$Test_Reject_OR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Power Odds Ratio  (unpooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0.65,1)
)  
points(db1$error, db1$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=2)  
points(db2$error, db2$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=3)  
points(db3$error, db3$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.65,1), col=4) 
abline(h=0.8, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 

################################################################ 
################################################################
# MISSPECIFICATION ERROR: Correlation within the category AND Correlation out of the category

# within category
data$within=0
data$within[data$rho_assumed>data$rho_true & data$v==1]=1 
data$within[data$rho_assumed>data$rho_true & (data$rho_assumed-data$rho.rank/3)<data$rho_true & data$v==2]=2 
data$within[data$rho_assumed>data$rho_true & (data$rho_assumed-data$rho.rank/3)<data$rho_true & data$v==3]=3 
data$within[data$v==0]=0 


data$within = as.factor(data$within)
summary(data)


# Correlation within the category
################################################################
# under weak correlation assumption
summary(subset(data, data$within==1))
# under moderate correlation assumption
summary(subset(data, data$within==2))
# under strong correlation assumption
summary(subset(data, data$within==3))

# Correlation out the category
################################################################ 
# under weak correlation assumption
summary(subset(data, data$v==1 & data$within!=1))
# under moderate correlation assumption 
summary(subset(data, data$v==2 & data$within!=2))
# under strong correlation assumption 
summary(subset(data, data$v==3 & data$within!=3))


