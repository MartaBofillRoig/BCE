
################################################################
# SIMULATIONS SAMPLE SIZE FOR COMPOSITE BINARY ENDPOINTS
# Results H0 true
# oct-2018  
# Marta Bofill and Guadalupe Gómez
################################################################

setwd("C:/Users/marta.bofill/Simulations")
data<-read.table(file="./RESULTS_H0True_CompleteUnpooled.csv",sep = ";", header = T, dec=",")

data$diff.comp_true = Bahadur.composite(data$p1.1, data$p2.1, data$rho_true)-Bahadur.composite(data$p1.0, data$p2.0, data$rho_true)

data$effect_error = data$diff.comp_true-data$diffcomp

summary(data) 
data$indicator = as.factor(data$indicator)
data$v = as.factor(data$v) 

################################################################

db0 = subset(data, data$v==0)
db1 = subset(data, data$v==1)
db2 = subset(data, data$v==2)
db3 = subset(data, data$v==3)

case_plot = c("Correct correlation", "Weak correlation", "Moderate correlation", "Strong correlation")

windows(height = 14, width = 21)
par(mfrow = c(2, 3))

# RISK DIFFERENCE
plot(db0$error, db0$Test_Reject, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Alpha Risk Difference  (unpooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0, 0.035))  
points(db1$error, db1$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=2)  
points(db2$error, db2$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=3)  
points(db3$error, db3$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=4) 
abline(h=0.05, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=1, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


# RISK RATIO
plot(db0$error, db0$Test_Reject_RR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Alpha Risk Ratio  (unpooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0, 0.035))  
points(db1$error, db1$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=2)  
points(db2$error, db2$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=3)  
points(db3$error, db3$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=4) 
abline(h=0.05, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=1, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


# ODDS RATIO
plot(db0$error, db0$Test_Reject_OR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Alpha Odds Ratio  (unpooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0, 0.035))  
points(db1$error, db1$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=2)  
points(db2$error, db2$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=3)  
points(db3$error, db3$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=4) 
abline(h=0.05, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=1, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 

################################################################ 

setwd("C:/Users/marta.bofill/Simulations")
data<-read.table(file="./RESULTS_H0True_CompletePooled.csv",sep = ";", header = T, dec=",")


summary(data) 
data$indicator = as.factor(data$indicator)
data$v = as.factor(data$v) 

################################################################

db0 = subset(data, data$v==0)
db1 = subset(data, data$v==1)
db2 = subset(data, data$v==2)
db3 = subset(data, data$v==3)

summary(db0)
summary(db1)
summary(db2)
summary(db3)

################################################################

case_plot = c("Correct correlation", "Weak correlation", "Moderate correlation", "Strong correlation")
# windows(height = 7, width = 21)
# par(mfrow = c(1, 3))

# RISK DIFFERENCE
plot(db0$error, db0$Test_Reject, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Alpha Risk Difference  (pooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0, 0.035))  
points(db1$error, db1$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=2)  
points(db2$error, db2$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=3)  
points(db3$error, db3$Test_Reject, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=4) 
abline(h=0.05, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


# RISK RATIO
plot(db0$error, db0$Test_Reject_RR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Alpha Risk Ratio (pooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0, 0.035))  
points(db1$error, db1$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=2)  
points(db2$error, db2$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=3)  
points(db3$error, db3$Test_Reject_RR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=4) 
abline(h=0.05, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 


# ODDS RATIO
plot(db0$error, db0$Test_Reject_OR, cex = 1, pch=20, 
     xlab = expression(paste(rho)-paste(rho["true"])), 
     ylab="Empirical Alpha Odds Ratio (pooled variance)",
     xlim=c(-0.5,0.8), ylim=c(0, 0.035))  
points(db1$error, db1$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=2)  
points(db2$error, db2$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=3)  
points(db3$error, db3$Test_Reject_OR, cex = 1, pch=20, xlim=c(-0.5,0.8), ylim=c(0.04, 0.06), col=4) 
abline(h=0.05, lty=3)
abline(v=0, cex = 0.5,  lty=3)
legend("topright", legend = case_plot, col=c(1,2,3,4), ncol=2 , lwd=2, cex=0.8, pt.cex=1, inset=0.03, x.intersp=0.3, y.intersp=0.9) 

################################################################ 




