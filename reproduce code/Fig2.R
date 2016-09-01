## make sure to set the directory correctly 
## i.e. directory contains the "bathtub.R" source code
setwd("C:/Users/wfu/Desktop/biostatistics/code")
source("bathtub.R")

set.seed(1)
bath1 <- Bathtub(n=100000, a = 0.01)
bath2 <- Bathtub(n=100000, a = 0.05)
bath3 <- Bathtub(n=100000, a = 0.1)
bath4 <- Bathtub(n=100000, a = 0.7)



##############################################################
####=====Code below gives the Desnity plot of survival time T=========
#################################################################
pdf("struct_trunc_plot2.pdf",width=10)

par(mfrow=c(2,3))

plot(density(rexp(100000,0.1)),main="Exponential",ylab="",xlab="T",ylim=c(0,0.75),xlim=c(0,20),lwd=2)
points(density(rexp(100000,0.23)),type="l",col="red",lty=2,lwd=2)
points(density(rexp(100000,0.4)),type="l",col="blue",lty=6,lwd=2)
points(density(rexp(100000,0.9)),type="l",col="purple",lty=5,lwd=2)
legend("topright", c(expression(tilde(T)[1]),expression(tilde(T)[2]),expression(tilde(T)[3]),expression(tilde(T)[4])), col = c("black","red","blue","purple"),text.col = c("black","red","blue","purple"), lty = c(1, 2, 6,5))



plot(density(rweibull(100000, 0.9, 7)),main="Weibull decreasing hazard",ylab="",xlab="T",ylim=c(0,0.9),xlim=c(0,20),lwd=2)
points(density(rweibull(100000, 0.9, 3)),type="l",col="red",lty=2,lwd=2)
points(density(rweibull(100000, 0.9, 2.5)),type="l",col="blue",lty=6,lwd=2)
points(density(rweibull(100000, 0.9, 1)),type="l",col="purple",lty=5,lwd=2)
legend("topright", c(expression(tilde(T)[1]),expression(tilde(T)[2]),expression(tilde(T)[3]),expression(tilde(T)[4])), col = c("black","red","blue","purple"),text.col = c("black","red","blue","purple"), lty = c(1, 2, 6,5))


plot(density(rweibull(100000, 3, 10)),main="Weibull increasing hazard",ylab="",xlab="T",ylim=c(0,0.6),xlim=c(0,20),lwd=2)
points(density(rweibull(100000, 3, 6.2)),type="l",col="red",lty=2,lwd=2)
points(density(rweibull(100000, 3, 4.3)),type="l",col="blue",lty=6,lwd=2)
points(density(rweibull(100000, 3, 2)),type="l",col="purple",lty=5,lwd=2)
legend("topright", c(expression(tilde(T)[1]),expression(tilde(T)[2]),expression(tilde(T)[3]),expression(tilde(T)[4])), col = c("black","red","blue","purple"),text.col = c("black","red","blue","purple"), lty = c(1, 2, 6,5))


plot(density(rlnorm(100000, meanlog = 2.0, sdlog = 0.3)),main="Lognormal",ylab="",xlab="T",ylim=c(0,0.6),xlim=c(0,20),lwd=2)
points(density(rlnorm(100000, 1.7, 0.2)),type="l",col="red",lty=2,lwd=2)
points(density(rlnorm(100000, 1.3, 0.3)),type="l",col="blue",lty=6,lwd=2)
points(density(rlnorm(100000, 0.5, 0.5)),type="l",col="purple",lty=5,lwd=2)
legend("topright", c(expression(tilde(T)[1]),expression(tilde(T)[2]),expression(tilde(T)[3]),expression(tilde(T)[4])), col = c("black","red","blue","purple"),text.col = c("black","red","blue","purple"), lty = c(1, 2, 6,5))

plot(density(bath1),main="Bathtub",ylab="",xlab="T",ylim=c(0,0.65),xlim=c(0,20),lwd=2)
points(density(bath2),type="l",col="red",lty=2,lwd=2)
points(density(bath3),type="l",col="blue",lty=6,lwd=2)
points(density(bath4),type="l",col="purple",lty=5,lwd=2)
legend("topright", c(expression(tilde(T)[1]),expression(tilde(T)[2]),expression(tilde(T)[3]),expression(tilde(T)[4])), col = c("black","red","blue","purple"),text.col = c("black","red","blue","purple"), lty = c(1, 2, 6,5))

dev.off()

###===============================================================