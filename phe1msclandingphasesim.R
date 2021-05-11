beta <- c(1.5,-0.5,0.2)
time <- seq(0,20,len=1000)

X <- cbind(1,log(time),log(time)^2)
length(X)
y <- as.numeric(X%*%beta)
length(y)
plot(time,y,type="l",lwd=2,col=2,ylab = "Log10 Biomarker",sub="Yi=Bo + function log(time) + log(time) ^2 ")
  text(17,6,labels="coef beta=1.5,-0.5,0.2",col=2)
