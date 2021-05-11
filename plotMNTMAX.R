read.csv("phemdy1",header=TRUE)
pheMod=read.csv(file="phemdy1.csv",header=T)
pheMod
summary(pheMod)
par(mfrow=c(1,2))
plot(pheMod$MaxMTT,col=1,ylim=c(0,40),main="Max Value of Biomarkers",xlab="measures")
legend("center",legend=c("MTT","MNT","NMNT"),col=c(1,2,3),lwd=1,cex=0.5,bty="n")
points(phe$MaxMNT,col=2,ylim=c(0,40))
points(phe$MaxNMNT,col=3,ylim=c(0,40))
plot(log(pheMod$MaxMTT),col=1,ylim=c(0,4),main="Max Value of LOG[e] Biomarkers",xlab="measures")
legend("center",legend=c("MTT","MNT","NMNT"),col=c(1,2,3),lwd=1,cex=0.5,bty="n")
points(log(phe$MaxMNT),col=2,ylim=c(0,4))
points(log(phe$MaxNMNT),col=3,ylim=c(0,4))

