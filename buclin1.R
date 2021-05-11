phe=read.csv(file="phe.csv",header=T,sep=";")
phe
str(phe)
dim(phe)
sum(is.na(phe))
summary(phe)

#############################FORMATTING
born=as.Date(phe$DNPatient)
born

library(lubridate)
today()
phe$age=(today()-born)/365
phe$age=round(phe$age,0)
phe$age
plot(table(phe$age))
hist(as.numeric(phe$age))
phe$arr=as.Date.factor(phe$DateArr)
str(phe$arr)
###A NEW DATA FRAME#########
phe1=data.frame(phe$IDPatient,born,phe$age,phe$arr,phe$Sexe,phe$MaladieGent,phe$Malignite,phe$Anamnese,phe$Variable,phe$Conc)
##Dealing with measurment concentration
hist(phe1$phe.Conc)
phe1$logconc=log(phe1$phe.Conc      )
hist(phe1$logconc)
describe(phe1$logconc)

boxplot(phe1$phe.Conc~as.factor(phe1$phe.Variable))
boxplot(phe1$logconc~as.factor(phe1$phe.Variable),cex.names=0.8,col=c(1,2,3,4,5,6),ylab="log Value")
abline(h=6,lwd=1,lty=3,col=3)
abline(h=7.5,lwd=1,lty=3,col=2)

legend("top",legend=c("UpLimitpheo","Error"),lty = 3,col = c(3,2),bty="n")
phe1
phe1$phe.Conc
##A wide format***
library(reshape2)
phe2=data.frame(phe1$phe.IDPatient,phe1$phe.Variable,phe1$phe.Conc)
str(phe2)
wide <- reshape(phe2, v.names = "phe1.phe.Conc", idvar = "phe1.phe.IDPatient",timevar = "phe1.phe.Variable", direction = "wide")
str(wide)
log(wide)
wide$`phe1.phe.Conc.Métanéphrine L`
7*516
library(psych)
pairs.panels(log(wide[,2:7]),main="all 1st measurment  of biomarker wFMT in log scale")
###LOG CONC VS AGE
length(phe1$age)
length(phe1$logconc)
plot(phe1$logconc[phe1$phe.Variable=="Normétanéphrine L"]~phe1$phe.age[phe1$phe.Variable=="Normétanéphrine L"])
line=lm(phe1$logconc[phe1$phe.Variable=="Normétanéphrine L"]~phe1$phe.age[phe1$phe.Variable=="Normétanéphrine L"])
abline(line,col=2)

list(phe1$phe.Variable)
met=subset(phe1,phe1$phe.Variable=="Normétanéphrine L")
library(lattice)
met$phe.age=as.factor(met$phe.age)
bwplot(met$logconc~met$phe.age)
bwplot(met$logconc~met$phe.age|met$phe.Sexe)
plot(met$logconc[met$phe.Anamnese=="Suivi phéochromocytome opéré"]~as.numeric(met$phe.age[met$phe.Anamnese=="Suivi phéochromocytome opéré"]))
reg1age=lm(met$logconc[met$phe.Anamnese=="Suivi phéochromocytome opéré"]~as.numeric(met$phe.age[met$phe.Anamnese=="Suivi phéochromocytome opéré"]))
abline(reg1age,col=2,lwd=3)
exp(0.00194)
##on doit cree une classe
summary(reg1age)

hist(met$phe.Conc)
summary(met$phe.Conc)
x=log(wide[,2:7]    )
summary(x)
xx=scale(x)
hist(xx,xlab = "SD of log scale N(0,1)",main="scaled version of all meta 1st measurements [.LOG]")
V=na.omit(xx)
sum(is.na(V))
KK=princomp(V)
summary(KK)
###2 biomarker explains 90%
biplot(KK,cex=0.5,main="PC1 & PC2 Biomarkers 1st measurments")
plot(KK,type="l")

library(lattice)
str(phe1$phe.Variable)
histogram( ~ logconc | phe.Variable, data = phe1,
           xlab = "Biomarker in log", type = "density",
           panel = function(x, ...) {
             panel.histogram(x, ...)
             panel.mathdensity(dmath = dnorm, col = "black",
                               args = list(mean=mean(x),sd=sd(x)))
           } )
histogram( ~ logconc | phe.Sexe+phe.Variable, data = phe1,
           xlab = "Biomarker in log and Sexe", type = "density",
           panel = function(x, ...) {
             panel.histogram(x, ...)
             panel.mathdensity(dmath = dnorm, col = "black",
                               args = list(mean=mean(x),sd=sd(x)))
           } )
str(phe1)
str(phe1$phe.Sexe)
#so the distribution shape of all Meta is independant of age
summary((phe1$phe.Sexe))
barplot(table(phe1$phe.Sexe))
##quzilibre male female et ?=12

#valeur de 604 dans Normétanéphrine L

##AMANESE
met$phe.Anamnese
barplot(table(met$phe.Anamnese),cex.axis=0.7,cex.names = 0.5)

######################Trajectory
summary(met)
install.packages("epicalc", repos = "http://medipe.psu.ac.th/epicalc")
library(epicalc)
install.packages("nmle")
library(nlme)
library(MASS)
library(epicalc)
library(lattice)
##fourrfodl plot
phe1
met
.data <- as.data.frame(met)
use(.data)
des()
.data
.data$phe.age=as.numeric(.data$phe.age)
followup.plot(id=as.factor(phe.IDPatient) , time=phe.arr, outcome=logconc, line.col="multicolor",lwd=1)

title(main="Patients trajectory on LOG [Normétanéphrine L] Biomarker") 
xyplot(logconc~phe.arr|phe.MaladieGent,groups=as.factor(phe.IDPatient),      type="l"       ,aspect=0.7,data=.data)
