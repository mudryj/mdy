library(lattice)
getwd()
setwd("~/buclin")
md1 <- read.table("md1.csv",header=TRUE,sep=",")         # charge donnÃ©es
str(md1)
load("md1.Rdata")
names(md1)

md1$DateArr <- as.Date(md1$DateArr,format="%Y-%m-%d")    # converti DateArr en date
md1$Op_prec <- as.Date(md1$Op_prec,format="%Y-%m-%d")    # converti Op_prec en date
md1$Time <- as.numeric(md1$DateArr - md1$Op_prec)/365.25 # Nb annÃ©es depuis opÃ©ration
md1$Time

XI=md1$Time>0  ##°enlevé les obs 0-3moilength(XI)
XI
XI*1
##selcting positive time ID
md11=subset(md1,XI==1)
md11$Time#check only pos time occur
length(md11$Time)
XI025=md11$Time<0.25  ##°enlevé les obs 0-3moilength(XI)
XI025*1##indicatrice opf data which are timne below 3 months
sum(XI025*1)
md11$XI025=XI025*1
md11
md11=subset(md11,XI025==0)##removing XI with valu 1 = interval 0-0.25Y
dim(md11)
md11$Time
md11$Time##only patient with positive timne greater than 0.25 Y

md11$gen <- as.factor(md11$gen)                           # convertie gen en facteur
str(md11$gen)#ok 5 levels factor
# Supprimer les NA pour MNL
md2 <- subset(md11, !is.na(MNL))
sum(is.na(md2$MNL)) ##Plus de NA ds MNL
dim(md2)
dotplot(md2$IDPatient~md2$Time,main="502 longitudinal measures | ID 1ops only positive TIME",col=rainbow(91))#nobrs of longitidunal data per ID
# Intercepts et pentes individuelles (sur donnÃ©es post-op uniquement!)
md2$Time
save(md2,file="md2.Rdata")
id.postop2 <- as.numeric(names(which(tapply(md2$Time,md2$IDPatient,function(t){sum(t>0)>1})))) # ID patient ayant au moins 2 data post-op
id.postop2#patients with at least 2 data
N.postop2 <- length(id.postop2)
N.postop2
intercepts <- slopes <- numeric(N.postop2)

for(i in 1:N.postop2){
  lmod <- lm(log10(MNL)~Time,data=subset(md2,IDPatient==id.postop2[i] & Time>0)) # lm que sur les donnÃ©es post-op: Time>0
  intercepts[i] <- coef(lmod)[1]
  slopes[i] <- coef(lmod)[2]
}


slopes
cor(intercepts,slopes)


###is correlation of an lm centerd vs uncentered different from eachoth?
intercepts <- slopes <- numeric(N.postop2)
mean(md2$Time)
md2$Time
for(i in 1:N.postop2){
  lmod <- lm(log10(MNL)~I(Time-2),data=subset(md2,IDPatient==id.postop2[i] & Time>0)) # lm que sur les donnÃ©es post-op: Time>0
  intercepts[i] <- coef(lmod)[1]
  slopes[i] <- coef(lmod)[2]
}
slopes
intercepts
cor(slopes,intercepts)
dev.off()
par(pty="s")
plot(c(0,15),c(-5,3),type="n",xlab="Time after 1st surgery (positive Time)",ylab="log10 MNL",asp = 1.3)
for(i in 1:N.postop2){
  abline(intercepts[i],slopes[i],col=rainbow(N.postop2)[i])
}
cor(intercepts,slopes)

# Profils individuels after Time>0.25
id <- unique(md2$IDPatient)

id
N <- length(id)
N
plot(range(md2$Time),range(log10(md2$MNL),na.rm=T),type="n",xlab="Nb annÃ©es depuis opÃ©ration",ylab="log10(MNL)",xlim=c(-0.5,20))
for(i in 1:N){
  tmp <- subset(md2,IDPatient==id[i])
  lines(tmp$Time,log10(tmp$MNL),col=rainbow(N)[i],type="o")
}
abline(v=0,lty=2)
abline(h=log10(0.83),lty=2,col="red")
mean(log10(md2$MNL))

#some basic analysis
par(mfrow=c(1,2))
dotplot(with(md2,tapply(log(MNL),gen,mean)),sub="89 ID",main="Mean of log10MNL |gen")

dotplot(with(md2,tapply(log10(MNL),IDPatient,mean)),main="Mean of log10MNL | IDPatient",panel=function(x,y,...){
  panel.dotplot(x,y,...)
  meanval=mean(log10(md2$MNL))
  panel.abline(v=meanval,col.line = 2)
  panel.abline(v=-1.53,col.line = 2,lty=3)
  panel.text(x=290,y=250,labels="grandmean",col="red",cex=0.6)
  panel.text(x=150,y=250,labels="lower lim",col="red",cex=0.6)
  })
dotplot(with(md2,tapply(log(MNL),gen,var)),sub="89 ID ",main="Sample Variance of log10MNL | gen")
K=with(md2,tapply(log10(MNL),IDPatient,var))
K
dotplot(with(md2,tapply(log10(MNL),factor(IDPatient),var)),cex.axis=0.4,sub="Outliers ID 1499 & 9644 are MEN Genes",main="Sample Variance of log10MNL | IDPatient",panel=function(x,y,...){
  panel.dotplot(x,y,...)
  meanval=var(log10(md2$MNL))
   })

var(log10(md2$MNL))
par(mfrow=c(1,2))
plot.design(log10(md1$MNL)~md1$gen,main="mean of log10 MNL all time",cex.main=0.8)
plot.design(log10(md2$MNL)~md2$gen,main="mean log10 MNL Time >0.25 Year",cex.main=0.8)
library(lattice)
# Time format
length(md2$Time)
md2$Time0 <- (md2$Time>0)*1   
(md2$Time0)
# indicatrice qui vaut 1 si Time>0 et 0 sinon
md2$Time0
md2$Timep <- md2$Time>=0.25  # temps positif: Ã©gal Ã  Time si Time>0, sinon Ã©gal Ã  0 Temps après ops
md2$Timep-md2$Time##same modeified on 27.11 call aziz in md11=ms2
md2$Timep=md2$Timep*md2$Time
md2$Timep
save(md2,file="md2.Rdata")

plot(as.numeric(md2$IDPatient),md2$Timep,main="suivi patient post ops",xlab="ID PATIENT",ylab="Temps années",pch=19,col=4)
abline(h=10,col=2,lty=2)
legend(20000,20,legend=("ops après 10 ans"),col=2,lwd=2,lty=2,cex=0.7,bty="n")
# LMM
##slope pop by genre lm
library(nlme)
library(lmerTest)

##slope and int families for post ops data patient
md3=subset(md2,md2$Timep>0)####temps positive > 0.35 3 mois
md3=subset(md3,md3$Timep>0.25)

plot(intervals(lmID))
plot(confint(lmID,na.rm=T))
confint(lmID$`1499`)##why so large ?
md2_1499=md2$IDPatient!="1499"
md2_1499
md21499=md2[md2_1499,]
mean(md21499$Timep)
unique(md21499$IDPatient)##patient
lmIDwitjout1499=lmList(log(MNL) ~ Timep | IDPatient ,na.action=na.fail,data=md21499)
plot(confint(lmIDwitjout1499),scale=list(y=list(y.cex=0.4)),ylab="ID ",pch=22,xlim=c(-6,1))
lmIDwitjout1499centered=lmList(log(MNL) ~ I(Timep-4.29) | IDPatient ,na.action=na.fail,data=md21499)
plot(confint(lmIDwitjout1499centered),scale=list(y=list(y.cex=0.4)),ylab="ID ",pch=22,xlim=c(-6,1))

res=data.frame(Subejct=rownames(coef(lmIDwitjout1499)),coef(lmIDwitjout1499))
res=na.omit(res)
cor(res$X.Intercept.,res$Timep)
?cor
res1=data.frame(Subejct=rownames(coef(lmIDwitjout1499centered)),coef(lmIDwitjout1499centered))
res1=na.omit(res1)
res1
cor(res1$X.Intercept.,res1$I.Timep...4.29.)
ID=unique(md1$IDPatient)
IDs=sample(ID,20)
IDs
library(lattice)
MaxMNL
maxlog10mnl=log10(MaxMNL)
xyplot(log10(MNL)~Timep,groups=IDPatient,data=subset(md3,md3$IDPatient%in%IDs),ylim=c(-1.4,0),type=c("o"),lwd=1.5,pch=19,aspect=0.5,auto.key = list(space="right",cex=0.5,title="Patient ID"),xlab="Time in year post surgery",
  panel=function(x,y,...){
    panel.xyplot(x,y,...)
  panel.abline(h=maxlog10mnl,col.line = 2,lty=3)
  panel.text(0.5,-0.1,"URL",col=2)
})
##MNT biomarker
MaxMNT
maxlog10mnt=log10(MaxMNT)
md2$
md3
xyplot(log10(MNT2014)~Timep,groups=IDPatient,data=subset(md2,md2$IDPatient%in%IDs),type=c("o"),lwd=1.5,pch=19,aspect=0.5,auto.key = list(space="right",cex=0.5,title="Patient ID"),xlab="Time in year post surgery",
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h=maxlog10mnt,col.line = 2,lty=3)
         panel.text(0.5,-0.1,"URL",col=2)
       })
##NMNL
MaxNMNL
maxlog10nmnl=log10(MaxNMNL)       

xyplot(log10(NMNL)~Timep,groups=IDPatient,data=subset(md3,md3$IDPatient%in%IDs),ylim=c(-1.4,1),type=c("o"),lwd=1.5,pch=19,aspect=0.5,auto.key = list(space="right",cex=0.5,title="Patient ID"),xlab="Time in year post surgery",
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h=maxlog10nmnl,col.line = 2,lty=3)
         panel.text(0.5,-0.1,"URL",col=2)
       })


MaxNMNT
maxlog10nmnt=log10(MaxNMNT)    
xyplot(log10(NMNT2014)~Timep,groups=IDPatient,data=subset(md2,md3$IDPatient%in%IDs),type=c("o"),lwd=1.5,pch=19,aspect=0.5,auto.key = list(space="right",cex=0.5,title="Patient ID"),xlab="Time in year post surgery",
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h=maxlog10nmnt,col.line = 2,lty=3)
         panel.text(0.5,-0.1,"URL",col=2)
       })
##MTL
MaxMTL
maxlog10mtl=log10(MaxMTL)           
xyplot(log10(MTL)~Timep,groups=IDPatient,data=subset(md3,md3$IDPatient%in%IDs),ylim=c(-2.8,0),type=c("o"),lwd=1.5,pch=19,aspect=0.5,auto.key = list(space="right",cex=0.5,title="Patient ID"),xlab="Time in year post surgery",
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h=maxlog10mtl,col.line = 2,lty=3)
         panel.text(0.5,-1,"URL",col=2)
       })
MaxMTT
maxlog10mtt=log10(MaxMTT) 
xyplot(log10(MTT2014)~Timep,groups=IDPatient,data=subset(md3,md3$IDPatient%in%IDs),type=c("o"),lwd=1.5,pch=19,aspect=0.5,auto.key = list(space="right",cex=0.5,title="Patient ID"),xlab="Time in year post surgery",
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h=maxlog10mtt,col.line = 2,lty=3)
         panel.text(0.5,-1,"URL",col=2)
       })
save(md2,file="md2.Rdata")
save(md3,file="md3.Rdata")
##threshold cens check
sum(md2$MNL<0.03)
dim(md2)
md2$MNL[which(md2$MNL<=0.03)]=0.03
my=xyplot(log10(MNL)~Timep,groups=IDPatient,ylim=c(-2,0),xlab="time after ops Y",ylab="LOG10 MNL",main="MNL trajectories of 89 ID patients",data=md2,lwd=2,type="o",col=1:10,auto.key=list(size=0.7))#reg line per ID a iwthin subkevt reg
my##censoting visible at 0.03








########EFFECTIF############

cluster <- as.numeric(as.factor(md2$IDPatient))
cluster
age <- tapply(md2$age,cluster,unique)
age
tapply(md2$gen,cluster,unique)
gen <- levels(md2$gen)[tapply(md2$gen,cluster,unique)]
gen
table(gen)

library(lme4)
library(lmerTest)
library(nlme)
###LMIXED MODELS
str(md2)
##checking covariates: GENE AGE
unique(md2$IDPatient)##89 PD
X=groupedData(log10(MNL)~cut(age,5)|gen,data=as.data.frame(md2),order.groups = T)
plot(X)
Y=groupedData(log10(MNL)~gen|IDPatient,data=as.data.frame(md2),order.groups = T)
plot(Y,outer=~cut(age,5),main="89ID valeur Log10MNL post ops 1 seul ops Temps Pos",ylab=list(label="IDPatient",cex=0.7),scales=(y=list(cex=0.4)))
plot(X,outer=~gen,main="valeur Log10MNL post ops 1 seul ops 89IDP",ylab=list(label="IDPatient",cex=0.7),scales=(y=list(cex=0.4)))
Z=groupedData(log10(NMNL)~cut(age,5)|IDPatient,data=as.data.frame(md2),order.groups = T)
ZZ=groupedData(log10(NMNL)~gen|IDPatient,data=as.data.frame(md2),order.groups = T)
plot(Z,outer=~cut(age,5),main="valeur NMNL post ops 1 seul ops 91IDP",ylab=list(label="IDPatient",cex=0.7),scales=(y=list(cex=0.4)))
plot(ZZ,outer=~cut(age,5),main="valeur NMNL post ops 1 seul ops 91IDP",ylab=list(label="IDPatient",cex=0.7),scales=(y=list(cex=0.4)))

#######MIXED MODELS 
lmeC=lmeControl(msMaxIter=100,opt="optim")
#just mean and random intercept
mix00RI=lme(log10(MNL)~1,data=md2,random=~1|IDPatient,control=lmeC)###random intercept only on postive time
summary(mix00)
mix00RIAS=lme(log10(MNL)~1,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")###random intercept only on postive time
summary(mix00RIAS)
0.26^2/(0.26^2+0.12^2)##ICC
anova.lme(mix00RI,mix00RIAS)
Mean=with(md2,sapply(split(log10(MNL),IDPatient),mean))
dotplot(Mean)##mean by ID
summary(mix00RIASnocor)
anova(mix00RIAS,mix00RIASnocor)
##only random intercept on Time pos Timep 
##exact LRT
install.packages("RLRSim")
require(RLRSim)
exactRLRT(finalmodremlmer)


xyplot(log10(MNL)~Timep|md2$gen,groups=IDPatient,data=md2,type="b",pch=20)
#selon tel avec Aziz 30.11 faire un MIXX sans predicecteur et voir les raneff vs covar.

##centrons le raneff aussi
mixx2=lme(log10(MNL)~Timcenpos,data=md2,random=~Timcenpos|IDPatient,control=lmeC)###random intercept only on postive time
summary(mixx2)##ne change preseque rien

##square tersm
Timepc=md2$Timep-0.25
Timepc
mixRIASsq=lme(log10(MNL)~(Timep)+I(Timep^2),data=md2,random=~I(Timep^2)+Timep|IDPatient,control=lmeC,method="ML")###random intercept only on postive time
summary(mixRIASsq)
anova(mix00RIAS,mixRIASsq)
                 ##Grouped data
##random RIS all cov
mixxgen=lme(log10(MNL)~Timep*gen,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")#)
summary(mixxgen)#not valid intercation gen
mixxgenage=lme(log10(MNL)~Timep*gen*age,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")#)
summary(mixxgenage)
mixxgenagerandomgen=lme(log10(MNL)~Timep*gen*age,data=md2,random=~Timep|IDPatient/gen,control=lmeC,method="ML")#)
summary(mixxgenagerandomgen)
mixrandomgencr=lmer(log10(MNL)~Timep*gen*age+(Timep|IDPatient)+(1|gen),data=md2)
summary(mixxgenagerandomgen)
##lmer
mixxar=update(mixx,cor=corAR1(),data=md2)
mix00RIAS
summary(mix00RIAS)
stargazer(mix00RIAS,mixRIASsq,mixxgenage,mixxgenagerandomgen,mixrandomgencr,type="html",out="Stargazermodsum.doc",digits=2)
anova(mixRIASsq,mix00RIAS,mixxgenage,mixxgenagerandomgen)
########################33
########################33MODEL NMNL##############33

########################################3MODEL MNL
################################################
mixgenagenointer3and2=lme(log10(MNL)~Timep+gen+age+gen*age,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")#)
mixgenagenointer3and2
anova(mixgenagenointer3and2,type="marginal")
finalmod=lme(log10(MNL)~Timep+gen+age+gen:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
finalmodreml=lme(log10(MNL)~Timep+gen+age+gen:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="REML")
finalmodremlmer=lmer(log10(MNL)~Timep+gen+age+gen:age+(Timep|IDPatient),data=md2)

summary(finalmod)
summary(finalmodreml)
fixef(finalmod)
fixef(finalmodreml)##REML ML no diff
round(anova(finalmodreml,type="marginal",test="LRT"),2)
round(anova(finalmodremlmer,ddf="Kenward-Rogger"),2)
summary(finalmodreml)

#other contsrats
library(nlme)
contrasts(md2$gennew)
finalmodreml2=lme(log10(MNL)~Timep+gennew+age+gennew:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="REML")
summary(finalmodreml2)
stargazer(finalmodreml2,type="html",out="finalmodsum.doc",digits=2,cor=T)
##relvel for sporadic
(md2$gen)
md2$gennew=relevel(md2$gen,ref="sporadi")
contrasts(md2$gennew)
finalmodrelevelspor=lme(log10(MNL)~Timep+gennew+age+gennew:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="REML")
summary(finalmodrelevelspor)
anova(finalmodrelevelspor,type="marginal")

anova(mixgenagenointer3and2,finalmod)

finalmodspo=lme(log10(NMNL)~Timep+gen+age,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
summary(finalmodspo)
anova(finalmodspo,type="marginal")
anova(finalmod,finalmodspo)###are the same iwth cotrast change
md2$Timep
##sum contrast for typeIII Anova

finalmodsumcont=lme(log10(MNL)~Timep+gennew2+age+age:gennew2,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
anova(finalmodsumcont,type="marginal")
summary(finalmodsumcont)##no gene sign?
finalmod
##in such model with sign intercept should not be removed
#RMLE fit

###fixef plot
library(car)
plot(allEffects(finalmod),cex.text=0.4)
plot(allEffects(groups="age",finalmod),cex.text=0.4)
str(md2$age)
library(effects)
plot(effect("age",finalmod))


##post HOC test after anova(finalmod,type = "marginal")
install.packages("lme4")
library(lme4)
library(lmerTest)
lsmeans(finalmodreml,adjust="tukey")
install.packages("multicomp")
library(multicomp)
head(fitted(finalmod,level=0:1))
head(fitted(finalmod,level=0))
fitval=fitted(finalmod,level=0)
plot(log10(md3$MNL),fitval,ylim=c(-2,0),xlim=c(-2,0))
abline(0,1)
plot(log10(md3$MNL))

length(fitval)
dim(md3)
PT=pairwise.t.test(fitval,md3$gen,method=bonferroni)
PT
library(stargazer)
print(PT)
summary(finalmod)
###############################resid check
par(mfrow=c(1,1))
C=boxplot(resid(finalmod,type="p")~md2$IDPatient,main="Pearson type residuals plot ",col=8,pch=19,cex.main=0.7,xlab="ID")
abline(h=0,col=2,lty=2)
K=boxplot(resid(finalmod,type="p")~md2$gen,main="Pearson type residuals plot ",col=8,pch=19,cex.main=0.7,xlab="ID")
abline(h=0,col=2,lty=2)
V=boxplot(resid(finalmod,type="p")~cut(md2$age,5),main="Pearson type residuals plot ",col=8,pch=19,cex.main=0.7,xlab="ID")
C$out
K$out
library(lattice)
xyplot(resid(finalmod,type="p")~md2$Timep|md2$Sexi,type=c("p","smooth"),ylab="Pearson residuals",xlab="Time in Years",pch=20,panel=function(...){
  panel.abline(h=2,lty=2)
  panel.xyplot(...)
  })
xyplot(resid(finalmod,type="p")~md2$Timep|cut(md2$age,5),type=c("p","smooth"),ylab="Pearson residuals",xlab="Time in Years | Age classes",pch=20,panel=function(...){
  panel.abline(h=2,lty=2)
  panel.xyplot(...)
})
xyplot(resid(finalmod)~fitted(finalmod)|md2$gen,type=c("p","r"),pch=20)
qqnorm(resid(finalmod,type="p"),pch=19)
qqline(resid(finalmod),pch=19)
qqnorm(finalmodreml,~resid(.)|gen,pch=20)
hist(resid(finalmod),col=8,xlab=" raw residuals")
plot(Variogram(finalmod,pch=19,form=~Timep|IDPatient),xlab="Time",col.line=c(2),pch=20)
plot(Variogram(finalmod,pch=19,form=~Timep|IDPatient),xlab="Time",col.line=c(3),resType="p",pch=20,robust=T,resType="n")
finalmodcorrel=update(finalmod,corr=corExp(form=~Timep,nugget=T))
finalmodcorrel2=update(finalmod,correlation=corAR1())

summary(finalmodcorrel)
summary(finalmodcorrel2)
plot(Variogram(finalmodcorrel,form=~Timep|IDPatient),pch=19,xlab="Time",col.line=c(3),resType="n")
plot(Variogram(finalmodcorrel2,form=~Timep|IDPatient),pch=19,xlab="Time",col.line=c(3),resType="n",abline=(1))

jack=resid(finalmod,type="pearson")
jack
jack[which.max(abs(jack))]
range(jack)
library(outliers)
grubbs.test(resid(finalmod))
##removing 1499
md21499=md2[md2_1499,]
md21499$IDPatient##1499 removed
finalmod1499=lme(log10(MNL)~Timep+gen+age+gen:age,data=md21499,random=~Timep|IDPatient,control=lmeC,method="ML")
summary(finalmod1499)
library(nlmeU)
library(stargazer)
K=cbind("model"=fixef(finalmod),"without 1499"=fixef(finalmod1499))
stargazer(K,type="html",out="m499.doc")
getwd()
####RANEF ##################

plot(ranef(mix0),main="Log10(MNL)~Timep,data=md2,random= ~1|IDPatient 91ID ops 1 fois",cex.main=0.2)
qqnorm(mix0,~ranef(.))##ri
qqnorm(mixx,~ranef(.))#RIS
dev.off()
plot(finalmod,id=0.05)
get=getME(finalmodremlmer,"Z")
head(get)


##sexe
sexops=read.csv("sexops.csv",header=T,sep=",")
dim(sexops)
sexops2=md2

MEME=merge(sexops,sexops2)
dim(MEME)
order(MEME$Sexi)
MEME1=MEME[order(MEME$IDPatient),]

    md2$Sexi= MEME1$Sexi
md2$Sexi

finalmodsexreml=lme(log10(MNL)~Timep+Sexi+gen+age+gen:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
summary(finalmodsexreml)
IDD=sample(unique(md2$IDPatient),5)
IDD
finalmodsex5ID=lme(log10(MNL)~Timep+gen+age+gen:age,data=subset(md2,IDPatient%in%IDD),random=~Timep|IDPatient,control=lmeC,method="ML")
fit=predict(finalmodreml)
xyplot(fit~md2$Timep,groups=md2$IDPatient,data=md2,type=c("b"),pch=19,col=rainbow(89),xlim=c(10,35),grid=TRUE,ylab="predicted value",xlab="Time 10-35 Y",panel = function(x,y,...,fit){
  panel.xyplot(x,y,...)
    panel.points(log10(MNL)~md2$Timep,groups=md2$IDPatient,type="p",col=4,pch="o")
})
y=predict(finalmodreml,level=1)
log10(MNL)~Timep
points(md2$Timep,fitted(finalmod,level=1),col=rainbow(89),type=c("p"),pch="o")
anova(finalmodsexreml,type="marginal")
drop1(finalmodsexreml,test="Chisq")###type II sum ssq
plot(resid(finalmodsexreml),fitted(finalmodsexreml))##pattern still exist
library(lattice)
xyplot(fitted(finalmodsexreml)~resid(finalmodsexreml)|md2$Sexi)
###################################RANEF
parmfrow=c(1,1)
plot(ranef(finalmodremlmer)$IDPatient,pch=20)
qqnorm(finalmodreml,~ranef(.),pch=20)##ri
qqmath(~ranef(finalmodremlmer)|gen,md2)
names(ranef(finalmod))
dev.off()
V=ranef(finalmodreml)[1]
V
par(mfrow=c(1,2))
boxplot(V$`(Intercept)`~gen,col=8,pch=19,ylab="Ranef Intercept")
abline(h=0,col=2)
boxplot(ranef(finalmodreml)$Timep~gen,col=8,pch=19,ylab="Ranef Slope")
abline(h=0,col=2)

W=ranef(finalmodreml)
par(mfrow=c(1,2))
boxplot(V$`(Intercept)`~cut(age,5),col=8,pch=19,xlab="Ages class",ylab="Ranef Intercept",cex.axis=0.7)
abline(h=0,col=2)
boxplot(W$Timep~cut(age,5),col=8,pch=19,ylab="Ranef Slope",xlab="Ages class",cex.axis=0.7)

abline(h=0,col=2)
dev.off()

##ranef correlatef
cor(ranef(finalmod))
###############################
############################3
#####################
#######################
install.packages("merTools")
library(merTools)
B=bootstrap(finalmodremlmer,B=100)
confint(finalmodremlmer,method="boot")
summary(finalmodremlmer)
formula(finalmodremlmer)
X=model.matrix(finalmodremlmer)
bet=fixef(finalmodremlmer)
X%*%(bet)
plot(fitted(finalmodreml)~md2$Timep,col=as.integer(md2$IDPatient),pch=20)
sim1=sample(fitted(finalmodreml),1000,replace=TRUE)

md=sample(md2,1000,replace=TRUE)
dim(md)
md
X=simulate(finalmod,nsim=100)
str(X)






fit=fitted(finalmodreml)
res=resid(finalmodreml)
xyplot(res~fit|cut(age,5),ylab="residuals",xlab="fitted value",abline=c(0,col=2),pch=20)


#################################33
##bootstapping RANEF
mod1=lme(log10(MNL)~Timep+gen+age+age:gen,data=md2,random=~1|IDPatient,control=lmeC,method="ML")
mod2=lme(log10(MNL)~Timep+gen+age+age:gen,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
anova(mod1,mod2)
set.seed(1)
n.r=300
data.r=md2
my.y.nr=simulate(mod1,nsim=n.r)
my.y.nr
lrt.r=NULL
for (rw in 1:n.r){
  data.r$dist=my.y.nr$REML[,rw]
    mod1=lme(log10(MNL)~Timep+gen+age+age:gen,data=md2,random=~1|IDPatient,control=lmeC,method="ML")
  mod2=lme(log10(MNL)~Timep+gen+age+age:gen,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
  lrt.r[rw]=as.numeric(2*(logLik(mod2)-logLik(mod1)))
}
summary(finalmod)


#Random intercept assume pdSymm()
finalmodreml$varFix
fixef(finalmodreml)



##AUGPRED
md2$MNL10lg=log10(md2$MNL)
md2$MNL10lg
groupedData(MNL10lg~Timep|IDPatient,data=md2)
finalaug=lme(MNL10lg~Timep+gen+age+gen:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="REML")
summary(finalaug)
newdata=data.frame(Timep=md2$Timep,gen=md2$gen,age=md2$age)
newdata
predict(finalmodreml,level=1,newdata = newdata)
K=predict(finalmodreml,level=0:1)
K
K=data.frame(K)





mixxint=intervals(finalmod)##�fiexeffe intervall
install.packages("stargazer")
library(stargazer)
stargazer(mixxint,type="html",out="confint.doc",digits=2)
mixxint
hist(finalmod$residuals,type="pearson")
#################EFFECTPLOT CAR J FOX

library(effects)
par(mfrow=c(2,4))
plot(allEffects(finalmod))


##efectplot for linear age
finalmodmod=lme(log10(MNL)~Timep+gen+age+age:gen,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
summary(finalmodmod)
plot(effect("gen",finalmodmod))
?effect
finalmod2=lme(log10(MNL)~age+Timep+gen+Timep:age,random=~Timep|IDPatient,control=lmeC,method="ML",data=md2)
finalmod=lme(log10(MNL)~Timep+gen+age+gen:age,data=md2,random=~Timep|IDPatient,control=lmeC,method="ML")
summary(finalmod2)
anova(finalmod,finalmod2)##although keep lowest AIC
finalmod2
library(gam)
md2$Timep025=md2$Timep>=0.3
md3$Timep
dev.off()
par(xpd=TRUE)
library(gam)
fi=gam(log10(MNL)~bs(Timep,df=4,degree=3),data=md3,se=T)
dev.off()
plot(log10(md2$MNL)~md2$Timep,col=as.numeric(md2$gen),xlab="time post ops",ylab="log10 MNL",pch=3,cex=0.6,bty="n")
legend(26,0.8,txt,cex=0.6,col=1:5,pch=3,bty="n")
str(md2$gen)
legend(25,1.1,"Spline3",cex=0.9,col=2,pch=19,bty="n",xpd=TRUE)

points(fitted(fi)~md2$Timep,pch=20,col=2,cex=0.8)
txt=c("MEN","NF1","SDH","Sporadi","VHL")
plot(fi,se=T,lwd=2,ylab="Cubic Spline knot  at 4 quartile")
points(log10(md2$MNL)~md2$Timep,col=as.numeric(md2$gen),pch=3,cex=0.6,bty="n")

xnew=seq(0,4,length.out=200)
pred=predict(fi,data.frame(x=xnew))
length(pred)
plot(xnew,pred)



id=0.05
outlier=within(md2 ,{
  residp<-resid(finalmod,type="pearson")
  idx<-abs(residp>-qnorm(id/2))
})
 idx
out=subset(outlier,idx==1)
unique(out$IDPatient)#####Outlier id in mix0

9/89

md2[md2$IDPatient==outmixx$IDPatient,]##all men genes
head(ranef(finalmod),10)
KP=data.frame(ranef(finalmod))
par(mfrow=c(1,2))
plot(density(KP$Timep,xlim=c(-0.05,0.05),breaks=100),main="RANEF Time+")
plot(density(KP$X.Intercept.,xlim=c(-0.05,0.05),breaks=100),main="RANEF Intercept")
library(car)
library(effects)
plot(predictorEffect("Timep", finalmod),lines=list(multiline=TRUE, lwd=3), symbols=list(cex=1.5, pch=16),
     confint=list(style="none"), axes=list(y=list(lab="", type="response")),
      lattice=list(key.args=list(     cex=1, cex.title=1.2)),        grid=TRUE)





pairs(finalmod,~ranef(.),id=~IDPatient=="1499",main="RIS:Scatterplot of randeffect for 89ID ops1")
md2[md2$IDPatient=="5456",]##outliers all men


table(ID)
VarCorr(finalmodreml)
K=getVarCov(finalmod,type="marginal",individual="773")
round(cov2cor(K[[1]]),2)
intervals(finalmod,which="var-cov")
fix1=summary(finalmodreml)$tTable
fix1
printCoefmat(fix1,digits=3)
fix2=summary(finalmod)$tTable
fix2
library(stargazer)
stargazer(fix1,type="html",out="fix.doc")
printCoefmat(fix2)

lmIDwitjout1499
coef(lmIDwitjout1499)
coef(finalmod)

intervals(finalmodreml)
##predictions Fitted if nothing withMIXX RiS
dev.off()
plot(finalmodreml,log10(MNL)~fitted(.)|IDPatient,type="r",abline=c(0,1))#RI
#############PREDICTIONS FOR FINAL MODEL mixgenagenointer3and2############
fitted(finalmodreml)
ID=unique(md2$IDPatient)
N=length(ID)
N
MaxMTT=4.19
MaxNMNT=36.65
MaxMNT=13.45
MaxMTL=0.06
MaxNMNL=1.39
MaxMNL=0.85


fit=(fitted(finalmodreml))##id levels
  fit
fitpopu=fitted(finalmodreml,level=0)##group level
md2$fitfinal=fit
md2$fitpopu=fitpopu
pdf(file="C:/Users/mdy/Documents/buclin/predmictMNL.pdf")
for (i in 1:N){
  dat=subset(md2,IDPatient==ID[i])
  plot(dat$Timep,log10(dat$MNL),type="n",xlim=c(0,15),ylim = c(-1.60,0),xlab="time in Years after 1 ops",ylab="log10 MNL")
  title(paste0("Predic.MNL.ID",ID[i]))
  lines(dat$Timep,log10(dat$MNL),col=4,lty=2)
  points(dat$Timep,log10(dat$MNL),col=4,pch=19)
  points(dat$Timep,dat$fitfinal,pch=19,col=2)
  lines(dat$Timep,dat$fitfinal,col=2,lty=1)
  lines(dat$Timep,dat$fitpopu,col=3,lty=1,lwd=2)
  text(12,-0.75,labels="Observed",col=4)
  text(12,-1.0,labels="Predicted LMM",col=2)
  text(12,-1.3,labels="Pop fixef",col=3)
  abline(h=log10(MaxMNL),col=6,lty=2)
} 
dev.off()



############################################






beta+re
getwd()




mykey=list(tittle="ID",text=paste(levels(md2$IDPatient)))
MY=unique(md2$IDPatient)
MY
dev.off()
fit=predict(finalmodreml)
mixx=fit
xyplot(log10(md2$MNL)~md2$Timep|md2$IDPatient,layout=c(5,1),xlim=c(0,10), data=md2,fit=mixx,col="blue",type=c("p","r"),lwd=2,
      panel=function(x,y,...,fit,subscripts){
              panel.xyplot(x,y,...)
panel.text(8,1,labels=MY[panel.number()])
                             ypred=fitted(mixx)[subscripts]
      panel.lines(x,ypred,col=2,lwd=2,se.fit=T)
})

##make preddict function

##matric de predict for beta with lMEC
tpred=seq(0,20,len=1000)
tpred
Xpred=cbind(1,tpred)
str(Xpred)
md2$MNL[which(md2$MNL<=0.03)]<-0.03
md2$MNL
ID=unique(md2$IDPatient)
N=length(ID)
beta+re
getwd()
pdf(file="C:/Users/mdy/Documents/buclin/predmixx.pdf")
for (i in 1:N){
  dat=subset(md2,IDPatient==ID[i])
  bi=beta+re[i,]
  ypred=as.numeric(Xpred%*%bi)
  plot(md2$Timep,log10(md2$MNL),type="n",xlab="TIME IN Y since ops",ylab="Log10 MNL")
  title(paste0("ID",ID[i]))
  points(dat$Temps,log10(dat$MNL),pch=19)
  lines(tpred,ypred,col="red")
}
dev.off()














###SHRINKAGE

mean=with(md2,sapply(split(log10(MNL),IDPatient),mean))##mean IDPATIENT MNL

length(mean)
dev.off()
shrinakgemean=with(md2,sapply(split(fitted(finalmodreml),IDPatient),mean))
shrinakgemean
shri=data.frame(meanID=mean,fittedpred=shrinakgemean)
plot(shri$meanID,shri$fittedpred,main="comparison of mean /predicted value MIXX Shrinkage effect")##pas de valeurs tres shrinkee?
abline(0,1,col=2)
dev.off()
####NORMALITY OF RANEF####################################
cluster <- as.numeric(as.factor(md2$IDPatient))
cluster
age <- tapply(md2$age,cluster,unique)
hist(age)
tapply(md2$gen,cluster,unique)
gen <- levels(md2$gen)[tapply(md2$gen,cluster,unique)]
gen
table(gen)









##Correlation structur modif
 fm2Ovar.lme <- update( finalmod, correlation = corAR1() )##ou avec structure d'une variable
 plot(Variogram(fm2Ovar.lme),pch=19,form=~Timep|IDPatient,xlim=c(0,10),ylim=c(0,1.2))
 anova(finalmod,fm2Ovar.lme)
 
 


      






###avec censures Ã  gauche

cens <- (md2$MNL<=0.03)*1
cens
sum(cens)##32 valeurs MNL au seuil###############CENSOREING LMEC
cens*md2$IDPatient
getwd()
MNL003##censure ID des patients avec valeurs 0.03
MNL003=xtabs(~md2$IDPatient+cens)[,2]####list patient valeur MNL=00.3
MNL003
MNL003
ID[MNL003>=1]
write.csv(MNL003,file="MNL003.csv")


###############CENSOREING LMEC


yL <- md2$MNL 
yL
yL[cens==1]<- 0.03 
yL
X <- cbind(1,md2$Timep-0.25)
Z <- cbind(1,md2$Timep-0.25)
age <- tapply(md2$age,cluster,unique)
age
gen <- levels(md2$gen)[tapply(md2$gen,cluster,unique)]
gen
cluster <- as.numeric(as.factor(md2$IDPatient))
cluster

library(lmec)
lmm3 <- lmec(log10(yL),cens,X,Z,cluster)
beta <- lmm3$beta    # effets fixes
re <- t(lmm3$bi)    # effets alÃ©atoires
re
sigma <- lmm3$sigma  # stdev des rÃ©sidus
sigma
G <- lmm3$Psi        # covariance des effets alÃ©atoires
G
V.beta <- lmm3$varFix # covariance des effets fixes
se.beta <- sqrt(diag(V.beta)) # standard error des effets fixes
t.beta <- beta/se.beta # valeur t pour test de significativitÃ© des effets fixes
R.beta <- solve(diag(se.beta)%*%solve(V.beta)%*%diag(se.beta)) # corrÃ©lations entre effets fixes
R.beta
t.beta
2*pnorm(1.40,0,1,lower.tail = F)
# LL <- lmm3$loglik
# npar <- length(beta)+1+sum(lower.tri(G,diag=T))
# -2*LL+2*npar # AIC
lmm3$beta
# Voir effets alÃ©atoires en fonction de prÃ©dicteurs
age <- tapply(md2$age,cluster,unique)
gen <- levels(md2$gen)[tapply(md2$gen,cluster,unique)]
####REGROUPING gene in 3 classe Not desirable from Clinician CHUV
gen1 <- gen

gen1[which(gen%in%c("Men","NF1"))]<- "Men+NF1" 
gen1
gen1[which(gen%in%c("SDH","VHL"))] <- "SDH+VHL"
table(gen1)
str(gen1)
gen1=factor(gen1)

layout(matrix(1:(2*ncol(re)),ncol=ncol(re)),heights=rep(1,2),widths=rep(1,ncol(re)))
par(mar=c(3,3,2,0.5),mgp=c(1.8,0.6,0))
for(k in 1:ncol(re)){
  plot(age,re[,k],xlab="Age",ylab=paste("Random effect",k),main="RAnef versus Age and Genes")
  lines(20:80,predict(loess(re[,k]~age),20:80),col="red")
  boxplot(re[,k]~gen1,xlab="",ylab=paste("Random",k))
}
rep(1,nrow(re))
t(beta)[rep(1,nrow(re)),]
b <- t(beta)[rep(1,nrow(re)),]+re
b
layout(matrix(1:(2*ncol(b)),ncol=ncol(b)),heights=rep(1,2),widths=rep(1,ncol(b)))
par(mar=c(3,3,2,0.5),mgp=c(1.8,0.6,0))

for(k in 1:ncol(b)){
  plot(age,b[,k],xlab="Age",ylab=paste("predict",k))
  lines(20:80,predict(loess(b[,k]~age),20:80),col="red")
  boxplot(b[,k]~gen1,xlab="",ylab=paste("predict",k))
}

pente.pre <- b[,2]
pente.post <- b[,2]+b[,4]
plot(pente.pre,pente.post)
getwd()

###comparing 3 vs 5 hgene group nlme
summary(finalmodreml)
md2$gene1
md2$gene1
formula(finalmod)
md2$Timep
finalmodreml3gene=lme(log10(MNL)~Timep+gen1+age+gen1:age,data=md3,random=~Timep|IDPatient,control=lmeC,method="REML")
summary(finalmodreml3gene)
##design plot
##mean
plot.design(log10(md1$MNL)~md1$gen,main="mean of log10 MNL Overall time",cex.main=0.8,ylim=c(-1.2,-0.7),ylab="Log10 MNL")
plot.design(log10(new$MNL)~new$gen,main="mean of log10 MNL Postops",cex.main=0.8,ylim=c(-1.2,-0.7),ylab="Log10 MNL")##positive time
##var
plot.design(log10(md1$MNL)~md1$gen,fun=var,main="VAR of log10 MNL all time|ops=Yes",cex.main=0.8,ylim=c(0,1))
plot.design(log10(new$MNL)~new$gen,fun=var,main="VAR of log10 MNL all time|ops=Yes",cex.main=0.8,ylim=c(0,1))
new<-subset(md2,Timepc>=0)
new
md2
save.image()




































###interaction
str(md2)
attach(md2)
par(mfrow=c(1,2))
interaction.plot(cut(age,5),gen,log10(MNL),col=c(1,2,3,4,5),lwd=1.5,main="intercation age vs gen Overall | ops yes",xlab="classe age (5)",cex.main=0.8)##for preops ops post ops
##only with posp opst 3mois (Timeep)
md2$Timepos90
md222=subset(md2,md2$Timepos90==1)
md222
interaction.plot(cut(md222$age,5),md222$gen,log10(md222$MNL),col=c(1,2,3,4,5),lwd=1.5,,main="intercation age vs gen POST OPS | ops yes",xlab="classe age (5)",cex.main=0.8)##for preops ops post ops
##only with posp opst 3mois (Timeep))##for preops ops post ops


interaction.plot(cut(age,5),gen,log10(NMNL),col=c(1,2,3,4,5),lwd=1.5,main="intercation age vs gen Overall | ops yes",xlab="classe age (5)",cex.main=0.8)##for preops ops post ops
##only with posp opst 3mois (Timeep)
md2$Timepos90
md222=subset(md2,md2$Timepos90==1)
md2$MNT2014
dev.off()
##total biomarkers
dev.off()
md2$age5


####################INTERCATIONS
setwd("C:/Users/mdy/Documents/R")














#######################################################
################################################
#AVRIl 21

load(file="phe.Rda")
names(phe)
md1
md1$Time
md1$MNT2014
xyplot(log10(MNT2014)~Time,groups=IDPatient,data=md1,type="b",col=rainbow(220),xlim=c(-4,20))


####GROUPED DATA#######
##FAIRE PoUR Tiempos
id=sample(unique(md2$IDPatient),15)
id
data=subset(md2,md2$IDPatient%in%id)
names(data)
pat=groupedData(log10(MNL)~Time|IDPatient,data=data)
plot(pat,type=c("p","r"),lwd=2,xlim=c(0,10),col.line=2)



###mean on post preops
save(md2,file="md2.Rda")#"C:/Users/mdy/Documents/R"
save(md1,file="md1.Rda")
getwd()
setwd("C:~/R")
load("md2.Rdata")#"C:/Users/mdy/Documents/R"
load("md1.Rda")#C:/Users/mdy/Documents/R"
md1
md2$Time##is only positive time
md2$Timepos90*md2$Time
postopsmean=subset(md2,md2$Time>=0)
postopsmean
library(lattice)
dotplot(tapply(log10(postopsmean$MNL),postopsmean$gen,mean))
##looking for negative timne preops
md1$Timeneg##these are neg time
md1$Time
preopsmean=subset(md1,md1$Time<0)
preopsmean$Time##check onl neg time preops
par(mfrow=c(2,1))
x1=dotplot(tapply(log10(postopsmean$MNL),postopsmean$gen,mean),xlim=c(-1.5,1),xlab="mean log 10 MNL postops Time")
x2=dotplot(tapply(log10(preopsmean$MNL),preopsmean$gen,mean),xlim=c(-1.5,1),xlab="mean log 10 MNL preops Time")
install.packages("gridExtra")
library(gridExtra)
grid.arrange(x2,x1,ncol=2)

x1sd=dotplot(tapply(log10(postopsmean$MNL),postopsmean$gen,sd),xlim=c(0,1),xlab="sd log 10 MNL postops Time")
x2sd=dotplot(tapply(log10(preopsmean$MNL),preopsmean$gen,sd),xlim=c(0,1),xlab="sd log 10 MNL preops Time",abline(v=0))
grid.arrange(x2sd,x1sd,ncol=2)
designPlot(log10(md3$MNL),md3$gen,var)
##################INLUFENCE PLOT
install.packages("influence.ME")##extension of lme4
install.packages("lme4")
library(influence.ME)
library(lme4)
finalmodremlmer=lmer(log10(MNL)~Timep+gen+age+(Timep|IDPatient),data=md2)
md2$Timep
finalmodremlmer
formula(finalmodremlmer)
influ=influence(finalmodremlmer,obs=TRUE)
treillis.colors=1
plot(influ,which="cook",cutoff = 4/89,cex=0.5,cex.y=0.3,sort=TRUE,scales=list(cex=0.5,alternating=1),xlab="cook.distance",ylab="IDPatient")
cook=cooks.distance(influ)
dim(cook)
V=which(cook>=0.05)
cook[V,]
