rm(list = ls())
#master project
#1st data delivered for EDA statistics
#file from KARIME CALLED FVS one copy original kept on mail
#rename on C/R// phep.csv to avoid confusion with preceeding work scrap
#the on xls delted the gras value from pre ops + format unilangage factor category make phenopreops
#MALIN exclus
#format long
getwd()
setwd("C:/Users/mudry/OneDrive/Documents/UNINE/Buclin/DATA")
phe=read.csv(file="Phekarimnonmalintempsarr.csv",header=T,sep=";")
dim(phe)
str(phe)
phe
#removing malinant form ID
levels(phe$Malignite)
dim(phe)
phe=phe[phe$Malignite == "" | phe$Malignite == "Non malin", ]#MALIN exclus
dim(phe)
1405-1581##176 measures enlevées de malins
1405*25 #pixels
1405*6#biomarkers
#8430 mesures d ebiomarkeurs
str(phe$DNPatient)


#date as date variable ts : AGE and OPERATION converted
str(phe$IDPatient)
factor(phe$IDPatient)
unique(phe$IDPatient)
####QC XLS CHECK AS BEFORE 220 levels of ID
str(phe$DNPatient)
library(lubridate)
born=as.Date.factor(phe$DNPatient,format="%d.%m.%Y")
born
boxplot(born,main="age distribution N =220 ")
phe$born=born
phe$born
phe$age=(today()-born)/365
phe$age=round(phe$age,0)
median(as.numeric(phe$age),na.rm=TRUE)
mean(as.numeric(phe$age),na.rm=TRUE)
##plot age
hist(as.numeric(phe$age),breaks=15,sub=  "N =220 ",col=8,xlab="Age",main="",cex=1.4)
abline(v=48.25,col=3)
abline(v=50,col=2)
legend("topright",legend=c("mean 48.25","median 50"),col=c(3,2),lwd=1)

#gender coding
str(phe)
str(phe$Sexe)
levels(phe$Sexe)
#RECODING SEXE INTO 2 CAT with2 ? as women
phe$Sexi[phe$Sexe=="?"]<-"w"
phe$Sexi[phe$Sexe=="w"]<-"w"
phe$Sexi[phe$Sexe=="f"]<-"w"
phe$Sexi[phe$Sexe=="F"]<-"w"
phe$Sexi[phe$Sexe=="m"]<-"m"
phe$Sexi[phe$Sexe=="M"]<-"m"
str(phe$Sexi)
phe$Sexe=factor(phe$Sexi)
table(phe$Sexe)
ID=unique(phe$IDPatient)
table((phe$Sexi[ID]))##more men but more measrues fo women
#recheck 
unique(phe$Sexi)
table(phe$Sexi[unique(phe$IDPatient)])

##age bs Sexe dist
barplot(table(phe$Sexe,phe$age),main="measurements distribution of biomarker vs age
        N= 220 Patients",xlab="Age",ylab="count",legend=TRUE,col=c(8,6),cex.axis=0.8,cex=0.7)
#why women have more measurments in cohort is it control and understood?

#gene coding 
str(phe$MaladieGent)
levels(phe$MaladieGent)
phe$gen=phe$MaladieGent
#RECODING SEXE INTO 2 CAT with2 ? as women
phe$gen
phe$gene[phe$gen=="Men1"]<-"Men"
phe$gene[phe$gen=="Men2"]<-"Men"
phe$gene[phe$gen=="Men2a"]<-"Men"
phe$gene[phe$gen=="SDHD"]<-"SDH"
phe$gene[phe$gen=="SDHB"]<-"SDH"
phe$gene[phe$gen=="SDHx"]<-"SDH"
phe$gene[phe$gen=="VHL"]<-"VHL"
phe$gene[phe$gen=="NF1"]<-"NF1"
phe$gene[phe$gen=="Testé sans mutations connues"]<-"sporadi"
phe$gene[phe$gen=="Non testé"]<-"sporadi"

phe$gen=factor(phe$gene)
phe$gen
table((phe$gen[ID]))##n count per genes category
names(phe)
biomarkers=phe[,c(9:15)]
biomarkers
par(mfrow=c(2,3))
for (i in 1:6)

boxplot(log10(biomarkers[,i])~gen,data=phe,ylim=c(-2.3,3),col=8,main=paste("log10",colnames(biomarkers[i])))

for (i in 1:6)
boxplot(log10(biomarkers[,i])~gen,data=phe,ylim=c(-2.3,3),col=8,main=paste("log10",colnames(biomarkers[i])))
library(lattice)
#boxplot Biomarker|genes
bw.theme <- trellis.par.get()
bw.theme$box.dot$pch <- "|"
bw.theme$box.rectangle$col <- "black"
bw.theme$box.rectangle$lwd <- 2
bw.theme$box.rectangle$fill <- "grey90"
bw.theme$box.umbrella$lty <- 1
bw.theme$box.umbrella$col <- "black"
bw.theme$plot.symbol$col <- "grey40"
bw.theme$plot.symbol$pch <- 19
bw.theme$plot.symbol$cex <- 0.7
bw.theme$strip.background$col <- "grey80"
bwplot(log10(phe$MNL)~phe$Sexi|phe$gen,auto.key=TRUE, par.settings=bw.theme)

#othr type
barchart(xtabs(~phe$gen+phe$Sexe),sub="measures count: partition of Genes vs Sexe",auto.key = list(space = "right" ))
##nbs of measurements|ID

barplot(D,col=8,ylab="nbs of measurments",ylim=c(0,20),sub="ID",cex.lab = 0.8)
par(mfrow=c(1,1))
##request by Buclin all by ID chart|nbs measure
library(lattice)
barchart(D[1:30],data=phe,main="nbrs of measures per patients")
barchart(D[31:60],data=phe,main="nbrs of measures per patients")
barchart(D[61:90],data=phe,main="nbrs of measures per patients")
barchart(D[91:120],data=phe,main="nbrs of measures per patients")
barchart(D[121:150],data=phe,main="nbrs of measures per patients")
barchart(D[151:180],data=phe,main="nbrs of measures per patients")
barchart(D[181:210],data=phe,main="nbrs of measures per patients")
barchart(D[211:220],data=phe,main="nbrs of measures per patients")
length(D)###Rechechek nombres total de patient

######
#codin arrival
phe$DateArr
phe$arrival=as.Date.factor(phe$DateArr,format="%d.%m.%Y")
str(phe$arrival)
str(phe)
phe$arrival
summary(phe$TempsArr)
which(is.na(phe$TempsArr))#pas de NA


###URL BIOMARKERS
MaxMTT=4.19
MaxNMNT=36.65
MaxMNT=13.45
MaxMTL=0.06
MaxNMNL=1.39
MaxMNL=0.85

#log10 for use
MaxMTT=log10(4.19)
MaxNMNT=log10(36.65)
MaxMNT=log10(13.45)
MaxMTL=log10(0.06)
MaxNMNL=log10(1.39)
MaxMNL=log10(0.85)

#######################
#adjusted value for value bebore 2014 linear factor 1.4
phe[1,]
MTT2014=ifelse(phe$IDDemande<=18100,1.4*phe$MTT,phe$MTT)
MTT2014
MNT2014=ifelse(phe$IDDemande<=18100,1.4*phe$MNT,phe$MNT)
MNT2014
NMNT2014=ifelse(phe$IDDemande<=18100,1.4*phe$NMNT,phe$NMNT)
str(MNT2014)

##profiles plot (Spaghettis plots ID trajectrory
xyplot(log(phe$MTL)~phe$TempsArr|phe$gen,groups=factor(phe$IDPatient), type="l",as.table = TRUE,aspect=2,par.settings=simpleTheme(col="red", col.line="grey"),main="log MTL all Patient ID/gene with pre ops no malin",panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
  panel.abline(h = log(MaxMTL), lty = 2,col=2) 
  cof <- lm(y~x)
  panel.abline(reg = cof,lwd=2, col="blue", ...)
})

xyplot(log(phe$MNL)~phe$TempsArr|phe$gen,groups=factor(phe$IDPatient), type="l",as.table = TRUE,aspect=2,par.settings=simpleTheme(col="red", col.line="grey"),main="log MNL all Patient ID/gene with pre ops no malin",panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
  panel.abline(h = log(MaxMNL), lty = 2,col=2) 
  cof <- lm(y~x)
  panel.abline(reg = cof,lwd=2, col="blue", ...)
})
##rechute après 6 ans
xyplot(log(phe$NMNL)~phe$TempsArr|phe$gen,groups=factor(phe$IDPatient), type="l",as.table = TRUE,aspect=2,par.settings=simpleTheme(col="red", col.line="grey"),main="log NMNL all Patient ID/gene with pre ops no malin",panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
  panel.abline(h = log(MaxNMNL), lty = 2,col=2) 
  cof <- lm(y~x)
  panel.abline(reg = cof,lwd=2, col="blue", ...)
})
#rechute 6 ans enlevé une année pour récup : pas assez de rechute pour y voir claire
#####TOTALES

##############
xyplot(log(MNT2014)~phe$TempsArr|phe$gen,groups=factor(phe$IDPatient), type="l",as.table = TRUE,aspect=2,par.settings=simpleTheme(col="red", col.line="grey"),main="log MNT all Patient ID/gene with pre ops no malin",sub="Green = Mean, Red=url,Blue=regression",panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
  panel.abline(h = log(MaxMNT), lwd=1,lty = 2,col=2)
  panel.abline(h = 1.50, lwd=1,lty = 2,col=3)
  cof <- lm(y~x)
  panel.abline(reg = cof,lwd=2, col="blue", ...)
  panel.loess(x, y, span = 2/3, degree = 2,col=6,lwd=2,...)
})

mean(log(MNT2014),na.rm=T)
#rechute VHL 4ans
xyplot(log(MTT2014)~phe$TempsArr|phe$gen,groups=factor(phe$IDPatient), type="l",as.table = TRUE,aspect=2,par.settings=simpleTheme(col="red", col.line="grey"),main="log MTT all Patient ID/gene with pre ops no malin",sub="Max MTT 4.19, Green=Mean, Red=url,Blue=regression",panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
  panel.abline(h = log(MaxMTT), lwd=1,lty = 2,col=2) 
  panel.abline(h = 1.61, lwd=2.5,lty = 3,col=3)
  cof <- lm(y~x)
  panel.abline(reg = cof,lwd=1.3, col="blue", ...)
  panel.loess(x, y, span = 2/3, degree = 2,col=6,lwd=2,...)
})
mean(log(MTT2014),na.rm=T)
###MTT un url correct? All mean are above the url Overdiagnostic
MaxMTT
###MTT is disciminant dominant in pathology
########################NMNT
mean(log(NMNT2014),na.rm=T)

xyplot(log(NMNT2014)~phe$TempsArr|phe$gen,groups=factor(phe$IDPatient), type="l",as.table = TRUE,aspect=2,par.settings=simpleTheme(col="grey", col.line="grey"),main="log NMNT all Patient ID/gene with pre ops no malin",sub="mean = green,red = url,blue= Regression",panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
  panel.abline(h = log(MaxNMNT), lwd=1.2,lty = 2,col=2) 
  panel.abline(h = 2.61, lwd=2.5,lty = 3,col=3)
  cof <- lm(y~x)
  panel.abline(reg = cof,lwd=2, col="blue", ...)
  panel.loess(x, y, span = 2/3, degree = 2,col=6,lwd=2,...)
})



#men ref group
##all in general attracted to mean 
bwtheme  <- canonical.theme()
bwplot(log(NMNT2014)~phe$gen,main="boxplot of NMNT [red url]",
       panel=function(...) {
         panel.abline(h=log(MaxNMNT), col="red")
         panel.bwplot(...)
       })
bwplot(log(MTT2014)~phe$gen,main="boxplot of MTT [red url]",
       panel=function(...) {
         panel.abline(h=log(MaxMTT), col="red")
         panel.bwplot(...)
       })
bwplot(log(MNT2014)~phe$gen,main="boxplot of MNT [red url]",
       panel=function(...) {
         panel.abline(h=log(MaxMNT), col="red")
         panel.bwplot(...)
       })
mycol <- rgb(220, 60, 155, max = 255, alpha = 100, names = "red")
boxplot(log(phe$MNL)~phe$gen,main="Box-Strip for Log MNL on 1405 measures| Genes ",cex.main=0.8,col="lightgrey",outline=F)
stripchart(log(phe$MNL)~phe$gen,method="jitter",vertical=TRUE,cex = 1,pch=19,col=mycol,add=TRUE)
##NMNL
boxplot(log(phe$NMNL)~phe$gen,main="Box-Strip for Log NMNL on 1405 measures| Genes ",cex.main=0.8,col="lightgrey",outline=F)
stripchart(log(phe$NMNL)~phe$gen,method="jitter",vertical=TRUE,cex = 1,pch=19,col=mycol,add=TRUE)
##MTL
boxplot(log(phe$MTL)~phe$gen,main="Box-Strip for Log MTL on 1405 measures| Genes ",cex.main=0.8,col="lightgrey",outline=F)
stripchart(log(phe$MTL)~phe$gen,method="jitter",vertical=TRUE,cex = 1,pch=19,col=mycol,add=TRUE)
###MTN
boxplot(log(phe$MNT2014)~phe$gen,main="Box-Strip for Log MNT on 1405 measures| Genes ",cex.main=0.8,col="lightgrey",outline=F)
stripchart(log(phe$MNT2014)~phe$gen,method="jitter",vertical=TRUE,cex = 1,pch=19,col=mycol,add=TRUE)
#####NMNT
boxplot(log(phe$NMNT2014)~phe$gen,main="Box-Strip for Log NMNT on 1405 measures| Genes ",cex.main=0.8,col="lightgrey",outline=F)
stripchart(log(phe$NMNT2014)~phe$gen,method="jitter",vertical=TRUE,cex = 1,pch=19,col=mycol,add=TRUE)
##MTT
boxplot(log(phe$MTT2014)~phe$gen,main="Box-Strip for Log MTT on 1405 measures| Genes ",cex.main=0.8,col="lightgrey",outline=F)
stripchart(log(phe$MTT2014)~phe$gen,method="jitter",vertical=TRUE,cex = 1,pch=19,col=mycol,add=TRUE)

###les vhl ont des valeurs basales plus faibles d'ou l'importance du diagnostic 

##biomarkers histogram
par(mfrow=c(2,3))
for (i in 1:length(biomarkers)) { 
hist(log10(biomarkers[,i]),freq=FALSE,ylim=c(0,2.3),col="grey83",xlab="",main=paste("log10",colnames(biomarkers[i])))
lines(density(log10(biomarkers[,i])),col=2)
}
##skewness of log10BIO
library(psych)
b=NULL
bio=na.omit(biomarkers)
for (i in 1:6){
b[i]=round(skewness(log10(bio[ ,i])),2)
  }
b

biomarkers=data.frame(biomarkers)
phe$gene
pairs(log10(biomarkers),lower.panel = panel.smooth,pch=21,bg = rainbow(2)[V],na.action=na.omit)

install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(log10(biomarkers), histogram = TRUE, method = "pearson",sub="Pearson correlation matrix",pch=19)
par(mfrow=c(1,1))
##publication Skawa age vs biomarker trend
plot(log10(phe$NMNL)~phe$age,pch=19,col="grey53",xlab="Age",ylab="log10NMNL")
text(90,1, "slope=1.0052",col=2)
text(90,0.7, "pval= 0.00",col=2)

abline(lm(log10(phe$NMNL)~phe$age),col=2)
exp(coef(lm(log10(phe$NMNL)~phe$age)))
summary(lm(log10(phe$NMNL)~phe$age))
####verifictaion with boxplot vs age
AGEF=cut(Age,5)
AGEF

par(mfrow=c(1,3))
boxplot(log(phe$MNL)~AGEF,col="grey",col.pts="red")
boxplot(log(phe$MTL)~AGEF,col="grey")
boxplot(log(phe$NMNL)~AGEF,col="grey")
phe$gen
par(mfrow=c(1,1))
bwplot(log(phe$NMNL)~AGEF|phe$gen,main="NMNL vs Age|gene no outlier",xlab=list("Age de 6 a 96 ans",cex=0.7),ylab="Log NMNL",scales=list(cex=0.7),col="grey",auto.key = TRUE,box.width=0.7,fill="grey",box.border=FALSE,varwidth=F,par.settings = list(plot.symbol = list(pch = "")))



save(data1, file = "data.RData")
# Save multiple objects
save(phe,  file = "phe.RData")



##profile plot for MNT log10 with time of arrival ##Jm note use phe1 data from TempsBuclin.R whre phe1 is coded if needed
##msc chapter 2.1.2 for 2 biomarker add + in outcome var
xyplot(log10(phe1$MNT)+log10(phe1$NMNT)~phe1$TempsArr, groups=as.numeric(phe1$IDPatient), data=phe, type="b",ylab="Log10 MNT",xlab="Year from Arrival Time points",xlim=c(0,8),lwd=3,cex=0.8,col.line=rainbow(30),auto.key = FALSE,par.settings = list(superpose.symbol = list(pch = 19,fill=heat.colors(30))),panel = function(...) {
  panel.abline(h = 1.128, lty = 6,lwd=1,col=2)
  panel.text(7,1.3, label="URL MNT", font=2,col=2)
  panel.xyplot(...)
})

# To load the data again
getwd()
setwd("C:/Users/mudry/OneDrive/Documents/UNINE/Buclin/DATA")
load("phe.RData")
phe

##missing vale
MISS=xtabs(is.na(phe$MNL)+is.na(phe$MNT)+is.na(phe$NMNL)+is.na(phe$NMNT)+is.na(phe$MTT)+is.na(phe$MTL)~phe$Sexi+phe$gene)
margin.table(MISS,1)
########################TEST TO BE REMOVED
#############################
#########################
table(is.na(phe))
4899+38656
