rm(list = ls())

#file from KARIME CALLED FVS one copy original kept on mail
#rename on C/R// phep.csv to avoid confusion with preceeding work scrap
#the on xls delted the gras value from pre ops + format unilangage factor category make phenopreops
#MALIN exclus
phe=read.csv(file="Phekarimnonmalintempsarr.csv",header=T,sep=";")

levels(phe$Malignite)
dim(phe)
phe=phe[phe$Malignite == "" | phe$Malignite == "Non malin", ]#MALIN exclus
dim(phe)
1405-1581##176 measures enlevées de malins

str(phe$DNPatient)
factor(phe$IDPatient)
####QC CHECK AS BEFORE 220 levels of ID
library(lubridate)
born=as.Date.factor(phe$DNPatient,format="%d.%m.%Y")
born
phe$born=born
phe$born
library(lubridate)
today()
phe$age=(today()-born)/365
phe$age=round(phe$age,0)
phe$age
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
2452+7320
9772/24

phe$IDPatient=factor(phe$IDPatient)
phe$IDPatient##220 patients????? to confirm [including No malignant]



str(phe$MaladieGent)
levels(phe$MaladieGent)
phe$gen=phe$MaladieGent
#RECODING SEXE INTO 2 CAT with2 ? as women
sum(phe$gen=="")##je ne sais pas d'ou viens ce gene"" je recherche dans excel Pas trouvé 1 point gen/1405?
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
sapply(phe,class)
phe$gen=factor(phe$gene)

DW#codin arrival
phe$DateArr
phe$arrival=as.Date.factor(phe$DateArr,format="%d.%m.%Y")
str(phe$arrival)
summary(phe$TempsArr)
which(is.na(phe$TempsArr))
phe[1405,]
phe

###URL BIOMARKERS
MaxMTT=4.19
MaxNMNT=36.65
MaxMNT=13.45
MaxMTL=0.06
MaxNMNL=1.39
MaxMNL=0.85
phe$MaxMTT=rep(MaxMTT,1405)
phe$MaxNMNT=rep(MaxNMNT,1405)
phe$MaxMNT=rep(MaxMNT,1405)
phe$MaxMTL=rep(MaxMTL,1405)
phe$MaxNMNL=rep(MaxNMNL,1405)
phe$MaxMNL=rep(MaxMNL)
phe[1,]
MTT2014=ifelse(phe$IDDemande<=18100,1.4*phe$MTT,phe$MTT)
MTT2014
MNT2014=ifelse(phe$IDDemande<=18100,1.4*phe$MNT,phe$MNT)
MNT2014
NMNT2014=ifelse(phe$IDDemande<=18100,1.4*phe$NMNT,phe$NMNT)
str(MNT2014)
phe$MTT2014=MTT2014
phe$NMNT2014=NMNT2014
phe$MNT2014=MNT2014
phe[51,]
##???verif ligne 51 MAxMNT et MTT par exmple
6.18*1.40
Age
Age=as.numeric(phe$age)
AGEF=cut(Age,5)
AGEF
phe$AGE5=AGEF
phe[51,]
names(phe)
phemdy1=phe[ ,c(1,2,3,29,5,8,9,10,11,12,37,38,39,16,17,27,31,32,33,34,35,36,40)]
phemdy1
names(phemdy1)
summary(phemdy1)
dim(phemdy1)

write.csv(phemdy1,"C:/R/phemdy1.csv", row.names = TRUE) 
