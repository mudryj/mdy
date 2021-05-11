getwd()
setwd("C:/Users/mudry/OneDrive/Documents/UNINE/Buclin/DATA")
load("phe.RData")
sum(is.na(phe))
dim(phe)
table(phe$operation.y.n)
#to answer the question of covariate and to be operated YES NO a logistic regression can be runned after some coding and check separability
str(phe$operation.y.n)
levels(phe$operation.y.n)
phelogistic=phe[phe$operation.y.n=="y" | phe$operation.y.n=="n"|phe$operation.y.n=="Y",]###out of 220 patients patients with "" /? removed
 phe$operation.y.n
 
phelogistic
phelogistic$opsnew[phelogistic$operation.y.n=="y"]<-"yes"
phelogistic$opsnew[phelogistic$operation.y.n=="Y"]<-"yes"
phelogistic$opsnew[phelogistic$operation.y.n=="n"]<-"no"
is.na(phelogistic$opsnew)
str(phelogistic$opsnew)
phelogistic$opsnew=factor(phelogistic$opsnew)
str(phelogistic$opsnew)
unique(phelogistic$opsnew)
unique(phe$IDPatient)
phelogistic$opsnewyes=ifelse(phelogistic$opsnew=="yes",1,0)
table(phelogistic$opsnew)##operer 898 out of 1405 lines
length(unique(phelogistic$IDPatient))



droplevels(unique(phelogistic$IDPatient))
opere=phelogistic[phelogistic$opsnewyes==1,]
dim(opere)
opere
ID2=unique(opere$IDPatient)
droplevels(ID2)
187-114
sum(!is.na(ID2))##114 patients opere 114 73 non opere
table(opere$Sexi[ID2])#amongst 59M 55 W
#plot with lm!!!
sex=opere$Sexi
sex
opere
save(opere,file ="oper.RData")

pli=load("oper.RData")
pli


plot(jitter(as.numeric(phelogistic$age),3),as.numeric(phelogistic$opsnewyes),pch=19,ylab="ops statut 1 = surgery",xlab="age in Years")#seprability ok in all points weighting operation yes and no
abline(lm(phelogistic$opsnewyes~phelogistic$age))
##to evalkuate effect of covoariates
phelogistic$opsnewyes
logimodel=glm(opsnewyes~Sexi+gen+age,family="binomial",data=phelogistic)
summary(logimodel)
vcov(logimodel)#SE check
##check IT devaiance accrd df X hisqre GOFIT ok SE error OK coef ok
anova(logimodel,type="F")
step(logimodel)##rmoving age dont change much coefficient
new=update(logimodel,.~.-age)
new
plot(logimodel)##GOFIT

##plot randomized res.

library(statmod)
qqnorm(qresid(new))##normality asusmed and no overdisp in other plot
plot(new$fitted.values~phelogistic$age)

#contsrat are treatment made for MEN Gene / gender men 
contrasts(phelogistic$gen)
#relevel for check
phelogistic$gen
genespor <- relevel(phelogistic$gen, "sporadi")
contrasts(genespor)
new2=update(new,.~.-gen+genespor+age)#in accordance tpo 1st contrast men
summary(new2)
new3=update(new,.~.-gen+genespor)
##R contrast give you a group comparison which is good in ORR
exp(coef(new))

is.singular(new)#no matrix defficiency or ill cond.
#modelling patient with operation will reduced the data in hand N
#we cannot conclude is age is able to predict biomarker we just say that in hand age is not significant for surgery in our datasample but gender and some gene.This is reversed by what found in LMM but they answwer two different questions
save(phelogistic,  file = "phelogistic.RData")
##reformatting 

##################we need a logistic with 1 obs measurement per ID and only one
names(phelogistic)
logi1=phelogistic[ ,c(1,10,27,28,30,33)]
names(logi1)
logi1
ID=unique(logi1$IDPatient)
names(logi1)##187 yes and no ops
##need to reshape in wide format having 1 statitic predictor point for each ID otherwise flaws in long format

data_reshape1 <- reshape(logi1,                               
                         idvar = c("gene","IDPatient","Sexi","opsnewyes","age"),                         timevar = "MNL",
                         direction = "wide")
dim(data_reshape1)##ok days=years
data_reshape1
data_reshape1$age=as.numeric(data_reshape1$age)
data_reshape1$gene=factor(data_reshape1$gene)
data_reshape1$Sexi=factor(data_reshape1$Sexi)
save(data_reshape1,file="wideformat.Rdata")
phewide=data_reshape1
plot(phewide$age,phewide$opsnew,pch=19)
summary(logimodel)
##univariate logistic no contrast spec and zero age no inter:
logimodel1=glm(opsnewyes~age+Sexi+gene,family="binomial",data=phewide)
summary(logimodel1)##treatm contrast men
exp(coef(logimodel1))##age is not significant so probably no intere
#3univariate logistic nointercation with centered age
logimodel2=glm(opsnewyes~I(age-48)+Sexi+gene,family="binomial",data=phewide)
summary(logimodel2)##treatm contrast men
exp(coef(logimodel2))##age is not significant so probably no intere
vcov(logimodel1)#SE check
logimodel3=glm(opsnewyes~age+Sexi+gene+age:gene,family="binomial",data=phewide)
dim(phewide)
summary(logimodel3)##treatm contrast men
exp(coef(logimodel3))#overdispersion occurs
250/186
plot(logimodel3)
par(mfrow=c(1,1))
plot(rstandard(logimodel3))
K=influence.measures(logimodel3)
plot(cooks.distance(logimodel3),pch=19)
identify(cooks.distance(logimodel3),pch=19) 
li=phewide[71,]
phewide2=phewide[c(-71,-95,-132,-170),]
phewide
phewide2##obs 743 removed
logimodel3.743=glm(opsnewyes~age+Sexi+gene+age:gene,family="binomial",data=phewide2)
summary(logimodel3.743)
##removing outliers doesnt cure the probleme


vcov(logimodel1)#SE check
##removing this patient
phewide1=subset(phewide,select==(li))
##check IT devaiance accrd df X hisqre GOFIT ok SE error OK coef ok
anova(logimodel1,type="Chisq")###gene take out major deviance 2xage
step(logimodel1)##rmoving age dont change much coefficient
new=update(logimodel,.~.-age)
summary(new)

anova(new)
plot(logimodel)##GOFIT
##plot randomized res.

library(statmod)
qqnorm(qresid(new),id=0.05,sub="randomized residuals")##normality asusmed and no overdisp in other plot
qqline(qresid(new))

anova(int9,type="Chisq")

##FINAL MODEL 
summary(new)
exp(coef(new))


###checking for intercations
##run univariate check if IR chage when adding a cov
int1=glm(opsnewyes~Sexi,family="binomial",data=phewide)
summary(int1)
options("contrasts")
int2=glm(opsnewyes~Sexi+gene,family="binomial",data=phewide)
summary(int1)
summary(int2)##gene has changed so interca sex and gene possibily for operations
confint(int1)
confint(int2)
##Ok possibliby a two wa int gene factor (5 levels) with gender(2 levels)

int3=glm(opsnewyes~gene,family="binomial",data=phewide)
int4=glm(opsnewyes~Sexi+gene,family="binomial",data=phewide)
summary(int3)
summary(int4)##dont chnage as much as other
exp(coef(int1))
exp(coef(int2))

confint(int1)
confint(int2)##well the loggosd confint is not very different so inter is light
##age is continuous cut age in 5 see OR in slope if changing
agecut=cut(phewide$age,5)#factor 5 leevels
agecut
int5=glm(opsnewyes~age,family="binomial",data=phewide)
int6=glm(opsnewyes~agecut,family="binomial",data=phewide)
summary(int5)
exp(coef(int5))
summary(int6)
exp(coef(int6))###of course OR is not constant through age for surgery
int7=glm(opsnewyes~Sexi+agecut,family="binomial",data=phewide)
summary(int7)
exp(coef(int7))##no change in coef age cut no inter age/sex
int8=glm(opsnewyes~gene+agecut,family="binomial",data=phewide)
summary(int8)
exp(coef(int8))##il y a clairement une inte gene age  operation
##ecideemnt les sporadi sont opere en urgence et vieux
##men corresponde a une genetique age dependante.
##can be turned into reserach questions
int10=glm(opsnewyes~age+gene,family="binomial",data=phewide)
summary(int9)
int10=glm(opsnewyes~age+gene+age:gene,family="binomial",data=phewide)
summary(int10) 
