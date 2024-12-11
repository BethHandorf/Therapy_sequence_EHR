################################################################################
# This file creates a synthetic dataset for subsequent microsimulations ########
# Transition times are simulated according to 1 of 9 different models ########## 
# Synthetic covariates and model parameters are read in from external files ####
# This code also calculates the "true" costs and effects if everyone 33 ########
# receives "sequence 1" or if everyone receives "sequence 2"####################


library(copula)
library(MASS)
library(psych)
library(readxl)
library(dplyr)
library(flexsurv)

set.seed(100)



rm(list=ls())


#### Different output types 

#out.type<-"ind"            #independence
out.type<-"clayton"         #clayton copula, weibull model
#out.type<-"GaussCS"        #Gauss compound symmetric copula, weibull model
#out.type<-"GaussUN"        #Gauss unstructured copula, weibull model
#out.type<-"TCS"            #T compound symmetric copula, weibull model
#out.type<-"TUN"            #T unstructured copula, weibull model
#out.type<-"smallSS"        #clayton copula, weibull model
#out.type<-"LogLogisClayton" #claton copula, log logistic model
#out.type<-"lnormClayton"   #clayton copula, log normal model

#Read in synthetic  covariate data
covariates<-read.csv("sim_dataset.csv")

#just take first 2000
covariates<-covariates[1:2000,]
N=2000

#for small SS, just take the first 500
if(out.type=="smallSS"){
  covariates<-covariates[1:500,]
  N=500}

#Recode Performance status covariate
covariates$ECOG1<-as.numeric(covariates$EcogValue==1)
covariates$ECOG2<-as.numeric(covariates$EcogValue>1)
#matrix of model covariates
model.covar<-covariates %>% dplyr:::select(Surgery, AgeAdvDx, ECOG1,ECOG2,albumin) 


weibull.coef<-read_xlsx("reg_param_est.xlsx",sheet = "Weibull")
llogis.coef<-read_xlsx("reg_param_est.xlsx",sheet = "llogis")
lnorm.coef<-read_xlsx("reg_param_est.xlsx",sheet = "lnorm")


#trt = 
#simulate treatments based on pre-defined propensity score
trt<-rep(NA, dim(covariates)[1])
for (i in 1:length(trt)) {
  trt[i]<- rbinom(1,1,prob=covariates$p.trt.CisL1[i])
}


###############################################################################
############ Set up regression parameters #####################################
###############################################################################

######## Set up regression parameters ###################

mod.coef<-weibull.coef
if(out.type=="LogLogisClayton"){mod.coef<-llogis.coef}

if(out.type!="lnormClayton"){
  ######Progression from L1 
  #coefficients
  beta.trt.progL1<-as.numeric(mod.coef[3,2]) #Weibull coef are already reversed
  beta.coef.progL1<-mod.coef[c(4:7,9),2] # No need for ECOG NA
  #effect of covariates
  XB.progL1<-as.matrix(model.covar)%*%as.matrix(beta.coef.progL1)
  #weibull parameters
  shape.progL1=as.numeric(mod.coef[1,2])
  scale.progL1=as.numeric(mod.coef[2,2])
  
  #scale for weibull regression
  scale.reg.progL1<-scale.progL1*(exp(-1*shape.progL1*(beta.trt.progL1*trt+XB.progL1)))^(-1/shape.progL1)
  scale.reg.progL1.alltreat<-scale.progL1*(exp(-1*shape.progL1*(beta.trt.progL1+XB.progL1)))^(-1/shape.progL1)
  scale.reg.progL1.notreat<-scale.progL1*(exp(-1*shape.progL1*(XB.progL1)))^(-1/shape.progL1)
  
  
  ######Death from L1 
  #coefficients
  beta.trt.deathL1<-as.numeric(mod.coef[3,3]) #Weibull coef are already reversed
  beta.coef.deathL1<-mod.coef[c(4:7,9),3]
  #effect of covariates
  XB.deathL1<-as.matrix(model.covar)%*%as.matrix(beta.coef.deathL1)
  #weibull parameters
  shape.deathL1=as.numeric(mod.coef[1,3])
  scale.deathL1=as.numeric(mod.coef[2,3])
  
  #scale for weibull regression
  scale.reg.deathL1<-scale.deathL1*(exp(-1*shape.deathL1*(beta.trt.deathL1*trt+XB.deathL1)))^(-1/shape.deathL1)
  scale.reg.deathL1.alltreat<-scale.deathL1*(exp(-1*shape.deathL1*(beta.trt.deathL1+XB.deathL1)))^(-1/shape.deathL1)
  scale.reg.deathL1.notreat<-scale.deathL1*(exp(-1*shape.deathL1*(XB.deathL1)))^(-1/shape.deathL1)
  
  ######Progression from L2 
  #coefficients
  beta.trt.progL2<-as.numeric(mod.coef[3,4]) #Weibull coef are already reversed
  beta.coef.progL2<-mod.coef[c(4:7,9),4] 
  #effect of covariates
  XB.progL2<-as.matrix(model.covar)%*%as.matrix(beta.coef.progL2)
  #weibull parameters
  shape.progL2=as.numeric(mod.coef[1,4])
  scale.progL2=as.numeric(mod.coef[2,4])
  
  
  #get rid of the effect of treatment - too close to zero
  scale.reg.progL2<-scale.progL2*(exp(shape.progL2*-1*(XB.progL2)))^(-1/shape.progL2)
  scale.reg.progL2.alltreat<-scale.progL2*(exp(shape.progL2*-1*(XB.progL2)))^(-1/shape.progL2)
  scale.reg.progL2.notreat<-scale.progL2*(exp(shape.progL2*-1*(XB.progL2)))^(-1/shape.progL2)
  
  ######Death from L2 
  #coefficients
  beta.trt.deathL2<-as.numeric(mod.coef[3,5]) #Weibull coef are already reversed
  beta.coef.deathL2<-mod.coef[c(4:7,9),5]
  #effect of covariates
  XB.deathL2<-as.matrix(model.covar)%*%as.matrix(beta.coef.deathL2)
  #weibull parameters
  shape.deathL2=as.numeric(mod.coef[1,5])
  scale.deathL2=as.numeric(mod.coef[2,5])
  
  #scale for weibull regression
  scale.reg.deathL2<-scale.deathL2*(exp(-1*shape.deathL2*(beta.trt.deathL2*trt+XB.deathL2)))^(-1/shape.deathL2)
  scale.reg.deathL2.alltreat<-scale.deathL2*(exp(-1*shape.deathL2*(beta.trt.deathL2+XB.deathL2)))^(-1/shape.deathL2)
  scale.reg.deathL2.notreat<-scale.deathL2*(exp(-1*shape.deathL2*(XB.deathL2)))^(-1/shape.deathL2)
  
  ######Death from ED 
  #coefficients
  beta.trt.deathED<-as.numeric(mod.coef[3,6]) #Weibull coef are already reversed
  beta.coef.deathED<-mod.coef[c(4:7,9),6]
  #effect of covariates
  XB.deathED<-as.matrix(model.covar)%*%as.matrix(beta.coef.deathED)
  #weibull parameters
  shape.deathED=as.numeric(mod.coef[1,6])
  scale.deathED=as.numeric(mod.coef[2,6])
  
  #scale for weibull regression
  scale.reg.deathED<-scale.deathED*(exp(-1*shape.deathED*(beta.trt.deathED*trt+XB.deathED)))^(-1/shape.deathED)
  scale.reg.deathED.alltreat<-scale.deathED*(exp(-1*shape.deathED*(beta.trt.deathED+XB.deathED)))^(-1/shape.deathED)
  scale.reg.deathED.notreat<-scale.deathED*(exp(-1*shape.deathED*(XB.deathED)))^(-1/shape.deathED)
  
}



###### Log normal regressions - need to set up params differently ###################################
if(out.type=="lnormClayton"){
  mod.coef<-lnorm.coef
  
  
  ######Progression from L1 
  #coefficients
  beta.trt.progL1<-as.numeric(mod.coef[3,2]) 
  beta.coef.progL1<-mod.coef[c(4:7,9),2] # No need for ECOG NA
  #effect of covariates
  XB.progL1<-as.matrix(model.covar)%*%as.matrix(beta.coef.progL1)
  #parameters
  meanlog.progL1=as.numeric(mod.coef[1,2])
  sdlog.progL1=as.numeric(mod.coef[2,2])
  
  
  meanlog.reg.progL1<- (meanlog.progL1+beta.trt.progL1*trt+XB.progL1)
  meanlog.reg.progL1.alltreat<-(meanlog.progL1+beta.trt.progL1+XB.progL1)
  meanlog.reg.progL1.notreat<-(meanlog.progL1+XB.progL1)
  
  ######Death from L1 
  #coefficients
  beta.trt.deathL1<-as.numeric(mod.coef[3,3])
  beta.coef.deathL1<-mod.coef[c(4:7,9),3]
  #effect of covariates
  XB.deathL1<-as.matrix(model.covar)%*%as.matrix(beta.coef.deathL1)
  #lnorm params
  meanlog.deathL1=as.numeric(mod.coef[1,3])
  sdlog.deathL1=as.numeric(mod.coef[2,3])
  
  
  meanlog.reg.deathL1<-meanlog.deathL1+(beta.trt.deathL1*trt+XB.deathL1)
  meanlog.reg.deathL1.alltreat<-meanlog.deathL1+(beta.trt.deathL1+XB.deathL1)
  meanlog.reg.deathL1.notreat<-meanlog.deathL1+(XB.deathL1)
  
  ######Progression from L2 
  #coefficients
  beta.trt.progL2<-as.numeric(mod.coef[3,4]) 
  beta.coef.progL2<-mod.coef[c(4:7,9),4] 
  #effect of covariates
  XB.progL2<-as.matrix(model.covar)%*%as.matrix(beta.coef.progL2)
  #weibull parameters
  meanlog.progL2=as.numeric(mod.coef[1,4])
  sdlog.progL2=as.numeric(mod.coef[2,4])
  
  
  meanlog.reg.progL2<- (meanlog.progL2+beta.trt.progL2*trt+XB.progL2)
  meanlog.reg.progL2.alltreat<-(meanlog.progL2+beta.trt.progL2+XB.progL2)
  meanlog.reg.progL2.notreat<-(meanlog.progL2+XB.progL2)
  
  ######Death from L2 
  #coefficients
  beta.trt.deathL2<-as.numeric(mod.coef[3,3]) 
  beta.coef.deathL2<-mod.coef[c(4:7,9),3]
  #effect of covariates
  XB.deathL2<-as.matrix(model.covar)%*%as.matrix(beta.coef.deathL2)
  #lnorm params
  meanlog.deathL2=as.numeric(mod.coef[1,3])
  sdlog.deathL2=as.numeric(mod.coef[2,3])
  
  
  meanlog.reg.deathL2<-meanlog.deathL2+(beta.trt.deathL2*trt+XB.deathL2)
  meanlog.reg.deathL2.alltreat<-meanlog.deathL2+(beta.trt.deathL2+XB.deathL2)
  meanlog.reg.deathL2.notreat<-meanlog.deathL2+(XB.deathL2)
  
  
  ######Death from L2 
  #coefficients
  beta.trt.deathED<-as.numeric(mod.coef[3,6]) 
  beta.coef.deathED<-mod.coef[c(4:7,9),6]
  #effect of covariates
  XB.deathED<-as.matrix(model.covar)%*%as.matrix(beta.coef.deathED)
  #weibull parameters
  meanlog.deathED=as.numeric(mod.coef[1,3])
  sdlog.deathED=as.numeric(mod.coef[2,3])
  
  
  meanlog.reg.deathED<-meanlog.deathED+(beta.trt.deathED*trt+XB.deathED)
  meanlog.reg.deathED.alltreat<-meanlog.deathED+(beta.trt.deathED+XB.deathED)
  meanlog.reg.deathED.notreat<-meanlog.deathED+(XB.deathED)

}

###############################################################################
############ Generate outcomes - all "potential" transition times #############
###############################################################################

######  Generate outcomes independently
if(out.type=="ind"){
  
  #Generate outcomes independently
  set.seed(125)
  pfs1<-pfs2<-death1<-death2<-deathED<-rep(NA,N)
  for (i in 1:N) {
    pfs1[i]<-rweibull(1, shape=shape.progL1, scale=scale.reg.progL1[i])
    pfs2[i]<-rweibull(1, shape=shape.progL2, scale=scale.reg.progL2[i])
    death1[i]<-rweibull(1, shape=shape.deathL1, scale=scale.reg.deathL1[i])
    death2[i]<-rweibull(1, shape=shape.deathL2, scale=scale.reg.deathL2[i])
    deathED[i]<-rweibull(1, shape=shape.deathED, scale=scale.reg.deathED[i])
  }
  
  cor(cbind(pfs1,death1,pfs2,death2,deathED))
  
  outcome<-cbind(pfs1,death1,pfs2,death2,deathED)
  
  
}


seed=125
set.seed(seed)

#Generate correlated outcomes
n<-dim(covariates)[1]

#################### Set the parameters for the copula

## unstructured correlation
##     pfs1 os1 pfs2 os2 osed
#pfs1  1    0.2 0    0   .3
#os1        1   0    0   .3
#pfs2           1    .2  .3
#os2                1    .3
#osed                    1  

#Note: correlations from observed data:
#pfs1 pfs2 = 0
#pfs1 osed =0
#pfs2 osed=0.07
#Include a zero for pfs1 pfs2 but larger correlations for osed


if(out.type=="clayton" | out.type=="smallSS" | out.type=="lnormClayton" |out.type=="LogLogisClayton"){
  my.copula<-archmCopula(family="clayton",dim=5,param=0.5)
}else if(out.type=="GaussCS"){
  my.copula<-normalCopula(dim=5,dispstr="ex",param=0.2)
}else if(out.type=="GaussUN"){
  my.copula<-normalCopula(dim=5,dispstr="un",param=c(0.2,0,0,.3,0,0,.3,.2,.3,.3))
}else if(out.type=="TCS"){
  my.copula<-tCopula(dim=5,param=0.2)
}else if(out.type=="TUN"){
  my.copula<-tCopula(dim=5,dispstr="un",param=c(0.2,0,0,.3,0,0,.3,.2,.3,.3))
}

####### simulate outcomes via the copula ###########
if(out.type!="ind"  ){
  margins.input = c("weibull","weibull","weibull","weibull","weibull")
  if(out.type=="lnormClayton"){margins.input = c("lnorm","lnorm","lnorm","lnorm","lnorm")}
  if(out.type=="LogLogisClayton"){margins.input = c("llogis","llogis","llogis","llogis","llogis")}
  
  
  
  
  outcome<-matrix(rep(NA,5*n),nrow=n)
  for (i in 1:n) {
    
    if(out.type!="lnormClayton"){paramMargins.input = list(list(shape=shape.progL1,
                                                                scale=scale.reg.progL1[i]),
                                                           list(shape=shape.deathL1
                                                                ,scale=scale.reg.deathL1[i]),
                                                           list(shape=shape.progL2
                                                                ,scale=scale.reg.progL2[i]),
                                                           list(shape=shape.deathL2
                                                                ,scale=scale.reg.deathL2[i]),
                                                           list(shape=shape.deathED
                                                                ,scale=scale.reg.deathED[i]))}
    
    if(out.type=="lnormClayton"){paramMargins.input = list(list(meanlog=meanlog.reg.progL1[i],
                                                                sdlog=sdlog.progL1),
                                                           list(meanlog=meanlog.reg.deathL1[i]
                                                                ,sdlog=sdlog.deathL1),
                                                           list(meanlog=meanlog.reg.progL2[i]
                                                                ,sdlog=sdlog.progL2),
                                                           list(meanlog=meanlog.reg.deathL2[i]
                                                                ,sdlog=sdlog.deathL2),
                                                           list(meanlog=meanlog.reg.deathED[i]
                                                                ,sdlog=sdlog.deathED))}
    
    my.mvd<-mvdc(copula=my.copula, margins = margins.input ,
                 paramMargins = paramMargins.input
                 
    )
    outcome[i,]<-rMvdc(1,my.mvd)
  }
  
}

#Check correlations of outcomes
cor(outcome)


outcome<-data.frame(outcome
)
names(outcome)<-c("PFS1", "OS1", "PFS2", "OS2", "OSED")

############################ Additional outcomes ############
### simulate going directly from L1 to ED

#based on probability, simulate whether they skip L2
skip.L2<-rep(NA,N) 

for (i in 1:N) {
  skip.L2[i]<-rbinom(1,1,covariates$Pr.L2.0[i])
}


## Simulate censoring times

cens.T<-runif(N, min=6, max=72)

###############################################################################
############ Create synthetic dataset ("observed" transitions) ################
###############################################################################


simdf<-covariates %>%
  rename(PatientID=ID)

simdf$L1Cis<-trt
simdf$AnyL2Trt<-as.numeric(!skip.L2)
simdf$cens.T=cens.T

#for each patient, find path and transition times

simdf$PFS_months<-simdf$Progression<-simdf$PFS_monthsL2<-simdf$ProgressionL2<-simdf$Died<-simdf$OS_months_L1<-
  simdf$OS_months_L2<-rep(NA,N)

simdf$PFS_months<-simdf$Progression<-
  simdf$PFS_monthsL2<-simdf$ProgressionL2<-
  simdf$Died<-simdf$OS_months_L1<-
  simdf$OS_months_L2<-simdf$OS_months_ED<-rep(NA,N)


for (i in 1:N){
  
  modtime=0
  #Find the next stage
  simtimes=c(outcome$PFS1[i],outcome$OS1[i],cens.T[i])
  L1time=min(simtimes)
  next.stage=c("L2","D","NA")[which.min(simtimes)]
  if (next.stage=="L2" & simdf$AnyL2Trt[i]==0){next.stage="ED"}
  
  #Assign L1 outcomes
  simdf$PFS_months[i]<-L1time
  #Progression =1 if they died OR if they progressed
  simdf$Progression[i]=as.numeric(next.stage!="NA")
  simdf$Died[i]=as.numeric(next.stage=="D")
  
  modtime=modtime+L1time
  cens.T.state=cens.T[i]-modtime
  #L2 - if applicable - find next stage
  if(next.stage=="L2"){
    simtimes=c(outcome$PFS2[i],outcome$OS2[i],cens.T.state)
    L2time=min(simtimes)
    next.stage=c("ED","D","NA")[which.min(simtimes)]
    
    #Assign L2 outcomes
    simdf$PFS_monthsL2[i]<-L2time
    #Progression =1 if they died OR if they progressed
    simdf$ProgressionL2[i]=as.numeric(next.stage!="NA")
    if(next.stage=="D"){simdf$Died[i]=1}
    
    modtime=modtime+L2time
    cens.T.state=cens.T[i]-modtime   ### TO FIX - subtract L2 time instead fo mod time
  }
  
  if(next.stage=="ED"){
    simtimes=c(outcome$OSED[i],cens.T.state)
    EDtime=min(simtimes)
    next.stage=c("D","NA")[which.min(simtimes)]
    
    #Death from ED
    simdf$OS_months_ED[i]<-EDtime
    if(next.stage=="D"){simdf$Died[i]=1}
    
    modtime=modtime+EDtime
  }
  
  #Assign final OS times
  simdf$OS_months_L1[i]<-modtime
  if(!is.na(simdf$PFS_monthsL2[i])){simdf$OS_months_L2[i]=modtime-simdf$PFS_months[i]}
  
}

#Add text for directory here
path<-""


filename=paste(out.type,seed,".csv", sep="")


write.csv(simdf, paste(path, filename, sep=""))




###############################################################################
############ Now find the "truth" - how much time is spent in each state ######
############ potential outcomes of treatment, no censoring ####################
###############################################################################


############ All sequence 1 ################

seed=125
set.seed(seed)


#Generate outcomes independently - all sequence 1
if(out.type=="ind"){
  set.seed(125)
  pfs1<-pfs2<-death1<-death2<-deathED<-rep(NA,N)
  for (i in 1:N) {
    pfs1[i]<-rweibull(1, shape=shape.progL1, scale=scale.reg.progL1.alltreat[i])
    pfs2[i]<-rweibull(1, shape=shape.progL2, scale=scale.reg.progL2.alltreat[i])
    death1[i]<-rweibull(1, shape=shape.deathL1, scale=scale.reg.deathL1.alltreat[i])
    death2[i]<-rweibull(1, shape=shape.deathL2, scale=scale.reg.deathL2.alltreat[i])
    deathED[i]<-rweibull(1, shape=shape.deathED, scale=scale.reg.deathED.alltreat[i])
  }
  
  cor(cbind(pfs1,death1,pfs2,death2,deathED))
  
  outcome.alltreat<-cbind(pfs1,death1,pfs2,death2,deathED)
  
}

# Generate outcome via copula - all sequence 1
if(out.type!="ind"  ){
  margins.input = c("weibull","weibull","weibull","weibull","weibull")
  if(out.type=="lnormClayton"){margins.input = c("lnorm","lnorm","lnorm","lnorm","lnorm")}
  if(out.type=="LogLogisClayton"){margins.input = c("llogis","llogis","llogis","llogis","llogis")}
  
  outcome.alltreat<-matrix(rep(NA,5*n),nrow=n)
  for (i in 1:n) {
    
    if(out.type!="lnormClayton"){paramMargins.input = list(list(shape=shape.progL1,
                                                                scale=scale.reg.progL1.alltreat[i]),
                                                           list(shape=shape.deathL1
                                                                ,scale=scale.reg.deathL1.alltreat[i]),
                                                           list(shape=shape.progL2
                                                                ,scale=scale.reg.progL2.alltreat[i]),
                                                           list(shape=shape.deathL2
                                                                ,scale=scale.reg.deathL2.alltreat[i]),
                                                           list(shape=shape.deathED
                                                                ,scale=scale.reg.deathED.alltreat[i]))}
    
    if(out.type=="lnormClayton"){paramMargins.input = list(list(meanlog=meanlog.reg.progL1.alltreat[i],
                                                                sdlog=sdlog.progL1),
                                                           list(meanlog=meanlog.reg.deathL1.alltreat[i]
                                                                ,sdlog=sdlog.deathL1),
                                                           list(meanlog=meanlog.reg.progL2.alltreat[i]
                                                                ,sdlog=sdlog.progL2),
                                                           list(meanlog=meanlog.reg.deathL2.alltreat[i]
                                                                ,sdlog=sdlog.deathL2),
                                                           list(meanlog=meanlog.reg.deathED.alltreat[i]
                                                                ,sdlog=sdlog.deathED))}
    
    my.mvd<-mvdc(copula=my.copula, margins = margins.input ,
                 paramMargins = paramMargins.input
                 
    )
    outcome.alltreat[i,]<-rMvdc(1,my.mvd)
  }
  
}

#### Potential transition times

outcome.alltreat<-data.frame(outcome.alltreat
)
names(outcome.alltreat)<-c("PFS1", "OS1", "PFS2", "OS2", "OSED")


simdf.alltreat<-covariates %>%
  rename(PatientID=ID)

simdf.alltreat$L1Cis<-trt
simdf.alltreat$AnyL2Trt<-as.numeric(!skip.L2)
#All patients censored at 60 months
simdf.alltreat$cens.T=rep(60, length(simdf.alltreat$PatientID))

### for each patient, find path and observed transitions

simdf.alltreat$PFS_months<-simdf.alltreat$Progression<-simdf.alltreat$PFS_monthsL2<-simdf.alltreat$ProgressionL2<-simdf.alltreat$Died<-simdf.alltreat$OS_months_L1<-
  simdf.alltreat$OS_months_L2<-rep(NA,N)

simdf.alltreat$PFS_months<-simdf.alltreat$Progression<-
  simdf.alltreat$PFS_monthsL2<-simdf.alltreat$ProgressionL2<-
  simdf.alltreat$Died<-simdf.alltreat$OS_months_L1<-
  simdf.alltreat$OS_months_L2<-simdf.alltreat$OS_months_ED<-rep(NA,N)


for (i in 1:N){
  
  modtime=0
  #Find the next stage
  simtimes=c(outcome.alltreat$PFS1[i],outcome.alltreat$OS1[i],simdf.alltreat$cens.T[i])
  L1time=min(simtimes)
  next.stage=c("L2","D","NA")[which.min(simtimes)]
  if (next.stage=="L2" & simdf.alltreat$AnyL2Trt[i]==0){next.stage="ED"}
  
  #Assign L1 outcome.alltreats
  simdf.alltreat$PFS_months[i]<-L1time
  #Progression =1 if they died OR if they progressed
  simdf.alltreat$Progression[i]=as.numeric(next.stage!="NA")
  simdf.alltreat$Died[i]=as.numeric(next.stage=="D")
  
  modtime=modtime+L1time
  cens.T.state=simdf.alltreat$cens.T[i]-modtime
  #L2 - if applicable - find next stage
  if(next.stage=="L2"){
    simtimes=c(outcome.alltreat$PFS2[i],outcome.alltreat$OS2[i],cens.T.state)
    L2time=min(simtimes)
    next.stage=c("ED","D","NA")[which.min(simtimes)]
    
    #Assign L2 outcome.alltreats
    simdf.alltreat$PFS_monthsL2[i]<-L2time
    #Progression =1 if they died OR if they progressed
    simdf.alltreat$ProgressionL2[i]=as.numeric(next.stage!="NA")
    if(next.stage=="D"){simdf.alltreat$Died[i]=1}
    
    modtime=modtime+L2time
    cens.T.state= cens.T.state-L2time
  }
  
  if(next.stage=="ED"){
    simtimes=c(outcome.alltreat$OSED[i],cens.T.state)
    EDtime=min(simtimes)
    next.stage=c("D","NA")[which.min(simtimes)]
    
    #Death from ED
    simdf.alltreat$OS_months_ED[i]<-EDtime
    if(next.stage=="D"){simdf.alltreat$Died[i]=1}
    
    modtime=modtime+EDtime
  }
  
  #Assign final OS times
  simdf.alltreat$OS_months_L1[i]<-modtime
  if(!is.na(simdf.alltreat$PFS_monthsL2[i])){simdf.alltreat$OS_months_L2[i]=modtime-simdf.alltreat$PFS_months[i]}
  
}

########### Now calculate costs and effects

#time spent in each state
time.state.alltreat<-simdf.alltreat %>%
  dplyr:::select(PFS_months, PFS_monthsL2, OS_months_ED)

time.state.alltreat[is.na(time.state.alltreat)]<-0

#split up into "on treatment" and "off treatment" times
time.state.alltreat$L1ontrt<-time.state.alltreat$PFS_months
time.state.alltreat$L1ontrt[time.state.alltreat$L1ontrt>6]<-6
time.state.alltreat$L1aftertrt<-time.state.alltreat$PFS_months-6
time.state.alltreat$L1aftertrt[time.state.alltreat$L1aftertrt<0]<-0

time.state.alltreat$L2ontrt<-time.state.alltreat$PFS_monthsL2
time.state.alltreat$L2ontrt[time.state.alltreat$L2ontrt>24]<-24
time.state.alltreat$L2aftertrt<-time.state.alltreat$PFS_monthsL2-24
time.state.alltreat$L2aftertrt[time.state.alltreat$L2aftertrt<0]<-0


time.state.alltreat<-time.state.alltreat %>%
  dplyr:::select(L1ontrt, L1aftertrt, L2ontrt, L2aftertrt, OS_months_ED)

#Costs

cost.management.L1<-4916
cost.management.L2<-2310
cost.management.ED<-4940

cost.L1Cis<-223
cost.L1Carbo<-229
cost.L1gem<-449
cost.L2IO<-9299 #cost for pembro 


#Summarize costs
trt.cost.state<-c(cost.management.L1+cost.L1Cis+cost.L1gem,
                  cost.management.L1,
                  cost.management.L2+cost.L2IO,
                  cost.management.L2,
                  cost.management.ED)

cost.per.pat<-as.matrix(time.state.alltreat) %*% trt.cost.state

cost.mean.alltreat<-mean(cost.per.pat)

#utilities

u.L1.Cis<-0.67
u.L1.Carbo<-0.7
u.L1.postTrt<-0.72
u.L1AE<-0.6 
#Line 2
u.L2.IO<-0.6
#Ed
u.ED<-0.52
u.state<- c(u.L1.Cis, u.L1.postTrt, u.L2.IO, u.L2.IO, u.ED)

#summarize outcomes
qaly.per.pat<-as.matrix(time.state.alltreat) %*% u.state

qaly.mean.alltreat<-mean(qaly.per.pat)/12




############ All sequence 2 ################

seed=125
set.seed(seed)

#Generate outcomes independently - all sequence 2
if(out.type=="ind"){
  
  #Generate outcomes independently
  set.seed(125)
  pfs1<-pfs2<-death1<-death2<-deathED<-rep(NA,N)
  for (i in 1:N) {
    pfs1[i]<-rweibull(1, shape=shape.progL1, scale=scale.reg.progL1.notreat[i])
    pfs2[i]<-rweibull(1, shape=shape.progL2, scale=scale.reg.progL2.notreat[i])
    death1[i]<-rweibull(1, shape=shape.deathL1, scale=scale.reg.deathL1.notreat[i])
    death2[i]<-rweibull(1, shape=shape.deathL2, scale=scale.reg.deathL2.notreat[i])
    deathED[i]<-rweibull(1, shape=shape.deathED, scale=scale.reg.deathED.notreat[i])
  }
  
  cor(cbind(pfs1,death1,pfs2,death2,deathED))
  
  outcome.notreat<-cbind(pfs1,death1,pfs2,death2,deathED)
  
}

#Generate outcomes via copula - all sequence 2
if(out.type!="ind"  ){
  margins.input = c("weibull","weibull","weibull","weibull","weibull")
  if(out.type=="lnormClayton"){margins.input = c("lnorm","lnorm","lnorm","lnorm","lnorm")}
  if(out.type=="LogLogisClayton"){margins.input = c("llogis","llogis","llogis","llogis","llogis")}
  
  outcome.notreat<-matrix(rep(NA,5*n),nrow=n)
  for (i in 1:n) {
    
    if(out.type!="lnormClayton"){paramMargins.input = list(list(shape=shape.progL1,
                                                                scale=scale.reg.progL1.notreat[i]),
                                                           list(shape=shape.deathL1
                                                                ,scale=scale.reg.deathL1.notreat[i]),
                                                           list(shape=shape.progL2
                                                                ,scale=scale.reg.progL2.notreat[i]),
                                                           list(shape=shape.deathL2
                                                                ,scale=scale.reg.deathL2.notreat[i]),
                                                           list(shape=shape.deathED
                                                                ,scale=scale.reg.deathED.notreat[i]))}
    
    if(out.type=="lnormClayton"){paramMargins.input = list(list(meanlog=meanlog.reg.progL1.notreat[i],
                                                                sdlog=sdlog.progL1),
                                                           list(meanlog=meanlog.reg.deathL1.notreat[i]
                                                                ,sdlog=sdlog.deathL1),
                                                           list(meanlog=meanlog.reg.progL2.notreat[i]
                                                                ,sdlog=sdlog.progL2),
                                                           list(meanlog=meanlog.reg.deathL2.notreat[i]
                                                                ,sdlog=sdlog.deathL2),
                                                           list(meanlog=meanlog.reg.deathED.notreat[i]
                                                                ,sdlog=sdlog.deathED))}
    
    my.mvd<-mvdc(copula=my.copula, margins = margins.input ,
                 paramMargins = paramMargins.input
                 
    )
    outcome.notreat[i,]<-rMvdc(1,my.mvd)
  }
  
}

#### Potential transition times
outcome.notreat<-data.frame(outcome.notreat
)
names(outcome.notreat)<-c("PFS1", "OS1", "PFS2", "OS2", "OSED")


### for each patient, find path and observed transitions
simdf.notreat<-covariates %>%
  rename(PatientID=ID)

simdf.notreat$L1Cis<-trt
simdf.notreat$AnyL2Trt<-as.numeric(!skip.L2)
#All patients censored at 60 months
simdf.notreat$cens.T=rep(60, length(simdf.notreat$PatientID))

#for each patient, find path and transition times

simdf.notreat$PFS_months<-simdf.notreat$Progression<-simdf.notreat$PFS_monthsL2<-simdf.notreat$ProgressionL2<-simdf.notreat$Died<-simdf.notreat$OS_months_L1<-
  simdf.notreat$OS_months_L2<-rep(NA,N)

simdf.notreat$PFS_months<-simdf.notreat$Progression<-
  simdf.notreat$PFS_monthsL2<-simdf.notreat$ProgressionL2<-
  simdf.notreat$Died<-simdf.notreat$OS_months_L1<-
  simdf.notreat$OS_months_L2<-simdf.notreat$OS_months_ED<-rep(NA,N)


for (i in 1:N){
  
  modtime=0
  #Find the next stage
  simtimes=c(outcome.notreat$PFS1[i],outcome.notreat$OS1[i],simdf.notreat$cens.T[i])
  L1time=min(simtimes)
  next.stage=c("L2","D","NA")[which.min(simtimes)]
  if (next.stage=="L2" & simdf.notreat$AnyL2Trt[i]==0){next.stage="ED"}
  
  #Assign L1 outcome.notreats
  simdf.notreat$PFS_months[i]<-L1time
  #Progression =1 if they died OR if they progressed
  simdf.notreat$Progression[i]=as.numeric(next.stage!="NA")
  simdf.notreat$Died[i]=as.numeric(next.stage=="D")
  
  modtime=modtime+L1time
  cens.T.state=simdf.notreat$cens.T[i]-modtime
  #L2 - if applicable - find next stage
  if(next.stage=="L2"){
    simtimes=c(outcome.notreat$PFS2[i],outcome.notreat$OS2[i],cens.T.state)
    L2time=min(simtimes)
    next.stage=c("ED","D","NA")[which.min(simtimes)]
    
    #Assign L2 outcome.notreats
    simdf.notreat$PFS_monthsL2[i]<-L2time
    #Progression =1 if they died OR if they progressed
    simdf.notreat$ProgressionL2[i]=as.numeric(next.stage!="NA")
    if(next.stage=="D"){simdf.notreat$Died[i]=1}
    
    modtime=modtime+L2time
    cens.T.state= cens.T.state-L2time
  }
  
  if(next.stage=="ED"){
    simtimes=c(outcome.notreat$OSED[i],cens.T.state)
    EDtime=min(simtimes)
    next.stage=c("D","NA")[which.min(simtimes)]
    
    #Death from ED
    simdf.notreat$OS_months_ED[i]<-EDtime
    if(next.stage=="D"){simdf.notreat$Died[i]=1}
    
    modtime=modtime+EDtime
  }
  
  #Assign final OS times
  simdf.notreat$OS_months_L1[i]<-modtime
  if(!is.na(simdf.notreat$PFS_monthsL2[i])){simdf.notreat$OS_months_L2[i]=modtime-simdf.notreat$PFS_months[i]}
  
}

#Now calculate costs and effects

#time spent in each state
time.state.notreat<-simdf.notreat %>%
  dplyr:::select(PFS_months, PFS_monthsL2, OS_months_ED)

time.state.notreat[is.na(time.state.notreat)]<-0

#split up into on treatment and off treatment times
time.state.notreat$L1ontrt<-time.state.notreat$PFS_months 
time.state.notreat$L1ontrt[time.state.notreat$L1ontrt>6]<-6
time.state.notreat$L1aftertrt<-time.state.notreat$PFS_months-6
time.state.notreat$L1aftertrt[time.state.notreat$L1aftertrt<0]<-0

time.state.notreat$L2ontrt<-time.state.notreat$PFS_monthsL2
time.state.notreat$L2ontrt[time.state.notreat$L2ontrt>24]<-24
time.state.notreat$L2aftertrt<-time.state.notreat$PFS_monthsL2-24
time.state.notreat$L2aftertrt[time.state.notreat$L2aftertrt<0]<-0


time.state.notreat<-time.state.notreat %>%
  dplyr:::select(L1ontrt, L1aftertrt, L2ontrt, L2aftertrt, OS_months_ED)


#Treatment costs

cost.management.L1<-4916
cost.management.L2<-2310
cost.management.ED<-4940

cost.L1Cis<-223
cost.L1Carbo<-229
cost.L1gem<-449
cost.L2IO<-9299 #cost for pembro 


trt.cost.state<-c(cost.management.L1+cost.L1Carbo+cost.L1gem,
                  cost.management.L1,
                  cost.management.L2+cost.L2IO,
                  cost.management.L2,
                  cost.management.ED)

#utilities

u.L1.Cis<-0.67
u.L1.Carbo<-0.7
u.L1.postTrt<-0.72
u.L1AE<-0.6 
#Line 2
u.L2.IO<-0.6
#Ed
u.ED<-0.52
u.state<- c(u.L1.Cis, u.L1.postTrt, u.L2.IO, u.L2.IO, u.ED)


### summarize costs and outcomes

cost.per.pat<-as.matrix(time.state.notreat) %*% trt.cost.state

cost.mean.notreat<-mean(cost.per.pat)
qaly.per.pat<-as.matrix(time.state.notreat) %*% u.state

qaly.mean.notreat<-mean(qaly.per.pat)/12


# "True" outcomes - costs and effects for sequence 1 and sequence 2


print("sequence 1 cost and QALY")
print(c(cost.mean.alltreat, qaly.mean.alltreat
))
print("sequence 2 cost and QALY")
print(c(cost.mean.notreat, qaly.mean.notreat
))

