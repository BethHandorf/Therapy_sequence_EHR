####################################################################
# This code pre-fits a multi-state model to the synthetic dataset ##
# It also includes a function to extract the  model-based ##########
# transition matrix for a specific patient conditional on model time
# Results are saved in a .RData file for later use #################
# This code uses the R package mstate to fit the models and ########
# obtain relevant predictions ######################################
####################################################################

library(MASS)
library(dplyr)
library(corrplot)
library(fGarch)
library(twang)
library(nnet)
library(riskRegression)

library(prodlim)
library(lava)
library(survival)
library(flexsurv)
library(mstate)

#Clear out prior data
rm(list=ls())




########################### Import and set up data ##########################################
#Read in synthetic data

#out.type<-"ind"
out.type<-"clayton"
#out.type<-"GaussCS"
#out.type<-"GaussUN"
#out.type<-"TCS"
#out.type<-"TUN"
#out.type<-"smallSS"
#out.type<-"LogLogisClayton"
#out.type<-"lnormClayton"

set.seed(50)

df.orig<-read.csv(paste(out.type,"125.csv", sep=""))


#Set up covariates.
#Filter out misingness (note: synthetic data has no missingness)
df.orig$MDRD_GFR[df.orig$MDRD_GFR>200]<-NA
df.orig$hemoglobin[df.orig$hemoglobin>200]<-NA

df.orig.covar<- df.orig %>%
  dplyr::select(Surgery, Male, RaceB, RaceO, Hispanic, AgeAdvDx, EcogValue,  
          MDRD_GFR, hemoglobin, albumin, CHF_known,        
           hear_known, neuropathy_known, viscmet_known)

cor(df.orig$PFS_months, df.orig$PFS_monthsL2, use="pairwise.complete.obs")

countNA<-function(x) {sum(is.na(x))}

apply(df.orig, 2, countNA)

trt1df<-df.orig
trt1df<- trt1df %>%
  filter(!is.na(MDRD_GFR)) %>%
  filter(!is.na(hemoglobin)) %>%
  filter(!is.na(albumin)) %>%
  filter(!is.na(CHF_known)) %>%
  filter(!is.na(OS_months_L1))

trt1df$EcogValue<-as.factor(trt1df$EcogValue)
trt1df$EcogValue<-addNA(trt1df$EcogValue)



df.msm<-df.orig %>%
  mutate(EcogValue=factor(EcogValue))
  

####### Set up mstate model ####################################


#transition matrix
#states: L1, L2, ED, D (1,2,3,4)
#From (rows) to (cols) - 6 total transitions possible

#   L1 L2 ED D
#L1 -  1  2  3
#L2 -  -  4  5
#ED -  -  -  6
#D  -  -  -  -

tmat<- transMat( x = list(c(2,3,4),c(3,4),c(4),c()),names=c("L1","L2","ED","D"))

tmat

#All possible paths through model
paths(tmat)

#Need an indicator for being in the ED cohort:
#Progressed on line 2
tmp<-!is.na(df.msm$PFS_monthsL2)
tmp[(df.msm$PFS_monthsL2+.5 >= df.msm$OS_months_L2)]<-0
tmp[is.na(tmp)]<-0
#OR came directly from line1
tmp2<-as.numeric(df.msm$Progression==1 & df.msm$AnyL2Trt==0)
tmp2[(df.msm$PFS_months+.5 >= df.msm$OS_months_L1)]<-0
df.msm$InEDCohort<-as.numeric(tmp==1 | tmp2==1)

#time in ED state
df.msm$TimeEDStart<-rowSums(cbind(df.msm$PFS_months,df.msm$PFS_monthsL2),na.rm=TRUE)
df.msm$TimeEDStart[df.msm$InEDCohort==0]<-df.msm$OS_months_L1[df.msm$InEDCohort==0]

#reformat the data so that it is set up the way mstate wants data set up.
#L2 = time entering into L2 or last follow-up otherwise
#L2.s = status variable indicating event of entering L2
tmpdf<- df.msm %>%
  rename(L2.s=AnyL2Trt) %>%
  mutate(L2.s = ifelse(is.na(L2.s), 0, L2.s)) %>%
  mutate(L2=PFS_months*L2.s+OS_months_L1*(1-L2.s)) %>%
  rename(ED.s=InEDCohort) %>%
  rename(ED=TimeEDStart) %>%
  rename(D.s=Died) %>%
  rename(D = OS_months_L1) %>%
    dplyr:::select(PatientID,L2, L2.s,ED, ED.s, D, D.s,L1Cis,Surgery,AgeAdvDx,       
     EcogValue,albumin)

#mstate function to transform data to long format
msdf<-msprep(data=tmpdf, trans=tmat, time=c(NA, "L2","ED","D"), status=c(NA,"L2.s","ED.s","D.s"),keep=c("PatientID","L1Cis","Surgery","AgeAdvDx","EcogValue","albumin"))

events(msdf)

covs<-c("L1Cis","Surgery","AgeAdvDx","EcogValue","albumin")

msdf<-expand.covs(msdf,covs, longnames=FALSE)

####### Run mstate model ####################################

#Model with covariates
cfull<-coxph(Surv(Tstart,Tstop,status)~L1Cis.1+L1Cis.2+L1Cis.3+L1Cis.4+L1Cis.5+L1Cis.6+Surgery.1+Surgery.2+Surgery.3+   
Surgery.4+Surgery.5+Surgery.6+AgeAdvDx.1+AgeAdvDx.2+  
AgeAdvDx.3+AgeAdvDx.4+AgeAdvDx.5+AgeAdvDx.6+EcogValue1.1+
EcogValue1.2+EcogValue1.3+EcogValue1.4+EcogValue1.5+EcogValue1.6+
EcogValue2.1+EcogValue2.2+EcogValue2.3+EcogValue2.4+EcogValue2.5+
EcogValue2.6+albumin.1+albumin.2+albumin.3+albumin.4+   
albumin.5+albumin.6+strata(trans), data=msdf, method="breslow")


####### Create fitted object for each subject ####################################


#Create a list with a msfit object for each subject   
#Predictions with sequence 1
pt.msfit.list.cis<-vector(mode="list", length=dim(df.msm)[1])
for (i in 1:dim(df.msm)[1]) {
  whA<-which(msdf$PatientID==df.msm$PatientID[i])
  patA<-rep(whA[1],6,)
  patA <- msdf[rep(whA[1], 6), 10:14]
  #Set all treatment to Yes
  patA$L1Cis<-rep(1,6)
  #New variable Trans
  patA$trans <- 1:6
  #reference the transition matrix
  attr(patA, "trans") <- tmat
  #Make dummy covariates
  patA <- expand.covs(patA, covs, longnames = FALSE)
  patA$strata <- patA$trans
  msfA <- msfit(cfull, patA, trans = tmat)
  
  pt.msfit.list.cis[[i]]<-msfA
  
}

#Predictions with sequence 2
pt.msfit.list.carbo<-vector(mode="list", length=dim(df.msm)[1])
for (i in 1:dim(df.msm)[1]) {
  whA<-which(msdf$PatientID==df.msm$PatientID[i])
  patA<-rep(whA[1],6,)
  patA <- msdf[rep(whA[1], 6), 10:14]
  #Set all treatment to Yes
  patA$L1Cis<-rep(0,6)
  #New variable Trans
  patA$trans <- 1:6
  #reference the transition matrix
  attr(patA, "trans") <- tmat
  #Make dummy covariates
  patA <- expand.covs(patA, covs, longnames = FALSE)
  patA$strata <- patA$trans
  msfA <- msfit(cfull, patA, trans = tmat)
  
  pt.msfit.list.carbo[[i]]<-msfA
  
}


####### Functions to extract relevant probabilities from patient-specific msfit objects ####################

# Function to extract the probability of transitioning from time t to t+1 given starting state
Prob.pat.t.state<-function(patientid, time, state, Cis) {
  

  #extract the msfit object from the list 
  if(Cis==1){msf.pat<-pt.msfit.list.cis[[which(df.msm$PatientID==patientid)]]
  }else if (Cis==0){msf.pat<-pt.msfit.list.carbo[[which(df.msm$PatientID==patientid)]]}
  

  start.time <- Sys.time()  
  #Get the predicted probabilities of transitioning
  pt.T <- probtrans(msf.pat, predt = time, variance=FALSE)
    end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  s.obj<-summary(pt.T, from=state)[[state]]

  #limit to time<t+1
  s.obj.tmp<-s.obj[s.obj$time<time+1,]
  
  #probabilities at the maximum time<t+1
  trans.probs.pt.t.s<-tail(s.obj.tmp,n=1)
  
  return(unlist(trans.probs.pt.t.s[2:5]))
  
}

start.time <- Sys.time()
#Test function
probs.tmp<-suppressWarnings(Prob.pat.t.state(patientid = 1, time=4, state=3, Cis=1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


#Function: given a starting time and state, calculate the conditional transition matrix
#with transitions for each discrete cycle

Prob.mat.pat<-function(patientid, time, state, Cis) {
  

  #extract the msfit object from the list 
  if(Cis==1){msf.pat<-pt.msfit.list.cis[[which(df.msm$PatientID==patientid)]]
  }else if (Cis==0){msf.pat<-pt.msfit.list.carbo[[which(df.msm$PatientID==patientid)]]}
  


  #Get the predicted probabilities of transitioning
  pt.T <- probtrans(msf.pat, predt = time, variance=FALSE)

  s.obj<-summary(pt.T, from=state)[[state]]


  #limit to time<t+1
  s.obj.tmp<-s.obj[s.obj$time<time+1,]
  
  #probabilities at the maximum time<t+1
  trans.probs.pt.t.s<-tail(s.obj.tmp,n=1)
  
  #repeat for each cycle until 61

  for (t.tmp in c((time+1):61)) {
    #limit to time<t+1
    s.obj.tmp<-s.obj[s.obj$time<t.tmp+1,]
  
    #probabilities at the maximum time<t+1
    trans.probs.pt.t.s<-rbind(trans.probs.pt.t.s,tail(s.obj.tmp,n=1))
  
  }
  
  surv.mat<-cbind(p.times=c(time:61), trans.probs.pt.t.s)

  return(surv.mat)
}




start.time <- Sys.time()
#test function
probs.tmp<-suppressWarnings(Prob.mat.pat(patientid = 1, time=4, state=3, Cis=1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



save.image(paste("msmFits_",out.type,"125.RData", sep=""))




