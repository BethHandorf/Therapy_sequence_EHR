
########################################################################################
# This code runs a microsimulation model based on observed patient trajectories ########
# After censoring occurs, transition probabilities come from a pre-fit #################
# multi-state model. mstate_setup_for_microsimulation.R MUST be run FIRST and the ######
# resulting .Rmd file is read in.  This code runs the  microsimulation, checks #########
# model calibration, and summarizes results.  It also re-fits the model, resampling ####
# the patient population and changing the random seed.  This should ideally be #########
# modified to run in parallel on a high-performance computing environment.   ##########
# Note there are placeholders in this code for AE states, but they are not used here. ## 
########################################################################################

library(MASS)
library(dplyr)

library(survival)
library(mstate)
library(questionr)
library(flexsurv)
library(ggplot2)
library(reshape2)



rm(list = ls())  # remove any variables in R's memory 

#Read in data and fit models

out.type<-"ind"
#out.type<-"clayton"
#out.type<-"GaussCS"
#out.type<-"GaussUN"
#out.type<-"TCS"
#out.type<-"TUN"
#out.type<-"smallSS"
#out.type<-"LogLogisClayton"
#out.type<-"lnormClayton"




load(paste("msmFits_",out.type,"125.RData",sep=""))

############################################################################################
################# mPCA microsimulation model ###################
############################################################################################

# Microsimulation model for metastatic prostate cancer

# Based on code from:


# 'Microsimulation modeling for health decision sciences using R: a tutorial' 
# Authors: Eline Krijkamp, Fernando Alarid-Escudero, 
#          Eva Enns, Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou
#
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial
#





##################################### Model input #########################################
# Model input
n.i.full   <- dim(df.msm)[1]                # number of simulated individuals
n.t   <-60                    # time horizon, 60 1 month cycles (5 years)

#6 states:
#LN1: line 1 therapy
#AEL1: No progression, Adverse event (line 1) 

#LN2: Line 2 (IO)
#AEL2: No progression, Adverse event (line 2) 

#ED: Extensive Disease
#D: Death

v.n   <- c("LN1","AEL1","LN2","AEL2", "ED", "Death")  
n.s   <- length(v.n)           # the number of health states
#v.M_1 <- rep("LN1", n.i)         # everyone begins in the healthy state 
d.c   <- d.e <- 0.03/12           # equal discounting of costs and QALYs 



#2 strategies:
# 1. Cis: Cis/Gem -> IO single agent
# 2. Carbo:  Carbo/Gem -> IO single agent

v.Trt <- c("Cis", "Carbo") # store the strategy names

############################ Define each patient's trajectory ##########################

#Transitions from L2

df.msm$transL1<-df.msm$Progression
df.msm$transL1[df.msm$Progression==1 & df.msm$AnyL2Trt==0]<-2
#If time in L2 is <1 month, skip and go to ED
df.msm$transL1[df.msm$Progression==1 & df.msm$AnyL2Trt==1 & df.msm$PFS_monthsL2<1]<-2
df.msm$transL1[df.msm$Progression==1 & (df.msm$PFS_months+.5 >= df.msm$OS_months_L1)]<-3 #Allow 2 week window due to deaths at month-level resolution
#Avoid collisions of deaths and progression events
df.msm$transL1[df.msm$Progression==1 & (ceiling(df.msm$PFS_months)== ceiling(df.msm$OS_months_L1))]<-3

#Transitions form L2

df.msm$transL2<-df.msm$ProgressionL2
df.msm$transL2[df.msm$ProgressionL2==1 & (df.msm$PFS_monthsL2+.5 >= df.msm$OS_months_L2)]<-2 #Allow 2 week window due to deaths at month-level resolution
#Avoid collisions of deaths and progression events
df.msm$transL2[df.msm$ProgressionL2==1 & (ceiling(df.msm$PFS_monthsL2+df.msm$PFS_months)== ceiling(df.msm$OS_months_L1))]<-2

#matrices for each patient/time, one matrix for every possible transition

tj.L1.L2.Carbo<-
  tj.L1.ED.Carbo<-
  tj.L1.Death.Carbo<-
  tj.L1.L2.Cis<-
  tj.L1.ED.Cis<-
  tj.L1.Death.Cis<-
  tj.L2.ED<-
  tj.L2.Death<-
  tj.ED.Death<-
  matrix(rep(0,60*dim(df.msm)[1]),ncol=60)


#Now fill in the transition matrices with the appropriate indicators

for (i in 1:dim(df.msm)[1]){
  #Code a 1 when a transition occurs
  #Transitions from line 1
  if(df.msm$transL1[i]==1 & df.msm$PFS_months[i]<60 ){
    tj.L1.L2.Carbo[i,ceiling(df.msm$PFS_months[i])]<-1
    tj.L1.L2.Cis[i,ceiling(df.msm$PFS_months[i])]<-1
  }
  
  if(df.msm$transL1[i]==2 & df.msm$PFS_months[i]<60 ){
    tj.L1.ED.Carbo[i,ceiling(df.msm$PFS_months[i])]<-1
    tj.L1.ED.Cis[i,ceiling(df.msm$PFS_months[i])]<-1
  }
  
  if(df.msm$transL1[i]==3 & df.msm$PFS_months[i]<60 ){
    tj.L1.Death.Carbo[i,ceiling(df.msm$PFS_months[i])]<-1
    tj.L1.Death.Cis[i,ceiling(df.msm$PFS_months[i])]<-1
  }
  #Transitions from line 2
  #They had a transition at all
  if(!is.na(df.msm$transL2[i])&!is.na(df.msm$OS_months_L2[i])){
    if( df.msm$OS_months_L2[i]>0.5){
      #model time for pfs line 2
      pfs2time<-df.msm$PFS_months[i]+df.msm$PFS_monthsL2[i]
      if(df.msm$transL2[i]==1 & pfs2time<60 ){
        tj.L2.ED[i,ceiling(pfs2time)]<-1
      }
      if(df.msm$transL2[i]==2 & pfs2time<60 ){
        tj.L2.Death[i,ceiling(pfs2time)]<-1
      }
    }
  }
  #transitions from ED state
  #first, did they get to the ED state at all
  if( ((!is.na(df.msm$transL2[i])&df.msm$transL2[i]==1) | df.msm$transL1[i]==2 )  &
      !is.na(df.msm$OS_months_L1[i]) & 
      df.msm$transL1[i]!=3) {
    if(df.msm$Died[i]==1 & df.msm$OS_months_L1[i]<60 ){
      tj.ED.Death[i,ceiling(df.msm$OS_months_L1[i])]<-1
      
    }
  }
  
}



#Now, for each patient, find their censoring time
cens<-rep(0,dim(df.msm)[1])
cens.time<-rep(1,dim(df.msm)[1])

for (i in 1:dim(df.msm)[1]){
  #patients who never left L1
  if(df.msm$transL1[i]==0){
    cens[i]<-1
    cens.time[i]<-floor(df.msm$PFS_months[i])
    #Next, pts who never left L2
  }else if(!is.na(df.msm$transL2[i]) & df.msm$transL2[i]==0 ){
    cens[i]<-1
    cens.time[i]<-floor(df.msm$PFS_months[i]+df.msm$PFS_monthsL2[i])
  }else if(((!is.na(df.msm$transL2[i])&df.msm$transL2[i]==1) | df.msm$transL1[i]==2 )& !is.na(df.msm$OS_months_L1[i]) & df.msm$Died[i]==0) {
    cens[i]<-1
    cens.time[i]<-floor(df.msm$OS_months_L1[i])
  }
}

#deaths

tj.death<-rowSums(tj.L1.Death.Carbo+tj.L2.Death+tj.ED.Death)


sum(cens==0  & df.msm$OS_months_L1<60)
sum(tj.death[cens==0  & df.msm$OS_months_L1<60])

tmp<-df.msm[(tj.death==0 & cens==0  & df.msm$OS_months_L1<60),]



#######
# Other transition probabilities (per cycle)
#######

#AE transition probabilities would go here 


########
#Costs
########
#per month costs

cost.management.L1<-4916
cost.management.L2<-2310
cost.management.ED<-4940

cost.L1Cis<-223
cost.L1Carbo<-229
cost.L1gem<-449
cost.L2IO<-9299 #cost for pembro 



########
#Utilities
########
#Line 1
u.L1.Cis<-0.67
u.L1.Carbo<-0.7
u.L1.postTrt<-0.72

u.L1AE<-0.6 

#Line 2
u.L2.IO<-0.6

#Ed
u.ED<-0.52



##################################### Functions ###########################################

MicroSim <- function(v.M_1, n.i=n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE,TR.out.wt=TRUE, Trt = 5, seed = 1, samIDs, wt=rep(1,length(samIDs))) {

  # Arguments:  
  # v.M_1:   vector of initial states for individuals 
  # n.i:     number of individuals
  # n.t:     total number of cycles to run the model
  # v.n:     vector of health state names
  # d.c:     discount rate for costs
  # d.e:     discount rate for health outcome (QALYs)
  # TR.out:  should the output include a Microsimulation trace? (default is TRUE)
  # TS.out:  should the output include a transition array between states? (default is TRUE)
  # Trt:     are the n.i individuals receiving treatment? (scalar with a Boolean value, default is FALSE)
  # seed:    starting seed number for random number generator (default is 1)
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  
  v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e
  
  # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point
  # Added matricies to keep track of time spent in the current state, and time spent in line2
  m.M <- m.C <- m.E <- m.stime<- m.l1time<-m.l2time<-m.pdtime<- matrix(nrow = n.i, ncol = n.t + 1, 
                                                                       dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                                                                       paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state   
  m.stime[, 1]<-1    #Initial time in cycle is 1
  m.l1time[, 1]<-0    #Initial time in line 1 is 0
  m.l2time[, 1]<-0    #Initial time in line 2 is 0
  m.pdtime[, 1]<-0    #Initial time in PD is 0
  
  v.L1.tot.time<-rep(0,n.i)
  
  for (i in 1:n.i) {
    
    ptransmat=NULL
    #print(paste("i=",i))
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m.C[i, 1] <- Costs(m.M[i, 1], Trt)  # estimate costs per individual of the initial health state conditional on treatment
    m.E[i, 1] <- Effs (m.M[i, 1], Trt)  # estimate QALYs per individual of the initial health state conditional on treatment
    
    for (t in 1:n.t) {
      #print(paste("time=",t,"ID=",samIDs[i]))
      p.tmp <- Probs(m.M[i, t],stime=m.stime[i,t],mtime=t, l2time=m.l2time[i,t], pdtime=m.pdtime[i,t], Trt,L1.tot.time=v.L1.tot.time[i], ID=samIDs[i], ptransmat=ptransmat)           # calculate the transition probabilities at cycle t 
      
      v.p<-p.tmp[[1]]
      ptransmat<-p.tmp[[2]]
      # NOTE ASSUMPTION OF EXPONENTIAL HAZARDS HERE - times from first cycle ONLY.  
      # This only applies after censoring occurs.
      
      m.M[i, t + 1] <- sample(v.n, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      m.C[i, t + 1] <- Costs(m.M[i, t + 1], l1time=m.l1time[i,t], l2time=m.l2time[i,t], Trt)   # estimate costs per individual during cycle t + 1 conditional on treatment
      m.E[i, t + 1] <- Effs( m.M[i, t + 1], l1time=m.l1time[i,t], l2time=m.l2time[i,t], Trt)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
      
      #Update state time, line 2 time
      m.stime[i, t+1] <- as.numeric(m.M[i,t]==m.M[i, t + 1])* (m.stime[i, t]+1) + #Add 1 cycle if cycle doesn't change
        as.numeric(m.M[i,t]!=m.M[i, t + 1])*1                    #reset if new cycle
      
      #Time spent in L1
      v.L1.tot.time[i]<-v.L1.tot.time[i]+as.numeric(m.M[i, t + 1] %in% c("LN1","AEL1"))
      
      #increment L1 counter for every cycle in LN2, AEL2, PostL2
      m.l1time[i, t+1] <- m.l1time[i, t] + as.numeric(m.M[i, t + 1] %in% c("LN1","AEL1")) 
      
      #increment L2 counter for every cycle in LN2, AEL2, PostL2
      m.l2time[i, t+1] <- m.l2time[i, t] + as.numeric(m.M[i, t + 1] %in% c("LN2","AEL2")) 
      
      #increment PD counter for every cycle in PD
      m.pdtime[i, t+1] <- m.pdtime[i, t] + as.numeric(m.M[i, t + 1] %in% c("PD")) 
      
      
    } # close the loop for the time points 
    if(i/100 == round(i/100,0)) {          # display the progress of the simulation
      cat('\r', paste(i/n.i * 100, "% done", sep = " "))
    }
  } # close the loop for the individuals 
  
  tc <- m.C %*% v.dwc       # total discounted cost per individual
  te <- m.E %*% v.dwe       # total discounted QALYs per individual 
  
  tc_hat <- mean(tc)        # average discounted cost 
  te_hat <- mean(te)        # average discounted QALYs
  
  tc_hat_wt=mean(tc*wt)/mean(wt)  #Weighted averaged discounted cost
  te_hat_wt=mean(te*wt)/mean(wt)  #weighted averaged discounted QALYs
  
  if (TS.out == TRUE) {  # create a  matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns 
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) { # create a trace from the individual trajectories
    TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
    TR <- TR / n.i                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
    colnames(TR) <- v.n                                  # name the columns 
  } else {
    TR <- NULL
  }
  
  
  #Weighted trace function
  if (TR.out.wt == TRUE) { # create a trace from the individual trajectories
    TR.wt <- t(apply(m.M, 2, function(x) wtd.table(factor(x, levels = v.n, ordered = TRUE), weights=wt)))
    TR.wt <- TR.wt / sum(wt)                                       # create a distribution trace
    rownames(TR.wt) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
    colnames(TR.wt) <- v.n                                  # name the columns 
  } else {
    TR.wt <- NULL
  }
  
  
  results <- list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR, tc_hat_wt=tc_hat_wt, te_hat_wt=te_hat_wt, TR_wt=TR.wt) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  




#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.


Probs <- function(M_it, stime=1, mtime=1, l2time=0, pdtime=0, Trt=5, L1.tot.time=0, ID=1, ptransmat=NULL) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  # stime: state time
  # mtime: model time
  # l2time: line2 time
  # pdtime: time in PD
  # Trt: treatment code
  # ID: Patient ID (row number)
  
  v.p.it <- rep(0, n.s)     # create vector of state transition probabilities - initialize to zero
  names(v.p.it) <- v.n       # name the vector
  
  #v.p.it is a vector of length 8, each representing the probability of
  #transitioning to (or remaining in) a given state 
  #Elements of v.p.it are
  #v.p.it[6] LN1 (line 1)
  #v.p.it[5] AEL1 (adverse event from line 1 treatment)
  #v.p.it[4] LN2 (line 2)
  #v.p.it[3] AEL2 (adverse event from line 2 treatment)
  #v.p.it[2] ED ("extensive disease" state)
  #v.p.it[1] Death
  
  

  ################# while observed, using observed transitions #########
  if(cens[ID]==0 | mtime<cens.time[ID]) {
    if(M_it=="LN1") {
      #Probability of death
      v.p.it[6]<-(tj.L1.Death.Cis[ID,mtime] * as.numeric(Trt==1) + 
                    tj.L1.Death.Carbo[ID,mtime] * as.numeric(Trt==2) )
      #ED state
      v.p.it[5]<-
        (tj.L1.ED.Cis[ID,mtime] * as.numeric(Trt==1) + 
           tj.L1.ED.Carbo[ID,mtime] * as.numeric(Trt==2) )
      #L2
      v.p.it[3]<-
        (tj.L1.L2.Cis[ID,mtime] * as.numeric(Trt==1) + 
           tj.L1.L2.Carbo[ID,mtime] * as.numeric(Trt==2))
      
      #LN1
      v.p.it[1]<-1-sum(v.p.it)
    }
    
    
    if(M_it=="AEL1") {
      #Probability of death
      v.p.it[6]<-(tj.L1.Death.Cis[ID,mtime] * as.numeric(Trt==1) + 
                    tj.L1.Death.Carbo[ID,mtime] * as.numeric(Trt==2) )
      #ED state
      v.p.it[5]<-
        (tj.L1.ED.Cis[ID,mtime] * as.numeric(Trt==1) + 
           tj.L1.ED.Carbo[ID,mtime] * as.numeric(Trt==2) )
      #L2
      v.p.it[3]<-
        (tj.L1.L2.Cis[ID,mtime] * as.numeric(Trt==1) + 
           tj.L1.L2.Carbo[ID,mtime] * as.numeric(Trt==2))
      
      #AE

      #LN1
      v.p.it[1]<-1-sum(v.p.it)
    }
    
    if(M_it=="LN2") {
      #Probability of death
      v.p.it[6]<-(tj.L2.Death[ID,mtime]  )
      #ED state
      v.p.it[5]<-(tj.L2.ED[ID,mtime] )
      v.p.it[3]<-1-sum(v.p.it)
    }
    
    
    if(M_it=="AEL2") {
      #Probability of death
      v.p.it[6]<-(tj.L2.Death[ID,mtime]  )
      #ED state
      v.p.it[5]<-(tj.L2.ED[ID,mtime] )
      v.p.it[3]<-1-sum(v.p.it)
      
      
    }
    
    
    
    if(M_it=="ED") {
      #Probability of death
      #Probability of death
      v.p.it[6]<-(tj.ED.Death[ID,mtime]  )
      #ED state
      v.p.it[5]<-(1-sum(v.p.it) )
      
      
    }
    
    if(M_it=="Death") {
      #Probability of death
      v.p.it[6]<-1
      
    }
    
  }  
  
  ###### After censoring #######
  
  #If they are censored, use model-based probabilities
  if(cens[ID]==1 & mtime>=cens.time[ID]) {
    
    if(M_it=="LN1"){st=1}
    if(M_it=="AEL1"){st=1}
    if(M_it=="LN2"){st=2}
    if(M_it=="AEL2"){st=2}
    if(M_it=="ED"){st=3}
    if(M_it=="Death"){st=4}
    
    
    Cistrtcode=as.numeric(Trt==1)
    #call function to get probs based on current state, mtime, and patientID
    #Only need to do this for non-death states
    #if a transition occurred in the last cycle or ptransmat is not  yet defined for the patient
    if(st<4 & (stime==1 | is.null(ptransmat))){
      ptransmat_tmp<-suppressWarnings(Prob.mat.pat(patientid = ID, time=mtime,   state=st, Cis=Cistrtcode))
      ptransmat=ptransmat_tmp
    }
  

    
   if(st<4){
      v.p.it[1]<-ptransmat[1,3]
      v.p.it[3]<-ptransmat[1,4]
      v.p.it[5]<-ptransmat[1,5]
      v.p.it[6]<-ptransmat[1,6]
    }else if(st==4){v.p.it[6]<-1}
    
    #Prevent errors - if probabilities are outside the range, send to death state
    if(!(isTRUE(all.equal(sum(v.p.it),1)) & isTRUE(all.equal(sum(v.p.it<0),0)))){
      v.p.it[1]<-0
      v.p.it[3]<-0
      v.p.it[5]<-0
      v.p.it[6]<-1
      
    }
    
    
  }
  
  #print(M_it)
  #print(v.p.it)
  # Use all.equal instead of == to prevent numeric errors
  ifelse( isTRUE(all.equal(sum(v.p.it),1)) & isTRUE(all.equal(sum(v.p.it<0),0)) , return(list(v.p.it,ptransmat)), print(paste("Probabilities do not sum to 1 or negative probs present",M_it,stime,v.p.it[1:6]))) # return the transition probabilities or produce an error
}       


### Costs function
# The Costs function estimates the costs at every cycle.

#c("LN1","AEL1","PostL1","LN2","AEL2","PostL2", "PD", "Death")



Costs <- function (M_it, Trt = 5, stime=0, l1time=0, l2time=0) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  c.it <- 0                                  # by default the cost for everyone is zero 
  
  #6 cycles of L1 therapy, plust cost of management
  c.it[M_it == "LN1"]  <- cost.management.L1+
    cost.L1Cis*(Trt==1 & l1time<=5)+
    cost.L1Carbo*(Trt==2 & l1time<=5)+
    cost.L1gem*(l1time<=5)
  c.it[M_it == "AEL1"] <-cost.management.L1+
    cost.L1Cis*(Trt==1 & l1time<=5)+
    cost.L1Carbo*(Trt==2 & l1time<=5)+
    cost.L1gem*(l1time<=5) 
  c.it[M_it == "LN2"]<- cost.management.L1+
    cost.L2IO*( l2time<=24)
  c.it[M_it == "AEL2"] <-  cost.management.L1+
    cost.L2IO*( l2time<=24)
  c.it[M_it == "ED"] <- cost.management.ED
  return(c.it)        		                   # return the costs
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, Trt = 5, cl = (1/12), stime=0, l1time=0, l2time=0) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  what is the treatment 
  # cl:   cycle length (3 week cycles in years)
  
  u.it <- 0                                  # by default the utility for everyone is zero 
  u.it[M_it == "LN1"]  <- u.L1.Cis * (Trt==1 & l1time<=5) + u.L1.Carbo * (Trt==2& l1time<=5)  +
    u.L1.postTrt *(l1time>5)
  u.it[M_it == "AEL1"] <- u.L1AE 
  
  #Line 2 utilties
  u.it[M_it == "LN2"]<- u.L2.IO
  u.it[M_it == "AEL2"] <- u.L2.IO
  
  u.it[M_it == "ED"] <- u.ED
  
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}


##################################### Run the simulation ##################################


#Only run in pts who receive sequence 1
#Weight results accordingly
df.msm$wtL1<- (df.msm$L1Cis==1)/df.msm$p.trt.CisL1 + (df.msm$L1Cis==0)/(1-df.msm$p.trt.CisL1)  # Weights
samIDs<-c(1:n.i.full)[df.msm$L1Cis==1 ]
v.M_1 <- rep("LN1", length(samIDs))         # everyone begins in the healthy state 
sim_Cis  <- MicroSim(v.M_1, n.i=length(samIDs), n.t, v.n, d.c, d.e, Trt = 1, samIDs=samIDs, wt=df.msm$wtL1[samIDs]) 

#Only run in pts who receive sequence 2
#weight results accordingly
samIDs2<-c(1:n.i.full)[df.msm$L1Cis==0 ]
v.M_1 <- rep("LN1", length(samIDs2))         # everyone begins in the healthy state 

sim_Carbo <- MicroSim(v.M_1, n.i=length(samIDs2), n.t, v.n, d.c, d.e, Trt = 2, samIDs=samIDs2,  wt=df.msm$wtL1[samIDs2])  






################################# Plot state proportions by cycle (for each strategy)
#Tables summarizing counts by state
sim_Cis$TR
sim_Carbo$TR

tr.Cis<-data.frame(sim_Cis$TR)
tr.Carbo<-data.frame(sim_Carbo$TR)

L1Cis<-tr.Cis$LN1+tr.Cis$AEL1
L1Carbo<-tr.Carbo$LN1+tr.Carbo$AEL1

plot(c(0:60), L1Cis, type="l", col="red")
lines(c(0:60), L1Carbo, col="blue")

# 1. sequence 1
tr.Cis<-data.frame(sim_Cis$TR)


tr.Cis_wt<-data.frame(sim_Cis$TR_wt)

tr.Cis$cycle<-c(0:60)

mdata <- melt(tr.Cis, id=c("cycle")) 
names(mdata)<-c("cycle","state","proportion")

ggplot(data=mdata, aes(x=cycle, y=proportion, color=state,group=state)) +
  geom_line()+
  geom_point()+
  scale_color_brewer(palette="Dark2")+
  ggtitle("1. Cis")


# 2. sequence 2
tr.Carbo<-data.frame(sim_Carbo$TR)
tr.Carbo$cycle<-c(0:60)

mdata <- melt(tr.Carbo, id=c("cycle")) 
names(mdata)<-c("cycle","state","proportion")

ggplot(data=mdata, aes(x=cycle, y=proportion, color=state,group=state)) +
  geom_line()+
  geom_point()+
  scale_color_brewer(palette="Dark2")+
  ggtitle("2. Carbo")







################################# Cost-effectiveness analysis #############################



# store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v.C  <- c(sim_Cis$tc_hat_wt, sim_Carbo$tc_hat_wt)

# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)

v.E  <- c(sim_Cis$te_hat_wt, sim_Carbo$te_hat_wt)

# Summarize results for each 
res<-cbind.data.frame(v.C,v.E)
rownames(res)<-c("Cis", "Carbo")
res<-res[order(res$v.C),]

res$delta.C<-c(NA, res$v.C[2:dim(res)[1]]-res$v.C[1:(dim(res)[1]-1)])
res$delta.E<-c(NA, res$v.E[2:dim(res)[1]]-res$v.E[1:(dim(res)[1]-1)])
res$ICER<-res$delta.C/res$delta.E


print(res)

#Net monetary benefit at 100,000/QALY
NMB = res$delta.E*100000-res$delta.C
print(NMB)


################# Assess calibration #########################################

### Calibration figure
### Model based results vs IPTW overall survival curve

#OS- sequence 1
mdata<-data.frame(sim_Cis$TR_wt)
os.Cis.mod<-data.frame(sim_Cis$TR_wt)$Death
s<-survfit(Surv(OS_months_L1,Died)~1, df.msm[df.msm$L1Cis==1,], weight=wtL1)

png(paste("Resample_calibration",out.type,".png",sep=""),width = 960, height = 480)

## Komolgorov-Smirnov test
sum.s<-summary(s, times=c(1:60))$surv
trace.s<-(1-os.Cis.mod)[-1]
ks.test(sum.s, trace.s)
sum((trace.s-sum.s))^2


par(mfrow=c(1,2))
ks<-ks.test(sum.s, trace.s)
p.cis<-paste("K-S p-val =", round(ks$p.value,3))

plot(s, conf.int=FALSE, lwd=2, xlim=c(0,60), main="Overall Survival (weighted) - Cis/Gem", xlab="Months", ylab="Proportion surviving"
)
lines(x=c(0:60), y=(1-os.Cis.mod), col="slategray4", lwd=2)
text(45,.75, p.cis)
legend("topright", c("Observed - IPTW KM", "Model-based"), col=c("black","slategray4"), lty=1)
#dev.off()


#OS- sequence 2
mdata<-data.frame(sim_Carbo$TR_wt)
os.Carbo.mod<-data.frame(sim_Carbo$TR_wt)$Death
s<-survfit(Surv(OS_months_L1,Died)~1, df.msm[df.msm$L1Cis==0,], weight=wtL1)
#png("OS_Carbo_msm.png")


## Komolgorov-Smirnov test
sum.s<-summary(s, times=c(1:60))$surv
trace.s<-(1-os.Carbo.mod)[-1]
ks<-ks.test(sum.s, trace.s)
sum((trace.s-sum.s))^2
p.carbo<-paste("K-S p-val =", round(ks$p.value,3))


plot(s, conf.int=FALSE, lwd=2, xlim=c(0,60), main="Overall Survival (weighted) - Carbo/Gem", xlab="Months", ylab="Proportion surviving"
)
lines(x=c(0:60), y=(1-os.Carbo.mod), col="slategray4", lwd=2)
text(45,.75, p.carbo)
legend("topright", c("Observed - IPTW KM", "Model-based"), col=c("black","slategray4"), lty=1)
dev.off()


### Average number of cycles in each state
apply(tr.Carbo, 2, sum)
apply(tr.Cis, 2, sum)

#Save results
save.image(paste("trajMicrosimRes_",out.type,"125.RData", sep=""))

########################################### resample and use different random seeds ################


start.time <- Sys.time()
seeds<-c(101:200)

for(k in 1:length(seeds)){
  print(k)
  seed.resam<-seeds[k]
  set.seed(seed.resam)
  
  #resample with replacement from those treated with sequence 1
  samIDs<-sample(c(1:n.i.full)[df.msm$L1Cis==1 ], replace=TRUE)
  v.M_1 <- rep("LN1", length(samIDs))         # everyone begins in the healthy state 
  sim_Cis_resam  <- MicroSim(v.M_1, n.i=length(samIDs), n.t, v.n, d.c, d.e, Trt = 1, samIDs=samIDs, wt=df.msm$wtL1[samIDs],seed = seed.resam) 
  
  
  saveRDS(sim_Cis_resam,paste(out.type,"Cis",seed.resam,".RData", sep=""))
  
  set.seed(seed.resam)
  samIDs2<-sample(c(1:n.i.full)[df.msm$L1Cis!=1], replace=TRUE)
  v.M_1 <- rep("LN1", length(samIDs2))         # everyone begins in the healthy state 
  
  sim_Carbo_resam <- MicroSim(v.M_1, n.i=length(samIDs2), n.t, v.n, d.c, d.e, Trt = 2, samIDs=samIDs2,  wt=df.msm$wtL1[samIDs2],seed = seed.resam)  
  
  saveRDS(sim_Carbo_resam,paste(out.type,"Carbo",seed.resam,".RData", sep=""))
  
}


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


############################## get SE from resamples ########################################


eff.Carbo<-eff.Cis<-cost.Carbo<-cost.Cis<-rep(NA,100)

#Read in each R object from each resample and save relevant parts to a matrix
for (i in 101:200) {
  
  carbo.tmp<-readRDS(paste(out.type,"Carbo",i,".RData", sep=""))
  cis.tmp<-readRDS(paste(out.type,"Cis",i,".RData", sep=""))
  
  eff.Carbo[i-100]<-carbo.tmp$te_hat_wt
  eff.Cis[i-100]<-cis.tmp$te_hat_wt
  
  cost.Carbo[i-100]<-carbo.tmp$tc_hat_wt
  cost.Cis[i-100]<-cis.tmp$tc_hat_wt
  
  
}

#calculate results of interest
res.resam<-cbind.data.frame(cost.Carbo,cost.Cis, eff.Carbo, eff.Cis)
res.resam$cost.diff<-res.resam$cost.Cis-res.resam$cost.Carbo
res.resam$eff.diff<-res.resam$eff.Cis-res.resam$eff.Carbo
res.resam$ICER<-res.resam$cost.diff/res.resam$eff.diff
res.resam$NMB100k<- 100000*res.resam$eff.diff-res.resam$cost.diff

#Calculate SD in bootstrap resamples to get SE
apply(res.resam, 2, sd, na.rm=TRUE)


