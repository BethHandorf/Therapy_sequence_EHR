########################################################################################
# This code runs a microsimulation model based on probabilities from a pre-fit #########
# multi-state model. mstate_setup_for_microsimulation.R MUST be run FIRST and the ######
# resulting .Rmd file is read in.  This code runs the  microsimulation, checks #########
# model calibration, and summarizes results.  It also re-fits the model, resampling ####
# the patient population and changing the random seed.  This should ideally be #########
# modified to run in parallel on a high-performance computing environment.  A ##########
# small number of replicates are shown here for demonstration purposes. ################
# Note there are placeholders in this code for AE states, but they are not used here. ## 
########################################################################################

library(MASS)
library(dplyr)
library(mstate)
library(survival)
library(flexsurv)
library(ggplot2)
library(reshape2)


rm(list = ls())  # remove any variables in R's memory 


############################################################################################
############################################################################################

# Microsimulation model incorporating probabilites from multi-state models

# Based on code from:


# 'Microsimulation modeling for health decision sciences using R: a tutorial' 
# Authors: Eline Krijkamp, Fernando Alarid-Escudero, 
#          Eva Enns, Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou
#
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial
#



#out.type<-"ind"
out.type<-"clayton"
#out.type<-"GaussCS"
#out.type<-"GaussUN"
#out.type<-"TCS"
#out.type<-"TUN"
#out.type<-"smallSS"
#out.type<-"LogLogisClayton"
#out.type<-"lnormClayton"

#load workspace with data, models, and functions
load(paste("msmFits_",out.type,"125.RData", sep=""))


##################################### Model input #########################################
# Model input
n.i   <- dim(df.msm)[1]                # number of simulated individuals
n.t   <-60                   # time horizon, 60 1 month cycles (5 years)

#8 states:
#LN1: line 1 therapy
#AEL1: No progression, Adverse event (line 1) 

#LN2: Line 2 (IO)
#AEL2: No progression, Adverse event (line 2) 

#ED: Extensive Disease
#D: Death

v.n   <- c("LN1","AEL1","LN2","AEL2", "ED", "Death")  
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("LN1", n.i)         # everyone begins in the healthy state 
d.c   <- d.e <- 0.03/12           # equal discounting of costs and QALYs 
#d.c   <- d.e <- 0  


#2 strategies:
# 1. Cis: Cis/Gem -> IO single agent
# 2. Carbo:  Carbo/Gem -> IO single agent

v.Trt <- c("Cis", "Carbo") # store the strategy names




#######
# Other transition probabilities (per cycle)
#######

#Placeholders for AE probs 
p.NP_Carbo = c(rep(.4/6,6),rep(0, 60-6))
p.NP_Cis = c(rep(.4/6,6),rep(0, 60-6))

p.FNP_Carbo = c(rep(.03/6,6),rep(0, 60-6))
p.FNP_Cis = c(rep(.03/6,6),rep(0, 60-6))


p.NP_FNP.Carbo.L1<-p.NP_Carbo
  p.NP_FNP.Cis.L1<-p.NP_Cis
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

#Note: these are simplified and slightly different than those used in the clinical example.  
#Most notably, we used different utilities for the two different therapies, which 
#was not the case in the clinical example.

#Line 1
u.L1.Cis<-0.67
u.L1.Carbo<-0.7
u.L1.postTrt<-0.72

u.L1AE<-0.6 

#Line 2
u.L2.IO<-0.6

#Ed
u.ED<-0.52



##################################### Microsimulation Functions ###########################################

MicroSim <- function(v.M_1, n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, Trt = 5, samIDs=c(1:2000), seed = 2) {
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
  
  #matrix of transition probabilities
  ptransmat=NULL
  for (i in 1:n.i) {
    print(paste("individualID=",samIDs[i]))
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m.C[i, 1] <- Costs(m.M[i, 1], Trt)  # estimate costs per individual of the initial health state conditional on treatment
    m.E[i, 1] <- Effs (m.M[i, 1], Trt)  # estimate QALYs per individual of the initial health state conditional on treatment
    
    
    for (t in 1:n.t) {
      #print(paste("time=",t))

      #Get the transition matrix if a transition occurred or if we don't have it
      #This step is updated to increase efficiency
      if (m.stime[i,t]==1 | is.null(ptransmat)) {
        probsrun <- Probs(m.M[i, t],stime=m.stime[i,t],mtime=t, l2time=m.l2time[i,t], pdtime=m.pdtime[i,t], Trt,L1.tot.time=v.L1.tot.time[i], ID=samIDs[i],ptransmat=ptransmat)           
      }
      
	    ptransmat<-probsrun[[2]]
      v.p<-probsrun[[1]]   # get the transition probabilities 
				   # NOTE ASSUMPTION OF EXPONENTIAL HAZARDS HERE - times from first cycle ONLY
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
  
  results <- list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR) # store the results from the simulation in a list  
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
  

  #Probabilities vary by cycle time, model time, line 2 time
  #Note: probabilities of progressing/dying are already conditional
  
  #Note - assigning values to zero unnecessary b/c vector was initialized to zero
  #Commented out zeros included only for clarity
  
  #State mapping to mstate model
  #LN1 = 1
  #LN2 = 2
  #ED = 3
  #Death = 4
  
  if(M_it=="LN1"){st=1}
  if(M_it=="AEL1"){st=1}
  if(M_it=="LN2"){st=2}
  if(M_it=="AEL2"){st=2}
  if(M_it=="ED"){st=3}
  if(M_it=="Death"){st=4}
  

  Cistrtcode=as.numeric(Trt==1)
  
  if(st<4 & (stime==1 | is.null(ptransmat))){
      ptransmat_tmp<-suppressWarnings(Prob.mat.pat(patientid = ID, time=mtime,   state=st, Cis=Cistrtcode))
      ptransmat=ptransmat_tmp
      #print(ptransmat)
  }
  

  
  #call function to get probs based on current state, mtime, and patientID
  #Only need to do this for non-death states
  if(st<4){
    Pt.id<-df.msm$PatientID[ID]
    #Get the predictions from the pre-defined function Prob.pat.t.state
    

    v.p.it[1]<-ptransmat[stime,3]
    v.p.it[3]<-ptransmat[stime,4]
    v.p.it[5]<-ptransmat[stime,5]
    v.p.it[6]<-ptransmat[stime,6]
  }else if(st==4){v.p.it[6]<-1}
  
  #Fix an error that happens occasionally - negative probabilities of staying in L2
  if(sum(v.p.it<0)>0){v.p.it<-c(0,0,0,0,0,1)}
  

# Use all.equal instead of == to prevent numeric errors
 ifelse( isTRUE(all.equal(sum(v.p.it),1)) & isTRUE(all.equal(sum(v.p.it<0),0)) , return(list(v.p.it,ptransmat)), print(paste("Probabilities do not sum to 1 or negative probs present",M_it,stime,v.p.it[1:6]))) # return the transition probabilities or produce an error
}       


### Costs function
# The Costs function estimates the costs at every cycle.

#c("LN1","AEL1","PostL1","LN2","AEL2","PostL2", "PD", "Death")


#Add costs of AEs here if using
Costs <- function (M_it, Trt = 5, stime=0, l1time=0, l2time=0) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  c.it <- 0                                  # by default the cost for everyone is zero 

  #6 cycles of L12 therapy, plus cost of management
  c.it[M_it == "LN1"]  <- cost.management.L1+
                          cost.L1Cis*(Trt==1 & l1time<=5)+
                          cost.L1Carbo*(Trt==2 & l1time<=5)+
                          cost.L1gem*(l1time<=5)

  c.it[M_it == "LN2"]<- cost.management.L2+
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
  u.it[M_it == "AEL1"] <- u.L1AE #Weighted utility of FNP and NP
   
  #Line 2 utilties
  u.it[M_it == "LN2"]<- u.L2.IO
  u.it[M_it == "AEL2"] <- u.L2.IO

  u.it[M_it == "ED"] <- u.ED
  
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}


##################################### Run the microsimulation ##################################
start.time <- Sys.time()
sim_Cis  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 1)  # sequence 1

sim_Carbo <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 2)  # sequence 2
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


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

# 1. Sequence 1
tr.Cis<-data.frame(sim_Cis$TR)
tr.Cis$cycle<-c(0:60)

mdata <- melt(tr.Cis, id=c("cycle")) 
names(mdata)<-c("cycle","state","proportion")



ggplot(data=mdata, aes(x=cycle, y=proportion, color=state,group=state)) +
  geom_line()+
  geom_point()+
  scale_color_brewer(palette="Dark2")+
  ggtitle("1. Cis")


# 2. Sequence 2
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
v.C  <- c(sim_Cis$tc_hat, sim_Carbo$tc_hat)
sd.C <- c(sd(sim_Cis$tc), sd(sim_Carbo$tc))/  sqrt(n.i)

# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)

v.E  <- c(sim_Cis$te_hat, sim_Carbo$te_hat)
sd.E <- c(sd(sim_Cis$te), sd(sim_Carbo$te))/ sqrt(n.i)


# Summarize results for each 
res<-cbind.data.frame(v.C,sd.C,v.E,sd.E)
rownames(res)<-c("Cis", "Carbo")
res<-res[order(res$v.C),]


res$delta.C<-c(NA, res$v.C[2:dim(res)[1]]-res$v.C[1:(dim(res)[1]-1)])
res$delta.E<-c(NA, res$v.E[2:dim(res)[1]]-res$v.E[1:(dim(res)[1]-1)])
res$ICER<-res$delta.C/res$delta.E


print(res)

NMB = res$delta.E*100000-res$delta.C
print(NMB)


################# Assess calibration #########################################

trtmod<-glm(L1Cis~Surgery+AgeAdvDx+ I(AgeAdvDx^2)+       
              +EcogValue+albumin, data=df.msm, family="binomial")

summary(trtmod)
#propensity scores and weights

df.msm$p.trt.CisL1<-stats:::predict(trtmod, data=trtmod$data, type="response")
df.msm$wtL1<-df.msm$L1Cis/df.msm$p.trt.CisL1 + (1-df.msm$L1Cis)/(1-df.msm$p.trt.CisL1)

df.msm$wtL1carbo<-df.msm$L1Cis/(1-df.msm$p.trt.CisL1) + (1-df.msm$L1Cis)/(df.msm$p.trt.CisL1)

#General calibration
s<-survfit(Surv(OS_months_L1,Died)~L1Cis, df.msm, weight=wtL1)
summary(s, times=c(24,60))

s<-survfit(Surv(PFS_months,Progression)~L1Cis, df.msm, weight=wtL1)
summary(s, times=c(24,60))




######## plots of calibration ########################

png(paste("MSM_calibration",out.type,".png",sep=""),width = 960, height = 480)
par(mfrow=c(1,2))

#OS- sequence 1

mdata<-tr.Cis
os.Cis.mod<-tr.Cis$Death
s<-survfit(Surv(OS_months_L1,Died)~1, df.msm[df.msm$L1Cis==1,], weight=wtL1)

## Komolgorov-Smirnov test
sum.s<-summary(s, times=c(1:60))$surv
trace.s<-(1-os.Cis.mod)[-1]
ks<-ks.test(sum.s, trace.s)
p.cis<-paste("K-S p-val =", round(ks$p.value,3))


plot(s, conf.int=FALSE, lwd=2, xlim=c(0,60), main="Overall Survival - sequence 1", xlab="Months", ylab="Proportion surviving"
     )
lines(x=c(0:60), y=(1-os.Cis.mod), col="slategray4", lwd=2)
text(45,.75, p.cis)
legend("topright", c("Observed - IPTW KM", "Model-based"), col=c("black","slategray4"), lty=1)


#OS - sequence 2

mdata<-tr.Carbo
os.Carbo.mod<-tr.Carbo$Death
s<-survfit(Surv(OS_months_L1,Died)~1, df.msm[df.msm$L1Cis==0,], weight=wtL1)


## Komolgorov-Smirnov test
sum.s<-summary(s, times=c(1:60))$surv
trace.s<-(1-os.Carbo.mod)[-1]
ks<-ks.test(sum.s, trace.s)
p.carbo<-paste("K-S p-val =", round(ks$p.value,3))

plot(s, conf.int=FALSE, lwd=2, xlim=c(0,60), main="Overall Survival - sequence 2", xlab="Months", ylab="Proportion surviving"
     )
lines(x=c(0:60), y=(1-os.Carbo.mod), col="slategray4", lwd=2)
text(45,.75, p.carbo)
legend("topright", c("Observed - IPTW KM", "Model-based"), col=c("black","slategray4"), lty=1)
dev.off()


##### Save .Rmd file with results

save.image(paste("msmMicrosimRes_",out.type,"125.RData", sep=""))

########################################### resample and use different random seeds ################


# using 5 for demonstration increase number of repetitions to 50-100.  
# Note: this should ideally be done in parallel

for (i in 1:5) {

  seed_msim=100+i
  set.seed(seed_msim)

  #resample with replacement
  samIDs<-sample(c(1:n.i), replace=TRUE)

  sim_Cis  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 1, samIDs=samIDs, seed=seed_msim)  # 1. DCT -> AB

  sim_Carbo <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 2, samIDs=samIDs, seed=seed_msim)

  saveRDS(sim_Cis,paste(out.type,"seq1",seed_msim,".RData", sep=""))
  saveRDS(sim_Carbo,paste(out.type,"seq2",seed_msim,".RData", sep=""))
  
  
}


############################## get SE from resamples ########################################

#Read in each R object from each resample and save relevant parts to a matrix
for (i in 101:105) { #Change if using more samples
  
  if(file.exists(paste(out.type,"seq1",i,".RData", sep=""))){
    carbo.tmp<-readRDS(paste(out.type,"seq1",i,".RData", sep=""))
    cis.tmp<-readRDS(paste(out.type,"seq2",i,".RData", sep=""))
    
    eff.Carbo[i-100]<-carbo.tmp$te_hat
    eff.Cis[i-100]<-cis.tmp$te_hat
    
    cost.Carbo[i-100]<-carbo.tmp$tc_hat
    cost.Cis[i-100]<-cis.tmp$tc_hat}

}

#calculate results of interest
res.resam<-cbind.data.frame(cost.Carbo,cost.Cis, eff.Carbo, eff.Cis)
res.resam$cost.diff<-res.resam$cost.Cis-res.resam$cost.Carbo
res.resam$eff.diff<-res.resam$eff.Cis-res.resam$eff.Carbo
res.resam$ICER<-res.resam$cost.diff/res.resam$eff.diff
res.resam$NMB100k<- 100000*res.resam$eff.diff-res.resam$cost.diff

#Calculate SD in bootstrap resamples to get SE
apply(res.resam, 2, sd, na.rm=TRUE)