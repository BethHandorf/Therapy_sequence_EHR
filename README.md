# Therapy_sequence_EHR
Mircrosimulations to analyze cost-effectiveness of therapy sequences informed by EHR data  
  

This repository contains code to support the manuscript:  
Modeling therapy sequence for advanced cancer: A microsimulation approach leveraging Electronic Health Record data 

Acknowledgement: this project uses code adapted from
 'Microsimulation modeling for health decision sciences using R: a tutorial' 
 Authors: Eline Krijkamp, Fernando Alarid-Escudero, 
          Eva Enns, Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou
 See GitHub for more information or code updates
 https://github.com/DARTH-git/Microsimulation-tutorial

# Contents

## 1. Synthetic EHR data
Data dictionary for synthetic patient data (created to mimic EHR-derived variables)	
Note: this is fully simulated data, no patient data is present	
	
ID	numeric identifier for simulated patient (sequential numbers)
	
Covariates: These are simulated covariates based on the observed distribution of covariates (and their correlations) in a real, observed dataset	
	
Surgery	- binary indicator for receipt of surgery (1= surgery, 0= no surgery) (binomial distribution)  
Male	- binary indicator for sex (1=Male, 0=Female)  (binomial distribution)  
RaceB	- binary indicator for race (1=Black, 0=White or Other)   (binomial distribution)  
RaceO	- binary indicator for race (1=Other, 0=Black or White)  (binomial distribution)  
Hispanic	- binary indicator for Hispanic ethnicity (1=Hispanic, 0=non-Hispanic)  (binomial distribution)  
AgeAdvDx	- Age (in years) at advanced diagnosis of urothelial cancer (transformed Beta distribution)  
EcogValue	- Eastern Cooperative Oncology Group performance status score (categorized as 0=0, 1=1, 2=2+) (multinomial distribution)  
MDRD_GFR	- Estimated Glomerular Filtration rate at baseline (updated equation from the Modification to Diet in Renal Disease study) (Gamma distribution)  
hemoglobin	- lab-based hemoglobin at baseline (Gamma distribution)  
albumin	- lab-based albumin at baseline (skew-Normal distribution)  
CHF_known - 	binary indicator for Chronic heart failure at baseline (1=known CHF diagnosis, 0= no known diagnosis)  (binomial distribution)  
hear_known	- binary indicator for hearing loss at baseline (1= known hearning loss, 0= no known hearing loss)  (binomial distribution)  
neuropathy_known	- binary indicator for neuropathy at baseline (1= known neuropathy, 0= no known neuropathy)  (binomial distribution)  
viscmet_known	- binary indicator for visceral metastases at baseline (1= known visceral mets, 0 = no known visceral mets)  (binomial distribution)  
	
Treatments: These are predicted probabilities of treatments based on the simulated covariates and a propensity-score model fit on the true underlying dataset	
	
p.trt.CisL1	- probability of receiving Cisplatin/Gemcitabine first-line therapy
pr.L2.0	- probability of receiving no second-line therapy
pr.L2.1	- probabiilty of receiving non-IO second-line therapy
pr.L2.2	- probability of receiving IO second-line therapy

## 2. Code to simulate outcomes

## 3. Pre-fit multi-state models

## 4. Run microsimulation with multi-state models

## 4. Run microsimulation with patient trajectories
