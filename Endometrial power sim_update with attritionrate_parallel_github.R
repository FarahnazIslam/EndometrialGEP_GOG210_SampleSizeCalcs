
# packages required
library(survival)
library(contingencytables)
library(epiR)
library(parallel) 

# clear R memory 
rm(list=ls())

# data generation function
sim.fun=function(parallel_index,num.subjects,recur_rate,Sen,Spe,hazard,drop.hazard,alpha) {
  #generate survival data
  time_R = rexp(num.subjects, rate=hazard)        #time to recurrence (time_R)
  time_D = rexp(num.subjects, rate=drop.hazard)   #attrition time (time_D)
  event = time_R<=time_D & time_R<48              #4Y-recurrence event indicator (event) 
  drop.happen = time_R>time_D & time_D<48         #attrition indicator
  FU = ifelse(event, time_R, pmin(time_D, 48) )   #Follow-up time (FU)
  data = data.frame(FU, time_D, event, drop.happen)
  
  #generate disease status and GEP status
  disease_pos = ifelse(data$event, 1, 0)      #Disease status (1 for Disease + and 0 for Disease -)
  D_col_total = sum(disease_pos)              #total Disease +
  ND_col_total = num.subjects - D_col_total   #total Disease -
  TP = rbinom(D_col_total, 1, Sen)    #MR-High patients for the 4Y-recurrent patients (TP)
  TN = rbinom(ND_col_total, 1, Spe)   #MR-Low patients for the 4Y-recurrent free patients (TN) 
  index.D = which(disease_pos==1)     #indicator for Disease +
  index.ND = which(disease_pos==0)    #indicator for Disease -
  test.pos = rep(0, num.subjects)
  test.pos[index.D] = TP          #GEP status for Disease + patients
  test.pos[index.ND] = -(TN-1)    #GEP status for Disease - patients
  data = data.frame(data, disease_pos, test.pos)
  Total_N = nrow(data)      #total number of patients in data
  
  # calculate Sensitivity and Specificity using patients with full 4 year follow-up
  data.new = data[-which(data$drop.happen==1),]
  Data_mat <- table(data.new$test.pos, data.new$disease_pos)[2:1,2:1]
  result1 = epi.tests(Data_mat, method = "exact", digits = 3, conf.level = 0.95) 
  sensitivity = result1$detail[3,2]
  specificity = result1$detail[4,2]

  # perform survival rate comparisons using Z-test and variance is estimated by Greenwood formula
  fit = survfit(Surv(FU, disease_pos) ~ test.pos, error=c("greenwood"), data=data)
  surv_est = summary(fit, time=48)$surv
  surv_err = summary(fit, time=48)$std.err
  var = sum(surv_err^2)
  Z = (surv_est[1] - surv_est[2]) / (sqrt(var))
  p_value = pnorm(Z, lower.tail=FALSE)
  sig = ifelse(p_value<alpha, 1, 0)
  
  #estimate overall 4Y recurrence rate
  fit1 = survfit(Surv(FU, disease_pos) ~ 1, error=c("greenwood"), data=data)
  surv_est1 = summary(fit1, time=48)$surv 
  
  #estimate hazard ratio
  fit2 <- coxph(Surv(FU, disease_pos) ~ test.pos, data=data)
  HR = exp(fit2$coefficients)
  
  #combine and output relevant components
  result.temp=cbind(Sensitivity=sensitivity, Specificity=specificity, Total_N=Total_N,
                    recur.KM=surv_est1, HR=HR, OS_H=surv_est[2], OS_L=surv_est[1], sig=sig)
  result_one=data.frame(result.temp)
  return(result_one)
}

# Simulation assumptions
num.subjects=650    #number of subjects
recur_rate=.057   #4-year distant recurrence rate (4YR)
hazard= -log(1-recur_rate)/48   #4-year recurrence hazard rate calculation based on month
attr_rate=.15     #yearly attrition rate
drop.hazard = -log(1-attr_rate)/12   #attrition hazard rate calculation based on month
Sen=.545    #Sensitivity
Spe=.785    #Specificity
alpha=0.05  #Level of significance

# Conduct 10000 simulations and save to dataset "result"
trials=mclapply(as.list(1:10000), sim.fun, num.subjects=num.subjects, recur_rate=recur_rate, 
                Sen=Sen, Spe=Spe, hazard=hazard, drop.hazard=drop.hazard, alpha=alpha, mc.cores=8)
trials=Reduce(rbind,trials) 
result=data.frame(trials)

# Calculate averages from 10000 simulations
recur.KM = 1- round(mean(result$recur.KM),3)  #anticipated 4Y recurrence rate 
avg.sen = round(mean(result$Sensitivity),3)   #anticipated sensitivity 
avg.spe = round(mean(result$Specificity),3)   #anticipated specificity 
Total.N = mean(result$Total_N)    #sample size, N
power = round(mean(result$sig),3)   #power
OS.H = 1 - round(mean(1-result$OS_H),3)   #4Y DRFS in MR-High 
OS.L = 1 - round(mean(1-result$OS_L),3)   #4Y DRFS in MR-Low
avg.HR = round(mean(result$HR),1)     #hazard ratio

# Print simulation results table
Table = data.frame(recur.KM, avg.sen, avg.spe, Total.N, power, OS.H, OS.L, avg.HR)
colnames(Table) = c("Recurrence rate", "Sensitivity", "Specificity", "N", 
                      "Power", "DRFS MR-High", "DRFS MR-Low", "HR")
Table


