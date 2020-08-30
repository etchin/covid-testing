#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Rscript microsims.R <risk> <delay> <sensitivity> <sensitivity  multiplier>

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Title: Universal testing to reduce the risk of returning to work: a simulation analysis 
# Code author: Elizabeth Chin (Stanford) etchin@stanford.edu
# Origin date: 4/17/20
# Last updated: 8/30/20
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Study background  
# 
# Objective: Compare universal testing strategies to reduce risk of returning to work 
# 
# Study outcomes:
# 1) Proportion of diagnosed infections 
# 2) Average number of infectious days without diagnosis 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Load packages
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
library(truncnorm)
library(matrixStats)
library(readr)
library(tidyr)
library(dplyr)
library(future.apply) # parallelizes
library(msm)
library(stringr)
library(deSolve)
source("R/initialize_params.R")
source("R/utils.R")
source("R/testing.R")
source("R/run_simulation.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Community matrix- microsimulation 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Columns
#   1 - ID: Person number
#   2 - InfectionType: If gets infection, what type of infection will it be: 0: Subclinical; 1: Clinical
#   3 - State: 0: Susceptible; 1: Incubation; 2: Symptomatic; 3: Recovered;
#   4 - DaysIncubation: # Days, Incubation
#   5 - DaysSymptomatic: # Days, Symptomatic Infectious
#   6 - NewInfection: Whether infected that day
#   7 - DayOfInfection: Day of initial infection
#   8 - DayOfDetection: Day of infection detection (PCR positive or symptom detection)
#   9 - DayInObs: Day in OBSERVED state (Number of off work days)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulation constants / distributions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
arg_names <- c("risk","delay","sensitivity","sens_multiplier")
risk.multi <- as.numeric(args[1])
n.delay <- as.numeric(args[2])
sens_type <- args[3]
alpha.sens <- as.numeric(args[4]) # sensitivity multiplier 

tf_array <- c(1:5, 7, seq(10,30,5), NA)
tfMatrix <- as.data.frame(matrix(0, ncol=length(tf_array)*2,nrow=n_reps))

get.seed.alpha <- function(x) {
  require("digest")
  hexval <- paste0("0x",digest(x,"crc32"))
  intval <- type.convert(hexval) %% .Machine$integer.max
  return(intval)
}

sys_vars <- Sys.getenv(c("SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID"))


plan(multiprocess, workers = availableCores()-1, gc = TRUE) ## Parallelize using 15 processes

#Infectious parameters
days.incubation.array <- incubation_fxn(n_reps)
days.symptomatic.array <- symptomatic_fxn(n_reps)

counter <- 0
for(tf in tf_array){
  print(paste0("tf: ",tf))
  x <- future_sapply(1:n_reps, function(r) run_sim(
    list(
      n.days = 300,
      beta.t = beta.t,
      alpha.a = alpha.a,
      n.contacts = n.contacts,
      days.quarantine = days.quarantine,
      test.freq = tf,
      p.a = p.a,
      risk.multi = risk.multi,
      n.delay = n.delay,
      N_pop = 100,
      p.sens = alpha.sens * load_test_sensitivity(sens_type),
      r = r,
      days.incubation = days.incubation.array[r],
      days.symptomatic = days.symptomatic.array[r],
      p.community = p.community
    )))
  tfMatrix[,(2*counter+1):(2*(counter+1))] <- t(x)
  colnames(tfMatrix)[(2*counter+1):(2*(counter+1))] <- paste0(c("Infectiousness","nInfected"),"_",tf)
  counter <- counter + 1
}

comment <- paste0("#", paste(arg_names,args[1:length(arg_names)], sep=":", collapse=","))

out_fi <- paste0("output/sim_pop_", sys_vars["SLURM_JOB_ID"], "_",
                 sys_vars["SLURM_ARRAY_TASK_ID"], ".csv")
print(paste0("Writing to ", out_fi))
con <- file(out_fi, open="wt")
writeLines(comment, con)
write.csv(tfMatrix, con)
close(con)


future:::ClusterRegistry("stop")
