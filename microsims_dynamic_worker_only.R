#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Rscript microsims_dynamic_worker_only.R <n_sims> <risk.multiplier> <testingDelay> <alpha.late> <sensitivity type> <offset>

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Title: Universal testing to reduce the risk of returning to work: a simulation analysis 
# Code author: Elizabeth Chin (Stanford) etchin@stanford.edu
# Origin date: 4/17/20
# Last updated: 6/2/20
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Community matrix- microsimulation 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Columns
#   1 - ID: Person number
#   2 - InfectionType: If gets infection, what type of infection will it be: 0: Subclinical; 1: Clinical
#   3 - State: 0: Susceptible; 1: Exposed; 2: Early Subclinical; 3: Late Subclinical; 4: Early Clincal; 5: Late Clinical; 6: Recovered
#   4 - DaysExposed: # Days, Exposure
#   5 - DaysEarlyInf: # Days, Presymtpomatic infectious
#   6 - DaysLateInf: # Days, Infectious
#   7 - NewInfection: Whether infected that day
#   8 - DayOfInfection: Day of initial infection
#   9 - DayInTrue: Day in TRUE state (Number of days in infections)
#   10 - DayInObs: Day in OBSERVED state (Number of off work days)
#   11 - TestDays: if Day - TestDays % freq == 0 --> individual is tested
#   12 - InfWorkDays: Number of infectious work days

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulation constants / distributions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

states <- c("S","E","ESC","LSC","EC","LC","R")
n.states <- length(states)

quarantineDays <- 14

n.contacts <- 10 # daily contacts
days.earlyInfection <- 2
p.a = 0.4 # proportion asymptomatic
n_reps <- as.integer(as.numeric(args[1]))
risk.multi <- as.numeric(args[2])
R0.ec <- 2.2
beta.ec <- R0.ec/(n.contacts*days.earlyInfection)
nDelay <- as.numeric(args[3])
alpha.a <- 0.5
alpha.late <- as.numeric(args[4])
sens_type <- args[5]
offset <- as.integer(args[6])

# incubation parameters from
# Li, Q. et al. Early transmission dynamics in Wuhan, China, of novel coronavirus-infected pneumonia. N. Engl. J. Med. 382, 1199–1207 (2020).
# 5.2 days (95% CI, 4.1–6.4 days) fit to a log normal distribution
incubation_fxn <- function(n) return(rtnorm(n, mean=5.2, sd=1, lower=4.1, upper=6.4)-2.3) # log normal parameters for days of incubation

# infectious before symptom onset parameters from
# He X et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nature Medicine. 15 April 2020.
# 2.3 days (95% CI, 0.8–3.0 days) fit to gamma distribution
early_infectious_fxn <- function(n) return(rtnorm(n, mean=2.3, sd=1, lower=.8, upper=3)) # gamma parameters for days before symptom onset

# late infection number of days
#  He X et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nature Medicine. 15 April 2020.
# 7 days. We use a truncated normal distribution between 3-10
late_infectious_fxn <- function(n) return(rtnorm(n, mean=7, sd=3, lower=3, upper=10))

# community seir model
seir <- function(times, Y, parms){
  with(parms,{
    S<-Y[1]
    E<-Y[2]
    ESC<-Y[3]
    LSC<-Y[4]
    EC<-Y[5]
    LC<-Y[6]
    R <- Y[7]
    dY <- vector(length=length(Y))
    dY[1] <- -beta.t*n.contacts*(alpha.a*(ESC + alpha.late*LSC) + EC + alpha.late*LC)/sum(Y)*S #S
    dY[2] <-  -dY[1] - E/days.incubation #E
    dY[3] <- p.a/days.incubation*E - ESC/days.earlyInfection #ESC
    dY[4] <- ESC/days.earlyInfection - LSC/days.lateInfection #LSC
    dY[5] <- (1-p.a)/days.incubation*E - EC/days.earlyInfection  #EC
    dY[6] <- EC/days.earlyInfection - LC/days.lateInfection #LC
    dY[7] <- LSC/days.lateInfection + LC/days.lateInfection #R
    return(list(dY))
  })
}

simulate_community <- function(cparams){
  with(cparams,{
    initials <- c(S = n.community-1, E = 1, ESC = 0, LSC = 0, EC = 0, LC = 0, R = 0)
    
    sol = lsoda(initials, seq(0,n.days-1), seir, cparams)
    return(sol)
  })
}

update_state <- function(x){
  x[["DayInTrue"]] <- x[["DayInTrue"]] + 1
  if(x[["State"]] == 1 & x[["DayInTrue"]] > x[["DaysExposed"]]){ #if exposed
    if(x[["InfectionType"]] == 0) return(c(2, x[["DayInTrue"]]-x[["DaysExposed"]])) #subclinical
    else return(c(4, x[["DayInTrue"]]-x[["DaysExposed"]])) #clinical
  } 
  if(x[["State"]] %in% c(2,4) & x[["DayInTrue"]]  > x[["DaysEarlyInf"]] ) return(c(x[["State"]]+1, x[["DayInTrue"]]-x[["DaysEarlyInf"]]))
  if(x[["State"]] %in% c(3,5) & x[["DayInTrue"]]  > x[["DaysLateInf"]]) return(c(6, x[["DayInTrue"]]-x[["DaysLateInf"]]))
  return(c(x[["State"]], x[["DayInTrue"]] ))
}

day_in_infection <- function(x){
  if(x[["State"]] %in% c(2,4)) return(x[["DayInTrue"]])
  if(x[["State"]] %in% c(3,5)) return(x[["DayInTrue"]] + x[["DaysEarlyInf"]])
  return(NA)
}

calculate_new_worker_infections <- function(c.inf, w.inf, params){
  #c.inf <- split_arrays(c.inf)
  w <- worker_calculate_infectiousness(c.inf[1], c.inf[2], c.inf[3], c.inf[4], c.inf[5],
                                       w.inf[1], w.inf[2], w.inf[3], w.inf[4], w.inf[5],
                                       params)
  return(w)
}

#calculate for all workers
worker_calculate_infectiousness <- function(community.esc, community.lsc, community.ec, community.lc, community.n,
                                            workers.esc, workers.lsc, workers.ec, workers.lc, workers.n,
                                            params){
  with(params,{
    #assume 75% of interactions come from patients and 25% comes from other workers
    asx <- sapply(1:m, function(i) 0.75*sum(lambda.cw*(community.esc + alpha.late*community.lsc)/community.n, na.rm = TRUE) + 
                    0.25*lambda.ww*(workers.esc + alpha.late*workers.lsc)/workers.n)
    sx <- sapply(1:m, function(i) 0.75*sum(lambda.cw*(community.ec + alpha.late*community.lc)/community.n, na.rm = TRUE) + 
                   0.25*lambda.ww*(workers.ec + alpha.late*workers.lc)/workers.n)
    asx[is.na(asx)] <- 0
    sx[is.na(sx)] <- 0
    p <- alpha.a*asx + sx
    return(p)
  })
}

testing_sim <- function(sim_pop,
                        t, #time point in simulation
                        params_sim
){
  with(params_sim,{
    #New infections
    sim_pop$NewInfection <- 0
    
    n.w.esc <- sim_pop %>% filter(State == 2, DayInObs == 0) %>% nrow() # Early subclinical
    n.w.lsc <- sim_pop %>% filter(State == 3, DayInObs == 0) %>% nrow() # Late subclinical
    n.w.ec <- sim_pop %>% filter(State == 4, DayInObs == 0) %>% nrow() # Early Clinical
    n.w <- sim_pop %>% filter(State != 5, DayInObs == 0) %>% nrow() # not symptomatic or off work
    
    workers.inf <- c(n.w.esc, n.w.lsc, n.w.ec, 0, n.w)
    
    community.inf <- unlist(community_state[t,c("ESC","LSC","EC","LSC")])
    community.inf <- c(community.inf, sum(community_state[t,]))
    
    lambda.ww <- beta.t*n.contacts
    lambda.cw <- lambda.ww*risk.multi
    
    parameters <- list(
      lambda.cw = lambda.cw,
      lambda.ww = lambda.ww,
      p.a = p.a,
      m = 1
    )
    
    # S --> E
    # Workers
    sWorkers <- sim_pop %>% dplyr::filter(State == 0, DayInObs == 0) %>% pull(ID) %>% as.character()
    pW <- calculate_new_worker_infections(community.inf, workers.inf, parameters)
    sim_pop[sWorkers, "NewInfection"] <- rbinom(length(sWorkers), 1, pW)
    
    # Update infectious timelines
    newInf <- sim_pop %>%
      filter(NewInfection == 1) %>%
      pull(ID) %>%
      as.character()
    sim_pop[newInf,"State"] <- 1
    sim_pop[newInf,"DayInTrue"] <- 0
    sim_pop[newInf,"DayOfInfection"] <- t
    
    # Update states
    tmp <- t(apply(sim_pop, 1, update_state))
    sim_pop$State <- tmp[,1]
    sim_pop$DayInTrue <- tmp[,2]
    
    # Test workers
    sTested <- sim_pop %>% 
      dplyr::filter((t - TestDays) %% testingFreq == 0,
                    DayInObs == 0)
    sTested_pos_sc <- sTested %>%
      dplyr::filter(State %in% c(2,3))
    sTested_pos_c <- sTested %>%
      dplyr::filter(State %in% c(4))
    sTested_neg <- sTested %>%
      dplyr::filter(State %in% c(0,1,6)) %>%
      pull(ID) %>% as.character()
    
    sTested_fp <- rbinom(length(sTested_neg), 1, 1-p_spec)
    sTested_tp <- c()
    sPositive_ids_sc <- c()
    sPositive_ids_c <- c()
    if(nrow(sTested_pos_sc) > 0){
      day_i_sc <- ceiling(apply(sTested_pos_sc, 1, function(y) day_in_infection(y))) #rbinom(length(sTested_pos), 1, p_sens))
      sTested_tp_sc <- sapply(day_i_sc, function(y) ifelse(y > 0, rbinom(1,1,p_sens[y]), 0)) 
      sPositive_ids_sc <- as.character(sTested_pos_sc$ID)[sTested_tp_sc == 1]
    }
    if(nrow(sTested_pos_c) > 0){
      day_i_c <- ceiling(apply(sTested_pos_c, 1, function(y) day_in_infection(y)))
      sTested_tp_c <- sapply(day_i_c, function(y) ifelse(y > 0, rbinom(1,1,p_sens[y]), 0)) 
      sPositive_ids_c <- as.character(sTested_pos_c$ID)[sTested_tp_c == 1]
    }
    sPositive_ids <- c(sPositive_ids_sc, sPositive_ids_c)
    
    # Update infectious working days
    sInfWorkers <- sim_pop %>%
      dplyr::filter(DayInObs == 0, 
                    State %in% c(2,3,4)) %>%
      pull(ID) %>% as.character()
    sInfWorkers_notPos <- setdiff(sInfWorkers, sPositive_ids)
    sim_pop[sInfWorkers_notPos,"InfWorkDays"] <- sim_pop[sInfWorkers_notPos,"InfWorkDays"] + 1
    sim_pop[sPositive_ids_sc, "InfWorkDays"] <- sim_pop[sPositive_ids_sc, "InfWorkDays"] + nDelay
    sim_pop[sPositive_ids_c, "InfWorkDays"] <- sim_pop[sPositive_ids_c, ] %>%
      mutate(inf.max = if_else(InfWorkDays + nDelay > DaysEarlyInf, InfWorkDays, InfWorkDays + nDelay)) %>%
      pull(inf.max)
    
    # Update quarantine day
    sim_pop <- sim_pop %>%
      mutate(DayInObs = ifelse(DayInObs > 0 & (DayInObs <= quarantineDays), DayInObs + 1, 0))
    sim_pop[c(sPositive_ids, sTested_neg[sTested_fp == 1]), "DayInObs"] <- 1 #update for newly tested people
    
    return(sim_pop)
  })
}

run_sim <- function(sparams){
  with(sparams,{
    #Sensitivity & Specificity
    #p_sens <-rtruncnorm(n=1, a=.68, b=.80, mean=.75, sd=.2)
    p_spec <- rtruncnorm(n=1, a=.98, b=1, mean=.99, sd=.05)
    
    community_state = simulate_community(list(n.days = n.days,
                                               days.incubation = days.incubation,
                                               days.earlyInfection = days.earlyInfection,
                                               days.lateInfection = days.lateInfection,
                                               beta.t = beta.t,
                                               n.contacts = n.contacts,
                                               alpha.a = alpha.a,
                                               alpha.late = alpha.late,
                                               p.a = p.a,
                                               n.community = n.community
                                               )) %>% as.data.frame() %>% select(-`time`)
    
    #C=S(1xex(Ipqt=AV) ) calculate per day
    
    col_fields <- c("ID", "InfectionType", "State",
                    "DaysExposed","DaysEarlyInf","DaysLateInf",
                    "NewInfection","DayOfInfection","DayInTrue", "DayInObs",
                    "TestDays", "InfWorkDays")
    
    sim_pop <- as.data.frame(matrix(0,nrow=N_pop, ncol=length(col_fields)))
    colnames(sim_pop) <- col_fields
    
    sim_pop$ID <- seq(1:N_pop)
    rownames(sim_pop) <- sim_pop$ID
    sim_pop$InfectionType <- rbinom(N_pop, 1, 1-p.a)
    sim_pop$State <- c(rep(0, N_pop))
    sim_pop$DaysExposed <- days.incubation
    sim_pop$DaysEarlyInf <- days.earlyInfection
    sim_pop$DaysLateInf <- days.lateInfection
    sim_pop$TestDays <- sample(0:(testingFreq - 1), N_pop, replace = TRUE)
    
    for(j in 1:n.days){
      sim_pop <- testing_sim(sim_pop,
                             j,
                             list(
                               beta.t = beta.t,
                               n.contacts = n.contacts,
                               alpha.a = alpha.a,
                               p.a = p.a,
                               quarantineDays = quarantineDays,
                               testingFreq = testingFreq,
                               p_sens = p_sens,
                               p_spec = p_spec,
                               nDelay = nDelay,
                               community_state = community_state
                             ))
    }
    write.csv(sim_pop, paste0("data/workers/sim_pop_", risk.multi, "_", nDelay, "_", alpha.a, "_", testingFreq,"_", sens_type, "_", r,".csv"))
    return(sum(sim_pop$InfWorkDays, na.rm = TRUE))
  })
}

tf_array <- c(1:5, 7, seq(10,30,5), 1000)
tfMatrix <- matrix(0, nrow=length(tf_array),ncol=n_reps)
rownames(tfMatrix) <- as.character(tf_array)

sens_by_day <- read.csv("data/sens_by_day_ci.csv", header = TRUE) %>%
  select(fnr_med, fnr_lb, fnr_ub)

sens_by_day[1:7,] <- sens_by_day[8,]
sens_by_day <- sens_by_day[3:nrow(sens_by_day),] # only able to get a positive test during early infectious stage at the earliest

if(sens_type == "upper") p_fn <- sens_by_day$fnr_ub/100
if(sens_type == "median") p_fn <- sens_by_day$fnr_med/100
if(sens_type == "lower") p_fn <- sens_by_day$fnr_lb/100
if(sens_type == "perfect") p_fn <- rep(0,nrow(sens_by_day))

plan(multiprocess, workers = availableCores(), gc = TRUE) ## Parallelize using 15 processes

#Infectious parameters
days.incubation.array <- incubation_fxn(n_reps)
days.earlyInfection.array <- early_infectious_fxn(n_reps)
days.lateInfection.array <- late_infectious_fxn(n_reps)

for(tf in tf_array){
  print(paste0("tf: ",tf))
  x <- future_sapply(1:n_reps, function(r) run_sim(
    list(
      n.days = 300,
      beta.t = beta.ec,
      alpha.a = alpha.a,
      n.contacts = n.contacts,
      quarantineDays = quarantineDays,
      testingFreq = tf,
      p.a = p.a,
      risk.multi = risk.multi,
      nDelay = nDelay,
      N_pop = 100,
      n.community = 100000,
      p_sens = 1-p_fn,
      sens_type = sens_type,
      r = r + n_reps*offset,
      days.incubation = days.incubation.array[r],
      days.earlyInfection = days.earlyInfection.array[r],
      days.lateInfection = days.lateInfection.array[r]
    )))
  tfMatrix[as.character(tf), ] <- x
}

future:::ClusterRegistry("stop")
