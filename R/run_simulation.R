require(matrixStats)
require(readr)
require(tidyr)
require(dplyr)
require(data.table)

source("R/utils.R")
source("R/transmission.R")
source("R/testing.R")

testing_sim <- function(sim_pop,
                        t, #time point in simulation
                        params_sim
){
  with(params_sim,{
    sim_pop[DayOfDetection == t, DayInObs := 1]
    
    #New infections
    sim_pop[, NewInfection := 0]
    
    n.w.asx <- 0
    n.w.sx <- 0
    n.w.i <- sim_pop[State %in% c(1,2) & DayInObs == 0]
    if(nrow(n.w.i) > 0){
      n.w.i$infectiousness <- sapply(n.w.i$DayOfInfection, function(x) calc_infectiousness(t - x + 1, if_weights))
      setindex(n.w.i, NULL)
      n.w.asx <- n.w.i[InfectionType == 0, infectiousness] # Asymptomatic
      n.w.sx <-  n.w.i[InfectionType == 1 & State == 1, infectiousness] # No symptomatic workers
      if(length(n.w.asx) > 0) n.w.asx <- sum(n.w.asx, na.rm = TRUE) else n.w.asx <- 0
      if(length(n.w.sx) > 0) n.w.sx <- sum(n.w.sx, na.rm = TRUE) else n.w.sx <- 0
    }
    n.w <- sim_pop[DayInObs == 0, .N] # not symptomatic or off work
    workers.inf <- c(n.w.asx, n.w.sx, n.w)
    
    #define infectious period as -4 days to 7 days
    c.days.symptomatic <- 7
    c.days.incubation <- 4
    start.incubation <- days.incubation - c.days.incubation
    if(start.incubation < 1) start.incubation <- 1
    if_weights_c <- if_weights[start.incubation:(start.incubation + days.incubation + c.days.symptomatic - 1)]
    if_weights_c <- if_weights_c / sum(if_weights_c)
    n.community <- 100000
    c.inf <- sample(if_weights_c, sum(rbinom(n.community, 1, p.community)), replace = TRUE)
    c.asx <- rbinom(length(c.inf), 1, p.a)
    community.asx <- sum(c.inf[c.asx == 1])
    community.sx <- c.inf[c.asx == 0]
    community.sx1 <- community.sx[grepl('e', names(community.sx))]
    community.sx2 <- community.sx[grepl('l', names(community.sx))]
    community.sx <- sum(c(sum(community.sx1), alpha.late.c*sum(community.sx2)), na.rm = TRUE)
    community.inf <- c(community.asx, community.sx, n.community) / n.community
    
    lambda.ww <- beta.t*n.contacts
    lambda.cw <- lambda.ww*risk.multi
    
    iparams <- list(
      lambda.cw = lambda.cw,
      lambda.ww = lambda.ww,
      p.a = p.a,
      alpha.a = alpha.a,
      m = 1
    )
    
    # S --> E
    # Workers
    sWorkers <- sim_pop[State == 0 & DayInObs == 0, ID]
    pW <- calculate_new_worker_infections(community.inf, workers.inf, iparams)
    sim_pop[sWorkers, NewInfection := rbinom(length(sWorkers), 1, pW)]
    
    # Update states
    sim_pop[, State := apply(sim_pop, 1, function(x) update_state(x, t))]
    
    # Update infectious timelines
    sim_pop[NewInfection==1, c("State","DayOfInfection") := list(1, t + 1)]
    
    # Update quarantine day
    sim_pop[, DayInObs := ifelse(DayInObs <= days.quarantine & DayInObs > 0, DayInObs + 1, 0)]
    
    return(sim_pop)
  })
}

run_sim <- function(sparams){
  with(sparams,{
    #Sensitivity & Specificity
    p.spec <- rtruncnorm(n=1, a=.98, b=1, mean=.99, sd=.05)
    
    col_fields <- c("ID", "InfectionType", "State",
                    "DaysIncubation","DaysSymptomatic",
                    "NewInfection","DayOfInfection","DayInObs",
                    "DayOfDetection","DayOfDetection_iwd")
    
    sim_pop <- as.data.frame(matrix(0,nrow=N_pop, ncol=length(col_fields)))
    colnames(sim_pop) <- col_fields
    
    sim_pop$ID <- as.character(seq(1:N_pop))
    rownames(sim_pop) <- sim_pop$ID
    sim_pop$InfectionType <- rbinom(N_pop, 1, 1-p.a)
    sim_pop$State <- c(rep(0, N_pop))
    sim_pop$DaysIncubation <- days.incubation
    sim_pop$DaysSymptomatic <- days.symptomatic
    
    if_weights <- infectious_fxn(days.incubation, days.symptomatic)

    if_e <- sum(if_weights[1:days.incubation])
    if_l <- sum(tail(if_weights,days.symptomatic))
    multi.i <- days.symptomatic*if_l / (days.incubation*if_e)
    
    sim_pop <- as.data.table(sim_pop)
    setkey(sim_pop, ID)
    
    if(!is.na(test.freq)){
      test.start <- sample(1:test.freq, 1)
      test.days <- seq(from = test.start, to = n.days, by = test.freq)
    } else{
      test.days <- NA
    }
    
    tparams <- list(
      p.sens = p.sens,
      p.spec = p.spec,
      n.delay = n.delay,
      days.incubation = days.incubation,
      days.symptomatic = days.symptomatic
    )
    
    for(j in 1:n.days){
      # symptom screening
      sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2, DayOfDetection := j]
      sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2 & DayOfDetection_iwd == 0,
              DayOfDetection_iwd := j]
      
      if(j %in% test.days){
        # Test workers
        wids <- sim_pop[DayInObs == 0 & DayOfDetection < j, ID]
        sim_pop[wids,] <- test_pop(sim_pop[wids,], j, tparams)
      }
      
      sim_pop <- testing_sim(sim_pop,
                             j,
                             list(
                               beta.t = beta.t,
                               n.contacts = n.contacts,
                               alpha.a = alpha.a,
                               p.a = p.a,
                               days.quarantine = days.quarantine,
                               if_weights = if_weights,
                               p.community = p.community,
                               days.incubation = days.incubation
                             ))
    }
    return(sim2infDays(sim_pop, if_weights, n.days))
  })
}
