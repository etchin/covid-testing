require(tidyr)
require(dplyr)
library(fitdistrplus)
source("R/initialize_params.R")

load_test_sensitivity_raw <- function(sens_type){
  #False negative rate given day of infection
  #Modified from Kucirka et al. to remove Patient 13 from Danis et al.

  sens_by_day <- read.csv("data/sens_by_day_ci.csv", header = TRUE) %>%
    dplyr::select(fnr_med, fnr_lb, fnr_ub)
  
  if(sens_type == "lower") p_fn <- sens_by_day$fnr_ub
  if(sens_type == "median") p_fn <- sens_by_day$fnr_med
  if(sens_type == "upper") p_fn <- sens_by_day$fnr_lb
  if(sens_type == "perfect") p_fn <- rep(0,nrow(sens_by_day))
  if(sens_type == "random") p_fn <- apply(sens_by_day, 1, function(x) rtruncnorm(1, x[["fnr_lb"]], x[["fnr_ub"]], x[["fnr_med"]], 
                                                                                 (x[["fnr_ub"]] - x[["fnr_med"]])/2))
  x <- as.numeric(rownames(sens_by_day))
  y <- 1-p_fn
  y[1:5] <- rev(y[seq(5,18,3)])

  fit <- lm(y ~ bs(x,knots = c(1,3,7,10,15,18,20)))
  y1 <- predict(fit)
  fit_extrap <- approxExtrap(x,y,(max(x)+1):max_days_PCR_pos)
  
  sens <- c(y1, fit_extrap$y)
  sens[sens < 0] <- 0
  return(sens)
}

load_test_sensitivity <- function(sens_type){
  #False negative rate given day of infection
  #Modified from Kucirka et al. to remove Patient 13 from Danis et al.

  sens_by_day <- read.csv("data/sens_by_day_corrected.csv", header = TRUE) %>%
    dplyr::select(lower, median, upper)

  if(sens_type == "lower") sens <- sens_by_day$lower
  if(sens_type == "median") sens <- sens_by_day$median
  if(sens_type == "upper") sens <- sens_by_day$upper
  if(sens_type == "perfect") sens <- rep(1,nrow(sens_by_day))
  if(sens_type == "random") sens <- apply(sens_by_day, 1, function(x) rtruncnorm(1, x[["lower"]], x[["median"]], x[["upper"]],
                                                                                 (x[["upper"]] - x[["median"]])/2))
  return(sens)
}


test_pos <- function(sim_pop, t, tparams){
  with(tparams,{
    p_i <- ceiling(apply(sim_pop, 1, function(y) day_in_infection(y, t)))
    toffset <- days.incubation - 5 # can only test positive 5 days prior to symptom onset
    p_i <- p_i - toffset
    sTested_tp <- sapply(p_i, function(y) ifelse(y > 0 & y <= length(p.sens), rbinom(1,1,p.sens[y]), 0)) 
    sim_pop[sTested_tp == 1, DayOfDetection := t + n.delay]
    return(sim_pop)
  })
}

test_neg <- function(sim_pop, t, tparams){
  with(tparams,{
    sTested_fp <- rbinom(nrow(sim_pop), 1, 1-p.spec)
    sim_pop[sTested_fp == 1, DayOfDetection := t + n.delay]
    return(sim_pop)
  })
}

test_pop <- function(sim_pop, t, tparams){
  pids <- sim_pop[State %in% c(1,2) & DayOfDetection < DayOfInfection,ID]
  sim_pop[pids,] <- test_pos(sim_pop[pids,], t, tparams)
  sim_pop[DayOfDetection == t + n.delay & DayOfDetection_iwd == 0 & DayOfInfection > 0, 
          DayOfDetection_iwd := DayOfDetection]
  sim_pop[State == 0,] <- test_neg(sim_pop[State == 0,], t, tparams)
  rids <- sim_pop[State == 3 & DayOfDetection < DayOfInfection,ID]
  sim_pop[rids,] <- test_neg(sim_pop[rids,], t, tparams)
  return(sim_pop)
}
