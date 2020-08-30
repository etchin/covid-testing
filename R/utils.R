require(truncnorm)
require(Hmisc)
require(splines)
source ("R/initialize_params.R")


# incubation parameters from
# Li, Q. et al. Early transmission dynamics in Wuhan, China, of novel coronavirus-infected pneumonia. N. Engl. J. Med. 382, 1199–1207 (2020).
# 5.2 days (95% CI, 4.1–6.4 days) fit to a log normal distribution. rounded to nearest day
incubation_fxn <- function(n) return(round(rlnorm(n, 1.434065, 0.6612))) # log normal parameters for days of incubation

# Number of infectious days after symptom onset
# Since we vary infectiousness by day, we set them to 20 days
symptomatic_fxn <- function(n) return(rep(20, n))

infectious_fxn <- function(ni1, ni2){
  mso <- rtruncnorm(1, mean = 12.27248, a = 9, b = 15, sd = 2)
  d <- dgamma(seq(-(ni1-1),ni2) + mso, 20.51651, 1.592124)
  d <- d / sum(d)
  names(d) <- c(paste0("e", 1:ni1), paste0("l", 1:ni2))
  return(d)
}

calc_infectiousness <- function(iwd, if_weights)sum(c(0,if_weights)[1:(iwd+1)])

sim2infDays <- function(sim_pop, if_weights, n.days){
  #get number of infectious and total number
  sim_pop_i <- sim_pop[DayOfInfection > 0]
  sim_pop_i[,InfWorkDays := ifelse(DayOfDetection_iwd >= DayOfInfection, 
                                 DayOfDetection_iwd - DayOfInfection,
                                 pmin(DaysIncubation + DaysSymptomatic + DayOfInfection, n.days) - DayOfInfection)]
  sim_pop_i[,InfWorkDays := pmin(InfWorkDays, DaysIncubation + DaysSymptomatic)]
  sim_pop_i$infectiousness <- sapply(sim_pop_i$InfWorkDays, function(x) calc_infectiousness(x, if_weights))
  sim_pop_i[, winfectiousness := ifelse(InfectionType == 1, infectiousness, alpha.a*infectiousness)]
  return(c(sum(sim_pop_i[,winfectiousness]), nrow(sim_pop_i)))
}

update_state <- function(x, t){
  cols <- c("State","DayOfInfection","DaysIncubation","DaysSymptomatic")
  x <- as.numeric(x[cols])
  names(x) <- cols
  if(x["State"]==1 & t >= (x["DaysIncubation"] + x["DayOfInfection"])) return(2) # Late Infectious period
  if(x["State"]==2 & t >= (x["DaysIncubation"] + x["DaysSymptomatic"] + x["DayOfInfection"])) return(3) # Recovered
  return(x["State"])
}

day_in_infection <- function(x,t){
  cols <- c("State","DayOfInfection")
  x <- as.numeric(x[cols])
  names(x) <- cols
  if(x["State"] %in% c(1,2)) return(t - x["DayOfInfection"] + 1)
  return(NA)
}



