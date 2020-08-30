states <- c("SU","IB","SX","R")
n.states <- length(states)

n_reps <- 5

max_days_PCR_pos <- 25

days.quarantine <- 14

n.contacts <- 10 # daily contacts

alpha.a <- 0.5
p.a <- 0.4
n_sims <- 100
p.community <- 5000/1e6  #proportion of community that is infected

hosp.delay <- 4 # days since symptom onset
icu.delay <- 4 # assume same as hospitalization
death.delay <- 10 # from symptom onset

p.a = 0.4 # proportion asymptomatic
alpha.a <- 0.5

p.community <- 0.005 # 0.5% of community infected
beta.t <- 0.3489 # Calculated using ngm
