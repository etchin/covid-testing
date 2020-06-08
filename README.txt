Frequency of testing in the workplace

To obtain a total R0 of 2.4, we assume that the R0 for the early clinical stage is 1.8.

The risk multiplier is defined as an increase between worker-community contacts. We calculate the risk multiplier necessary to bring the workplace R0 to 1.5, 2, and 2.5.

Rscript microsims.R <n_sims> <risk.multiplier> <testingDelay> <alpha.late> <sensitivity type> <R0 of early clinical> <replicate>

n_sims: number of simulations
risk.multiplier: additional risk for workers from community-worker contact
testingDelay: number of days for testing delay
alpha.late: discount factor for reduction in infectiousness for the late infectious stage
sensitivity type: "lower","upper","median" represent the estimates from the Kucirka et al. paper for RT-PCR test sensitivites by day of infection. "random" samples using a truncated normal distribution between the lower and upper confidence interval bounds from the Kucirka et al paper. "perfect" represents 100% sensitivity
R0 of early clinical: R0 for early clinical infectious stage. Assumed to be 1.8 in the paper
replicate: index of replicate simulation job
