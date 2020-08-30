Frequency of testing in the workplace

]The risk multiplier is defined as an increase between worker-community contacts. We calculate the risk multiplier necessary to bring the workplace R0 to 1.5, 2, and 2.5.

Rscript microsims.R <risk> <delay> <sensitivity> <sensitivity  multiplier>

risk.multiplier: additional risk for workers from community-worker contact
day: number of days for testing delay
sensitivity: "lower","upper","median" represent the estimates for RT-PCR test sensitivites by day of infection. Modified from the Kucirka et al. paper (see preprint for more details). "random" samples using a truncated normal distribution between the lower and upper confidence interval bounds from the Kucirka et al paper. "perfect" represents 100% sensitivity
sensitivity  multiplier: Discount factor to multiply sensitivity by. Assumed to be 0.8 for rapid PCR tests
