library(ggplot2)
library(dplyr)
library(tidyr)
source("R/utils.R")

# generate infectious profile
ni1 <- 12
ni2 <- 20
sx_onset <- 12.27248

infectiousness_df <- as.data.frame(matrix(nrow = 0, ncol = 3))

n_reps <- 1000
for(i in 1:n_reps){
  d <- infectious_fxn(ni1, ni2)
  tmp <- cbind(1:length(d), unlist(d), rep(i, length(d)))
  infectiousness_df <- rbind(infectiousness_df, tmp)
}

colnames(infectiousness_df) <- c("day","infectiousness","run")

idf <- infectiousness_df %>%
  group_by(day) %>%
  summarise(Mean = mean(infectiousness),
            Median = median(infectiousness),
            Q1 = quantile(infectiousness, .25),
            Q3 = quantile(infectiousness, .75))

ggplot(idf, aes(x = day, y = 100*Median, ymin = 100*Q1, ymax = 100*Q3)) + 
  geom_line() + geom_ribbon(fill = "red4",alpha = 0.5) +
  geom_vline(xintercept = sx_onset, linetype="dashed") +
  annotate("text",x = 19, y=15, label = "Symptom onset", size = 4.5) +
  theme_minimal(base_size = 20) + 
  xlab("Day of infection") + 
  ylab("Infectiousness (%)")

ggplot(infectiousness_df, aes(x = day, y = 100*infectiousness, group = factor(run))) + 
  geom_line(color = "red4", alpha = 0.01) + 
  geom_vline(xintercept = sx_onset, linetype="dashed") +
  annotate("text",x = 19, y=15, label = "Symptom onset", size = 4.5) +
  theme_minimal(base_size = 20) + 
  xlab("Day of infection") + 
  ylab("Infectiousness (%)")
