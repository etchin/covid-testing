library(stringr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)

tfArray <- c(1:5,7,seq(10,30,5), NA)

find_ci_norm <- function(x,y){
  f.t.test <- t.test(x, y)
  f.ci <- -f.t.test$conf.int/f.t.test$estimate[2]
  x.stats <- (f.t.test$estimate[2]-summary(x))/f.t.test$estimate[2]
  x.stats <- c(f.ci[1], mean(f.ci), f.ci[2], x.stats)
  names(x.stats)[1:3] <- c("lower.ci","mean.ci","upper.ci")
  return(x.stats)
}

find_reduction <- function(f1,f2){
  f3 <- as.data.frame(t(sapply(1:nrow(f1), function(i) find_ci_norm(f1[i,], f2[i,]))))
  f3$population <- rownames(f1)
  return(f3)
}

find_infections_days <- function(fis, alpha.a, alpha.late){
  infs <- lapply(fis, function(f) read.csv(f, header = TRUE, row.names = 1))
  infs <- rbindlist(infs)
  inf_array <- list()
  for(i in 1:length(tfArray)){
    tf_infs <- infs[,(3*(i-1)+1):(3*i)]
    tf_infs[,1] <- alpha.a*tf_infs[,1]
    tf_infs[,2] <- alpha.a*alpha.late*tf_infs[,2]
    tf_infs <- cbind(tf_infs, rowSums(tf_infs))
    colnames(tf_infs) <- c("asx_early","asx_late","sx","total")
    inf_array[[as.character(tfArray[i])]] <- t(tf_infs)
  }
  red_array <- list()
  for(tf in head(tfArray,-1)){
    tmp <- find_reduction(inf_array[[as.character(tf)]], inf_array[[as.character(tail(tfArray,1))]])
    tmp$tf <- tf
    red_array[[as.character(tf)]] <- tmp
  }
  red_array <- rbindlist(red_array)
  return(red_array)
}


fi_dir <- "data/new_sims/"
sim_sum_fis <- list.files(fi_dir, pattern="*.csv")

param_combos <- list(c(1.625, 1, "random",1.8), c(2.54, 1, "random",1.8), c(3.445, 1, "random",1.8),
                  c(1.625, 0, "perfect",1.8)
                  )

inf_red <- list()
all_files <- list.files(fi_dir,pattern = "*.csv")
counter <- 1
for(p in param_combos){
  fi_pattern <- paste("^sims",p[1],p[2],0.5,p[3],p[4],sep = "_")
  fis <- all_files[grep(fi_pattern,all_files)]
  f_counter_factual <- find_infections_days(paste0(fi_dir,fis), 0.4, 0.5)
  f_counter_factual$risk <- p[1]
  f_counter_factual$delay <- p[2]
  f_counter_factual$sensitivity <- p[3]
  f_counter_factual$R0 <- p[4]
  inf_red[[counter]] <- f_counter_factual
  counter <- counter + 1
}

inf_red <- rbindlist(inf_red)

convert_population <- function(x){
  if(x == "asx_early") return("Sub-clinical, early")
  if(x == "asx_late") return("Sub-clinical, late")
  if(x == "sx") return("Clinical")
  return('Total')
}

convert_risk_to_R0 <- function(x){
  if(x == 1.625) return(1.5)
  if(x == 2.54) return(2)
  if(x == 3.445) return(2.5)
}

inf_red <- inf_red %>% 
  dplyr::mutate(population = sapply(population, convert_population),
         risk = as.numeric(risk),
         R0 = sapply(risk, convert_risk_to_R0)
  )

write.csv(inf_red, "reduction.csv")

red_ci <- inf_red %>% 
  filter(population == "Total") %>% 
  mutate(CI = paste0(round(mean.ci*100, 1), "% (", 
                     round(upper.ci*100,1), 
                     ",", round(lower.ci*100,1), ")")) %>% 
  select(risk, delay, sensitivity, R0, tf, CI)

write.csv(red_ci, "confidence_intervals_all.csv")

convert_tf <- function(tf){
  if(tf == 1) return("Daily")
  if(tf == 30) return("Monthly")
  if(tf == 1000) return("No testing")
  return(paste0(tf, " days"))
}

tf_levels <- c("Daily",paste0(c(2:5,7,seq(10,25,5)), " days"), "Monthly")



estimatedR <- inf_red %>%
  dplyr::mutate(Re = R0*(1-mean.ci),
         Re.lower = R0*(1-upper.ci),
         Re.upper = R0*(1-lower.ci),
         Re.Q1 = R0*(1-`1st Qu.`),
         Re.Q3 = R0*(1-`3rd Qu.`)
         )

write.csv(estimatedR, "effective_R.csv")

ggplot(inf_red %>%
         filter(tf < 1000, population != "Total", R0 == 2.5, 
                sensitivity == 'random', delay == 1) %>%
         mutate("Reduction in infectious working days" = 100*mean.ci,
                ymin = sapply(100*`3rd Qu.`, function(x) min(100, max(0,x))),
                ymax = sapply(100*`1st Qu.`, function(x) min(100,max(0,x))), 
                "Testing Frequency" = tf,
                Infection = factor(population, 
                                   levels = c("Total","Clinical","Sub-clinical, early", "Sub-clinical, late")))
       , aes(x = `Testing Frequency`, y = `Reduction in infectious working days`, 
             group = Infection, fill = Infection)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.75) +
  geom_line() +
  theme_minimal(base_size = 14) +
  ylab("% Reduction in infectious working days") +
  xlab("Testing frequency (days)") +
  ggtitle("Percent reduction in infectious working days \nunder different testing frequencies") +
  ylim(0,100)

ggsave("figures/reduction_days_highrisk.pdf")

ggplot(estimatedR %>%
         filter(tf < 1000, population == "Total", R0 == 2.5, 
                sensitivity == 'random', delay == 1) %>%
         mutate("Reduction in infectious working days" = 100*mean.ci,
                ymin = sapply(100*`3rd Qu.`, function(x) min(100, max(0,x))),
                ymax = sapply(100*`1st Qu.`, function(x) min(100,max(0,x))), 
                "Testing Frequency" = tf)
       , aes(x = `Testing Frequency`, y = `Reduction in infectious working days`)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 1, fill = "gray") +
  geom_line() +
  theme_minimal(base_size = 20) +
  ylab("% Reduction in infectiousness") +
  xlab("Testing frequency (days)") +
  ylim(0,100)

ggsave("figures/reduction_infectiousness_highrisk_iqr.pdf")


ggplot(estimatedR %>%
         filter(delay == 1, sensitivity == "random", population == "Total"), 
       aes(x = tf, y = Re, group = R0, fill = R0)) +
  geom_ribbon(aes(ymin = Re.Q1, ymax = Re.Q3), alpha= 0.75) +
  geom_line() +
  geom_hline(yintercept = 1, size = 1.5, linetype = "dashed") + 
  scale_fill_gradient(low = "#355C7D",
                      high = "#F67280",
                      guide = "legend") + 
  theme_minimal(base_size = 20) +
  xlab("Testing frequency (days)") +
  ylab("Effective Reproduction Number")
ggsave("figures/estimated_Re.pdf")


