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
  x.stats <- c(f.ci[1], f.ci[2], x.stats)
  names(x.stats)[1:2] <- c("lower.ci","upper.ci")
  return(x.stats)
}

expand_sim <- function(v){
  df <- read_csv(v, skip=1, 
           progress = FALSE)[, -1]
  idf <- df %>%
    dplyr::select(starts_with("Infectiousness")) %>%
             pivot_longer(starts_with("Infectiousness"), names_prefix = "Infectiousness_", 
                      names_to = "Frequency", values_to = "Infectiousness")
  ndf <- df %>%
    dplyr::select(starts_with("nInfected")) %>%
    pivot_longer(starts_with("nInfected"), names_prefix = "nInfected_", 
                 names_to = "Frequency", values_to = "nInfected")
  df <- idf %>%
    left_join(ndf)
  return(df)
}

fi_dir <- "output/"
v_files<- Sys.glob("output/*.csv") #add file specs here

sims_all <- lapply(1:length(v_files), function(x) expand_sim(v_files[x]))
sims_all <- bind_rows(sims_all, .id="ModelRun")
sims_all$ModelRun <- as.numeric(sims_all$ModelRun)
  
#extract parameter values too so that we can avg over runs only if they have the same param value
params <- lapply(v_files, function(x) read_csv(x, n_max=1, col_names = FALSE, 
                                           col_types = cols(),
                                           progress = FALSE))
params <- bind_rows(params, .id="ModelRun")
col_names <- distinct(params %>% dplyr::select(-ModelRun) %>% mutate_all(~sub("\\:.*", "",.)) %>% mutate_all(~sub("^#", "",.)))
params <- params %>% mutate_all(~sub(".*:", "", .))
colnames(params) <- c("ModelRun", col_names)
params$ModelRun <- as.numeric(params$ModelRun)

convert_risk_to_R0 <- function(x){
  if(x == "0.645") return(1.5)
  if(x == "0.96") return(2)
  if(x == "1.28") return(2.5)
}

sims_all <- sims_all %>%
  left_join(params, by = "ModelRun") %>%
  mutate(R0 = sapply(risk, convert_risk_to_R0),
         Frequency = as.numeric(Frequency))

test.freqs <- unique(sims_all$Frequency)
cr_cols <- c("test_freq","outcome",tail(colnames(params), -1),
             "Lower","Upper","Min","Q1","Median","Mean","Q3","Max")
cum_red <- as.data.frame(matrix(NA, nrow=0, ncol=length(cr_cols)))
colnames(cum_red) <- cr_cols
ti <- 1
param_sets <- params %>% dplyr::select(-ModelRun) %>% unique()
for(r in 1:nrow(param_sets)){
  sims_params <- sims_all %>%
    right_join(param_sets[r,])
  for(tf in head(test.freqs,-1)){
    for(v in c("Infectiousness","nInfected")){
      x1 <- sims_params %>%
        filter(Frequency == tf) %>%
        pull(v)
      x2 <- sims_params %>%
        filter(is.na(Frequency)) %>%
        pull(v)
      y <- find_ci_norm(x1, x2)*100
      y[y < 0] <- 0
      cum_red[ti,] <- c(tf, v, unlist(param_sets[r,]), y)
      ti <- ti + 1
    }
  }
}

cum_red$R0 <- sapply(cum_red$risk, convert_risk_to_R0)
tf_levels <- c("Daily",paste0(c(2:5,7,seq(10,25,5)), " days"), "Monthly")

write.csv(cum_red, "reduction.csv")
cum_red <- read_csv("reduction.csv")[,-1]

red_ci <- cum_red %>%
  filter(outcome == "Infectiousness") %>%
  mutate(CI = paste0(round(Mean, 1), "% (", 
                     round(Upper,1), 
                     ",", round(Lower,1), ")")) %>% 
  dplyr::select(test_freq, R0, delay, sensitivity, sens_multiplier, CI)

write.csv(red_ci, "confidence_intervals_all.csv")

estimatedR <- cum_red %>%
  filter(outcome == "Infectiousness") %>%
  dplyr::mutate(R0_base = R0,
                R0 = R0_base*(100-Mean)/100,
                R0.upper = R0_base*(100-Upper)/100,
                R0.lower = R0_base*(100-Lower)/100,
                R0.Q3 = R0_base*(100-Q1)/100,
                R0.Q1 = R0_base*(100-Q3)/100,
                R0.Min = R0_base*(100-Max)/100,
                R0.Max = R0_base*(100-Min)/100
  ) %>%
  dplyr::select(test_freq, delay, sensitivity, sens_multiplier, starts_with("R0"))

write.csv(estimatedR, "R0.csv")
ggplot(estimatedR %>%
         filter(delay == 1, sensitivity == "random", sens_multiplier == 1), 
       aes(x = test_freq, y = R0, group = R0_base, fill = R0_base)) +
  geom_ribbon(aes(ymin = R0.Q1, ymax = R0.Q3), alpha= 0.75) +
  geom_line() +
  geom_hline(yintercept = 1, size = 1.5, linetype = "dashed") + 
  scale_fill_gradient(low = "#355C7D",
                      high = "#F67280",
                      guide = "legend") + 
  theme_minimal(base_size = 20) +
  xlab("Testing frequency (days)") +
  ylab("R_0") +
  scale_x_continuous(breaks = c(0,10,20,30), labels = c("Daily","10","20","Monthly"))
ggsave("figures/estimated_Re.pdf")


