require(matrixStats)
require(readr)
require(tidyr)
require(dplyr)
require(stringr)
require(data.table)
require(tidyverse)

source("R/initialize_params.R")
source("R/utils.R")

span <- .75
n.people <- 10

get_inf_work_days <- function(x, inf_col){
  cnames <- c("DayOfInfection", "DaysExposed", "DaysEarlyInf", "DaysLateInf", "DayOfDetection",
              "InfWorkDays", "InfEarlyWorkDays", "InfLateWorkDays")
  x <- as.numeric(x[cnames])
  names(x) <- cnames
  
  day.early <- ceiling(x["DayOfInfection"] + x["DaysExposed"])
  day.late <- ceiling(x["DayOfInfection"] + x["DaysExposed"] + x["DaysEarlyInf"])
  day.end <- ceiling(x["DayOfInfection"] + x["DaysExposed"] + x["DaysEarlyInf"] + x["DaysLateInf"])

  if(x["DayOfDetection"] > 0) {
    day.end <- min(day.end, x["DayOfDetection"])
    if(day.end > day.late) day.late <- day.end
    }
  
  early_work_days <- sample(day.early:(day.late-1), as.integer(x["InfEarlyWorkDays"]))
  late_work_days <- c()
  if(day.end > day.late) late_work_days <- sample(day.late:day.end, as.integer(x["InfLateWorkDays"]))
  
  return(list("early" = early_work_days, "late" = late_work_days))
}

get_weighted_inf_work_days <- function(x, inf_col){
  cnames <- c("DayOfInfection", "InfectionType", "DaysExposed", "DaysEarlyInf", "DaysLateInf", "DayOfDetection",
              "InfWorkDays", "InfEarlyWorkDays", "InfLateWorkDays")
  x <- as.numeric(x[cnames])
  names(x) <- cnames
  
  day.early <- ceiling(x["DayOfInfection"] + x["DaysExposed"])
  day.late <- ceiling(x["DayOfInfection"] + x["DaysExposed"] + x["DaysEarlyInf"])
  day.end <- ceiling(x["DayOfInfection"] + x["DaysExposed"] + x["DaysEarlyInf"] + x["DaysLateInf"])
  
  if(x["DayOfDetection"] > 0) {
    day.end <- min(day.end, x["DayOfDetection"])
    if(day.end > day.late) day.late <- day.end
  }
  
  early_work_days <- sample(day.early:(day.late-1), as.integer(x["InfEarlyWorkDays"]))
  late_work_days <- c()
  if(day.end > day.late) late_work_days <- sample(day.late:day.end, as.integer(x["InfLateWorkDays"]))
  alpha <- ifelse(x["InfectionType"]==1, 1, alpha.a)
  return(list("early" = alpha*early_work_days, "late" = alpha*alpha.late*late_work_days))
}

calc_daily_outcomes <- function(output, n.days) {
  day_start <- 1
  day_end <- day_start + n.days - 1
  #all daily infections (by day of infection)
  daily_cases <- output %>% filter(State!=0) %>%
    group_by(DayOfInfection) %>% summarise(Cases=n()) %>% mutate(DayOfInfection=ifelse(DayOfInfection==0, 1, DayOfInfection))
  daily_cases <- rename(daily_cases, Day=DayOfInfection)
  #all infections work days (by day of infection)
  output_iwd <- output %>% filter(InfWorkDays!=0)
  inf_work_days <- list(list("early"=c(), "late"=c()))
  inf_weighted_work_days <- list(list("early"=c(), "late"=c()))
  if(nrow(output_iwd) > 0) inf_work_days <- apply(output_iwd, 1, get_inf_work_days)
  if(nrow(output_iwd) > 0) inf_weighted_work_days <- apply(output_iwd, 1, get_weighted_inf_work_days)
  early_iwd <- do.call(c, lapply(inf_work_days, function(x) x$early))
  late_iwd <- do.call(c, lapply(inf_work_days, function(x) x$late))
  early_weighted_iwd <- do.call(c, lapply(inf_weighted_work_days, function(x) x$early))
  late_weighted_iwd <- do.call(c, lapply(inf_weighted_work_days, function(x) x$late))
  if(length(early_iwd) == 0) early_iwd <- integer(0)
  if(length(late_iwd) == 0) late_iwd <- integer(0)
  if(length(early_weighted_iwd) == 0) early_weighted_iwd <- integer(0)
  if(length(late_weighted_iwd) == 0) late_weighted_iwd <- integer(0)
  daily_early <- data.frame(Day=early_iwd) %>%
    group_by(Day) %>% summarise(EarlyInfWorkers = n())
  daily_late <- data.frame(Day=late_iwd) %>%
    group_by(Day) %>% summarise(LateInfWorkers = n())
  daily_weighted_early <- data.frame(Day=early_weighted_iwd) %>%
    group_by(Day) %>% summarise(WeightedEarlyInfWorkers = n())
  daily_weighted_late <- data.frame(Day=late_weighted_iwd) %>%
    group_by(Day) %>% summarise(WeightedLateInfWorkers = n())
  output %>% filter(DayOfInfection!=0) %>%
    group_by(DayOfInfection) %>% summarise(Cases=n())
  #daily symptomatic cases (by day of symptom onset - only those with symptoms by end of model period)
  daily_shid <- output %>% 
    mutate(DayOfSymptoms=ifelse(InfectionType==1 & DayOfInfection>0, 
                                round(DayOfInfection + DaysExposed + DaysEarlyInf), 0),
           DayOfHospital=ifelse(Hospitalization==1 &
                                  DayOfInfection>0, 
                                round(DayOfInfection + DaysExposed + DaysEarlyInf) + hosp.delay, 0),
           DayOfICU=ifelse(ICU==1 &
                             DayOfInfection>0, 
                           round(DayOfInfection + DaysExposed + DaysEarlyInf) + icu.delay, 0),
           DayOfDeath=ifelse(Death==1 &
                               DayOfInfection>0, 
                             round(DayOfInfection + DaysExposed + DaysEarlyInf) + death.delay, 0))
  daily_symptom <- daily_shid %>% filter(DayOfInfection!=0 &
                                           DayOfSymptoms <= day_end) %>%
    group_by(DayOfSymptoms) %>% summarise(SymptomaticCases=n())
  daily_symptom <- rename(daily_symptom, Day=DayOfSymptoms)
  #daily detected cases (by day of detection)
  daily_obs_cases <- output %>% filter(DayOfDetection!=0) %>%
    group_by(DayOfDetection) %>% summarise(ObsCases=n())
  daily_obs_cases <- rename(daily_obs_cases, Day=DayOfDetection)
  #daily hospitalizations (by day of hospitalization)
  daily_hosp <- daily_shid %>% filter(DayOfInfection!=0 &
                                        DayOfHospital <= day_end) %>%
    group_by(DayOfHospital) %>% summarise(Hospitalizations=n())
  daily_hosp <- rename(daily_hosp, Day=DayOfHospital)
  
  #daily icu (by day of icu)
  daily_icu <- daily_shid %>% filter(DayOfInfection!=0 &
                                       DayOfICU <= day_end) %>%
    group_by(DayOfICU) %>% summarise(ICU=n())
  daily_icu <- rename(daily_icu, Day=DayOfICU)
  
  #daily deaths (by day of death - covid only)
  daily_death <- daily_shid %>% filter(DayOfInfection!=0 &
                                         DayOfDeath <= day_end) %>%
    group_by(DayOfDeath) %>% summarise(Deaths=n())
  daily_death <- rename(daily_death, Day=DayOfDeath)
  
  #calculate population at sim start
  Pop <- nrow(output)
  Total_eiwd <- sum(ceiling(output$DaysEarlyInf))
  Total_liwd <- sum(ceiling(output$DaysLateInf))
  Total_iwd <- sum(ceiling(output$DaysEarlyInf + output$DaysLateInf))
  output <- output %>% 
    mutate(WeightedDaysEarlyInf = ifelse(InfectionType==0,
                                         alpha.a*ceiling(DaysEarlyInf),
                                         ceiling(DaysEarlyInf)),
           WeightedDaysLateInf = alpha.late*ifelse(InfectionType==0,
                                         alpha.a*ceiling(DaysLateInf),
                                         ceiling(DaysLateInf)),
           WeightedDaysInf = ifelse(InfectionType==0,
                                    alpha.a*(ceiling(DaysEarlyInf) + alpha.late*floor(DaysLateInf)),
                                    ceiling(DaysEarlyInf) + alpha.late*floor(DaysLateInf))
           )
  Total_w_eiwd <- sum(output$WeightedDaysEarlyInf)
  Total_w_liwd <- sum(output$WeightedDaysLateInf)
  Total_w_iwd <- sum(output$WeightedDaysInf)
  days <- data.frame(Day=seq(day_start, day_end))
  daily_all <- left_join(days, daily_cases, by=c("Day"))
  daily_all <- left_join(daily_all, daily_early, by=c("Day"))
  daily_all <- left_join(daily_all, daily_late, by=c("Day"))
  daily_all <- left_join(daily_all, daily_weighted_early, by=c("Day"))
  daily_all <- left_join(daily_all, daily_weighted_late, by=c("Day"))
  daily_all <- left_join(daily_all, daily_symptom, by=c("Day"))
  daily_all <- left_join(daily_all, daily_obs_cases, by=c("Day"))
  daily_all <- left_join(daily_all, daily_hosp, by=c("Day"))
  daily_all <- left_join(daily_all, daily_icu, by=c("Day"))
  daily_all <- left_join(daily_all, daily_death, by=c("Day"))
  daily_all[is.na(daily_all)] <- 0 #NAs are truly 0s
  daily_all <- daily_all %>% 
    mutate(InfWorkers = EarlyInfWorkers + LateInfWorkers,
           WeightedInfWorkers = WeightedEarlyInfWorkers + WeightedLateInfWorkers) %>%
    arrange(Day) %>%
    mutate(CumCases=cumsum(Cases), 
           CumEarlyInfWorkers=cumsum(EarlyInfWorkers), 
           CumLateInfWorkers=cumsum(LateInfWorkers),
           CumInfWorkers=cumsum(InfWorkers),
           CumWeightedEarlyInfWorkers=cumsum(WeightedEarlyInfWorkers), 
           CumWeightedLateInfWorkers=cumsum(WeightedLateInfWorkers),
           CumWeightedInfWorkers=cumsum(WeightedInfWorkers),
           CumObsCases=cumsum(ObsCases), CumSymptomaticCases=cumsum(SymptomaticCases),
           CumHospitalizations=cumsum(Hospitalizations), CumICU=cumsum(ICU), CumDeaths=cumsum(Deaths))
  daily_all <- daily_all %>% mutate(CumInfProp=CumCases*100/Pop, 
                                    CumInfWorkersProp=CumInfWorkers*100/Total_iwd, 
                                    CumEarlyInfWorkersProp=CumEarlyInfWorkers*100/Total_eiwd, 
                                    CumLateInfWorkersProp=CumLateInfWorkers*100/Total_liwd, 
                                    CumWeightedInfWorkersProp=CumWeightedInfWorkers*100/Total_w_iwd, 
                                    CumWeightedEarlyInfWorkersProp=CumWeightedEarlyInfWorkers*100/Total_w_eiwd, 
                                    CumWeightedLateInfWorkersProp=CumWeightedLateInfWorkers*100/Total_w_liwd, 
                                    CumObsProp=CumObsCases*100/Pop,
                                    CumSymptomProp=CumSymptomaticCases*100/Pop,
                                    CumHospitalizedProp=CumHospitalizations*100/Pop,
                                    CumICUProp=CumICU*100/Pop,
                                    CumDeadProp=CumDeaths*100/Pop)
  return(daily_all)
}

convert_test_freq <- function(x){
  if(is.na(x)) return("No routine testing")
  if(x == 1) return("Daily")
  if(x == 3) return("Twice weekly")
  if(x == 7) return("Weekly")
  if(x == 14) return("Twice monthly")
  if(x == 30) return("Monthly")
  return(NA)
}

convert_test_freq2 <- function(x){
  if(x == "Daily") return(1)
  if(x == "Twice weekly") return(3)
  if(x == "Weekly") return(7)
  if(x == "Twice monthly") return(14)
  if(x == "Monthly") return(30)
  return(NA)
}

find_ci_norm <- function(x,y){
  f.t.test <- t.test(x, y)
  f.ci <- -f.t.test$conf.int/f.t.test$estimate[2]
  x.stats <- (f.t.test$estimate[2]-summary(x))/f.t.test$estimate[2]
  x.stats <- c(f.ci[1], f.ci[2], x.stats)
  x.stats <- 1-x.stats
  names(x.stats) <- c("lower","upper","min","q1","median","mean","q3","max")
  return(x.stats)
}

#CI around mean
get_CI_half_width <- function(x, prob) {
  n <- length(x)
  z_t <- qt(1 - (1 - prob) / 2, df = n - 1)
  z_t * sd(x) / sqrt(n)
}

get_binom_CI_half_width <- function(x, prob) {
  n <- length(x)
  z_t <- qt(1 - (1 - prob) / 2, df = n - 1)
  #z_t * sd(x) / sqrt(n)
  z_t * sqrt(mean(x) * (1-mean(x)) / length(x))
}

lower <- function(x, prob = 0.95) {
  # do some inference. if all < 1, then use binom. otherwise use regular
  if(all(x <= 1)) return(mean(x) - get_binom_CI_half_width(x, prob))
  return(mean(x) - get_CI_half_width(x, prob))
}

upper <- function(x, prob = 0.95) {
  if(all(x <= 1)) return(mean(x) + get_binom_CI_half_width(x, prob))
  return(mean(x) + get_CI_half_width(x, prob))
}

#plot model output
plot_point <- function(output, x_var, y_var, col_var, x_lab, y_lab, col_lab, title) {
  ggplot(output, aes_string(x=x_var, y=y_var, color=col_var)) +
    geom_point() +
    labs(x=x_lab, y=y_lab, color=col_lab) +
    theme_bw() + ggtitle(title)
}
plot_line <- function(output, x_var, y_var, col_var, x_lab, y_lab, col_lab, title, smooth) {
  if(smooth){
    col_names <- c("x","y","color")
    df <- data.frame(matrix(nrow = 0, ncol = length(col_names), dimnames = list(c(), col_names)))
    for(c in unlist(unique(output[,col_var]))){
      output_c <- output %>% filter(get(col_var) == c)
      x <- seq(min(output_c[, x_var]), max(output_c[, x_var]), by=1)
      y_fit <- loess(as.formula(paste0(y_var, " ~ ", x_var)), output_c, span = span)
      tmp <- data.frame(matrix(nrow = length(x), ncol = length(col_names), dimnames = list(c(), col_names)))
      tmp$x <- x
      tmp$y <- predict(y_fit, x)
      tmp$color <- c
      df <- rbind(df, tmp)
    }
    df[df$y < 0, "y"] <- 0
    col_levels <- output %>% pull((col_var)) %>% levels()
    if(!is.null(col_levels)) df$color <- factor(df$color, col_levels)
    if(length(unique(df$color))==1){
      p <- ggplot(df, aes(x=x, y=y, color=color, group = 1))
    } else{
      p <- ggplot(df, aes(x=x, y=y, color=color))
    }
    p <- p +
      geom_line(size = 1.5, alpha = 0.75) +
      labs(x=x_lab, y=y_lab, color=col_lab) +
      theme_bw() + ggtitle(title) +
      coord_cartesian(ylim = c(0, NA))
  } else{
    if(length(unique(output[,col_var]))==1){
      p <- ggplot(output, aes_string(x=x_var, y=y_var, color=col_var, group = 1))
    } else{
      p <- ggplot(output, aes_string(x=x_var, y=y_var, color=col_var))
    }
    p <- p +
      geom_line(size = 1.5, alpha = 0.75) +
      labs(x=x_lab, y=y_lab, color=col_lab) +
      theme_bw() + ggtitle(title) +
      coord_cartesian(ylim = c(0, NA))
  }
  p
}
plot_grp_line <- function(output, x_var, y_var, col_var, grp_var, x_lab, y_lab, col_lab, title, smooth) {
  if(smooth){
    col_names <- c("x","y","color")
    df <- data.frame(matrix(nrow = 0, ncol = length(col_names), dimnames = list(c(), col_names)))
    for(c in unlist(unique(output[,grp_var]))){
      output_c <- output %>% filter(get(grp_var) == c)
      x <- seq(min(output_c[, x_var]), max(output_c[, x_var]), by=1)
      y_fit <- loess(as.formula(paste0(y_var, " ~ ", x_var)), output_c, span = span)
      tmp <- data.frame(matrix(nrow = length(x), ncol = length(col_names), dimnames = list(c(), col_names)))
      tmp$x <- x
      tmp$y <- predict(y_fit, x)
      tmp$group <- c
      tmp$color <- unique(output %>% filter(!!sym(grp_var)==c) %>% pull(!!sym(col_var)))
      df <- rbind(df, tmp)
    }
    df[df$y < 0, "y"] <- 0
    p <- ggplot(df, aes(x=x, y=y, color=color, group=group)) +
      geom_line(size = 1.5, alpha = 0.75) +
      labs(x=x_lab, y=y_lab, color=col_lab) +
      theme_bw() + ggtitle(title) +
      coord_cartesian(ylim = c(0, NA))
  } else{
    p <- ggplot(output, aes_string(x=x_var, y=y_var, color=col_var, group=grp_var)) +
      geom_line(size = 1.5, alpha = 0.75) +
      labs(x=x_lab, y=y_lab, color=col_lab) +
      theme_bw() + ggtitle(title) +
      coord_cartesian(ylim = c(0, NA))
  }
  p
}

plot_line_CIs <- function(output, x_var, y_var, y_lb_var, y_ub_var, col_var, x_lab, y_lab, col_lab, title, smooth) {
  if(smooth){
    col_names <- c("x","y","y_lb","y_ub","color")
    df <- data.frame(matrix(nrow = 0, ncol = length(col_names), dimnames = list(c(), col_names)))
    for(c in unlist(unique(output[,col_var]))){
      output_c <- output %>% filter(get(col_var) == c)
      x <- seq(min(output_c[, x_var]), max(output_c[, x_var]), by=1)
      y_fit <- loess(as.formula(paste0(y_var, " ~ ", x_var)), output_c, span = span)
      y_lb_fit <- loess(as.formula(paste0(y_lb_var, " ~ ", x_var)), output_c, span = span)
      y_ub_fit <- loess(as.formula(paste0(y_ub_var, " ~ ", x_var)), output_c, span = span)
      tmp <- data.frame(matrix(nrow = length(x), ncol = length(col_names), dimnames = list(c(), col_names)))
      tmp$x <- x
      tmp$y <- predict(y_fit, x)
      tmp$y_lb <- predict(y_lb_fit, x)
      tmp$y_ub <- predict(y_ub_fit, x)
      tmp$color <- c
      df <- rbind(df, tmp)
    }
    df[df$y < 0, "y"] <- 0
    df[df$y_lb < 0, "y_lb"] <- 0
    df[df$y_ub < 0, "y_ub"] <- 0
    col_levels <- output %>% pull((col_var)) %>% levels()
    if(!is.null(col_levels)) df$color <- factor(df$color, col_levels)
    if(length(unique(df$color)) == 1){
      p <- ggplot(df, aes(x=x, y=y, ymin=y_lb, ymax=y_ub, color=color, fill=color, group = 1))
    } else{
      p <- ggplot(df, aes(x=x, y=y, ymin=y_lb, ymax=y_ub, color=color, fill=color))
    }
    p <- p +
      geom_line(size = 1.5, alpha = 0.75) +
      geom_ribbon(alpha=0.5) +
      labs(x=x_lab, y=y_lab, color=col_lab, fill=col_lab) +
      theme_bw() + ggtitle(title) +
      coord_cartesian(ylim = c(0, NA))
  } else{
    if(length(unique(output[,col_var])) == 1){
      p <- ggplot(output, aes_string(x=x_var, y=y_var, ymin = y_lb_var, ymax=y_ub_var, color=col_var, fill=col_var, group = 1))
    } else{
      p <- ggplot(output, aes_string(x=x_var, y=y_var, ymin = y_lb_var, ymax=y_ub_var, color=col_var, fill=col_var))
    }
    p <- p +
      geom_line(size = 1.5, alpha = 0.75) + 
      geom_ribbon(alpha=0.5) +
      labs(x=x_lab, y=y_lab, color=col_lab, fill=col_lab) +
      theme_bw() + ggtitle(title) +
      coord_cartesian(ylim = c(0, NA))
  }
  p
}
