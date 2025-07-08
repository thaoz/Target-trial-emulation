# Obtain confidence intervals for the difference in survival probabilities
# using bootstraping method
# author: Thao Le

library(tidyverse)
library(here)
library(timereg)
library(survival)
library(furrr)
library(progressr)

source("AuxFun.R")

# FUNCTIONS
seq.fun <- function(dat){
  
  
  # 1. generate trials data ######################################################
  
  # define eligibility: a patient is eligible to be included in a trial at time t if they have not had any transfusion prior to t. 
  dat <- dat %>% 
    group_by(newid) %>% 
    mutate(eli = f.eligible(A))
  
  
  # define eligible id
  eli0 <- dat %>% filter(time == 0 & eli == 1) %>% distinct(newid)
  eli1 <- dat %>% filter(time == 1 & eli == 1) %>% distinct(newid)
  eli2 <- dat %>% filter(time == 2 & eli == 1) %>% distinct(newid)
  eli3 <- dat %>% filter(time == 3 & eli == 1) %>% distinct(newid)
  eli4 <- dat %>% filter(time == 4 & eli == 1) %>% distinct(newid)
  eli5 <- dat %>% filter(time == 5 & eli == 1) %>% distinct(newid)
  eli6 <- dat %>% filter(time == 6 & eli == 1) %>% distinct(newid)
  
  
  # trial data
  # censor a patient at time t if they start a new transfusion at a different threshold. 
  # in case of no transfusion, continue follow-up
  
  trial0 <- seq.trial.data.gen_boot_2(id.list = eli0$newid, time.s = 0, data = dat)
  trial1 <- seq.trial.data.gen_boot_2(id.list = eli1$newid, time.s = 1, data = dat)
  trial2 <- seq.trial.data.gen_boot_2(id.list = eli2$newid, time.s = 2, data = dat)
  trial3 <- seq.trial.data.gen_boot_2(id.list = eli3$newid, time.s = 3, data = dat)
  trial4 <- seq.trial.data.gen_boot_2(id.list = eli4$newid, time.s = 4, data = dat)
  trial5 <- seq.trial.data.gen_boot_2(id.list = eli5$newid, time.s = 5, data = dat) # no A0s = {1,2}
  trial6 <- seq.trial.data.gen_boot_2(id.list = eli6$newid, time.s = 6, data = dat) # no A0s = {1,2}
  
  
  # stack data; only use data from 0 to 4
  mega.d <- rbind(trial0, trial1, trial2, trial3, trial4, trial5, trial6) 
  
  # 2. Calculating weights #######################################################
  # estimate weight separately for those switching from no transfusion --> transfusion
  # and from transfusion --> no transfusion
  
  # denominator: cens ~ all time-fixed cov + current values of time-varying covs
  den.model.01 <- glm(Cens ~ age + gender + ECPR + rrt_pre + Surgical + diag.group +
                        L.lead.rrt + L.lead.anyms + L.lead.barc,
                      data = subset(mega.d, A0.s == 0),
                      family = "binomial")
  den.model.10 <- glm(Cens ~ age + gender + ECPR + rrt_pre + Surgical + diag.group +
                        L.lead.rrt + L.lead.anyms + L.lead.barc,
                      data = subset(mega.d, A0.s != 0),
                      family = "binomial")
  mega.d$wt.den = 1
  mega.d$wt.den[mega.d$A0.s == 0] <- predict(den.model.01,
                                             newdata = subset(mega.d, A0.s == 0),
                                             type = "response")
  mega.d$wt.den[mega.d$A0.s != 0] <- predict(den.model.10,
                                             newdata = subset(mega.d, A0.s != 0),
                                             type = "response")
  
  # numerator: cens ~ all time-fixed cov + time-varying covs at t = s
  ## Use rownum: 1 to 7
  num.model.01 <- glm(Cens ~ age + gender + ECPR + rrt_pre + Surgical + diag.group +
                        as.factor(rownum)*(L0s.rrt + L0s.anyMS + L0s.barc),
                      data = subset(mega.d, A0.s == 0),
                      family = "binomial")
  num.model.10 <- glm(Cens ~ age + gender + ECPR + rrt_pre + Surgical + diag.group +
                        as.factor(rownum)*(L0s.rrt + L0s.anyMS + L0s.barc),
                      data = subset(mega.d, A0.s != 0),
                      family = "binomial")
  mega.d$wt.num = 1
  mega.d$wt.num[mega.d$rownum != 7 & mega.d$A0.s == 0] <- predict(num.model.01, 
                                                                  newdata = mega.d[mega.d$rownum != 7 & mega.d$A0.s == 0,], 
                                                                  type = "response")
  mega.d$wt.num[mega.d$rownum != 7 & mega.d$A0.s != 0] <- predict(num.model.10, 
                                                                  newdata = mega.d[mega.d$rownum != 7 & mega.d$A0.s != 0,], 
                                                                  type = "response")
  
  # ipw, and cumulative weight
  mega.d <- mega.d %>% 
    group_by(trial, newid) %>% 
    mutate(ipw.stab = wt.num/wt.den,
           .cumprod.wt = cumprod(ipw.stab),
           ipw.stab.cum = lag(.cumprod.wt, default = 1)) %>% 
    ungroup() %>% 
    mutate(
      pc95 = quantile(ipw.stab.cum, probs = 0.95, na.rm = T),
      pc99 = quantile(ipw.stab.cum, probs = 0.99, na.rm = T),
      
      ipw.stab.cum.trunc.95 = ifelse(ipw.stab.cum> pc95, pc95, ipw.stab.cum),
      ipw.stab.cum.trunc.99 = ifelse(ipw.stab.cum> pc99, pc99, ipw.stab.cum)
      
    )
  
  # 3. Obtain estimated survival probabilities ##################################
  # Fit survival model for weighted data
  mega.d$A0.s <- as.factor(mega.d$A0.s)
  
  # prediction time
  t.hor = seq(0,6,0.01)
  d.time0 <- dat %>% filter(time == 0)
  
  # matrix of time-fixed covaraites and L at time 0
  Z <- model.matrix(~ age + gender + ECPR + rrt_pre + Surgical + diag.group +
                      RRT + anyMS + SevereBARC,
                    data = d.time0)
  
  # orignal weight
  orginal.weight.results <- surv.aalen.func(data = mega.d, 
                                            weight = mega.d$ipw.stab.cum, 
                                            t.hor = t.hor, 
                                            Z = Z)
  trunc95.weight.results <- surv.aalen.func(data = mega.d, 
                                            weight = mega.d$ipw.stab.cum.trunc.95, 
                                            t.hor = t.hor, 
                                            Z = Z)
  
  out <- data.frame(surv0.org = orginal.weight.results$surv0,
                    surv1.org = orginal.weight.results$surv1,
                    surv2.org = orginal.weight.results$surv2,
                    surv3.org = orginal.weight.results$surv3,
                    surv0.t95 = trunc95.weight.results$surv0,
                    surv1.t95 = trunc95.weight.results$surv1,
                    surv2.t95 = trunc95.weight.results$surv2,
                    surv3.t95 = trunc95.weight.results$surv3,
                    time = t.hor
  )
  return(out)
}
run1.boot <- function(x){
  # draw a random sample & create new unique id
  boot.d <- d2[index[,x],] %>% 
    group_by(id) %>% 
    mutate(nid = row_number()) %>% 
    mutate(newid = paste0(id,".",nid))
  
  # create long data format
  boot.dlong <- boot.d %>%
    pivot_longer(-c(id:Surgical, dead, diag.group, newid),  
                 cols_vary = "slowest",
                 names_to = c("day",".value"),
                 names_pattern = "day(.+)_(.+)"
    ) %>% 
    filter(OnECMO == 1) %>% 
    arrange(newid) %>% 
    mutate(time = as.numeric(day) - 1,
           A = ifelse(is.na(Hbtrigger), 0,
                      ifelse(Hbtrigger <=70, 1,
                             ifelse(Hbtrigger >=90, 2, 3)))
    ) %>%
    mutate(time.stop = time + 1) %>% 
    group_by(newid) %>% 
    mutate(event = ifelse(dead == 0, 0,
                          ifelse(dead == 1 & row_number() == n(), 1, 0)),
           time.stop = ifelse(event == 1, time + runif(1, min = 0, max = 0.7), time.stop))
  
  # impute missing
  boot.dlong <- boot.dlong %>% 
    group_by(newid) %>% 
    fill(RRT) %>% 
    mutate(RRT = ifelse(is.na(RRT), rrt_pre, RRT))
  boot.dlong <- boot.dlong %>% 
    group_by(newid) %>% 
    fill(SevereBARC)
  
  # perform sequential trials + MSM analysis
  r <- seq.fun(dat = boot.dlong)
  return(r)
} 
my_pkg_fn <- function(B) {
  alongB = 1:B
  p <- progressor(steps = length(alongB))
  
  future_map(.x = 1:B,
             .f = run1.boot,
             .options = furrr_options(seed = 1), 
             .progress = T)
}
# RUN bootstraping
B = 1000
Npat = n_distinct(d1.long$id)
index = sample(Npat, size = Npat*B, replace = T)
dim(index) = c(Npat, B)


plan(multisession)
system.time(tmp <- my_pkg_fn(B = B))

# result
bootstrap.results <- tmp %>% 
  bind_rows(.id = "boot.id")

save("bootstrap_result.Rdata")
