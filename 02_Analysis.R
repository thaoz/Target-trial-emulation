# R code for the manuscript: 
# "Liberal or restrictive transfusion for veno-arterial Extracorporeal Membrane 
# Oxygenation patients: a target trial emulation using the OBLEX study data", Thao Le et.al


library(tidyverse)
library(here)
library(timereg)
library(survival)

# load aux functions
source("AuxFun.R")

# 1. generate trials data ######################################################
# d1.long: long format of the OBLEX data set. This data set contains:
# - id: patient ID
# - baseline and time-varying confounders
# - A: transfusion strategy
# - event, time, time.stop: event status, time start and time stop 

# define eligibility: a patient is eligible to be included in a trial at time t if they have not had any transfusion prior to t. 
d1.long.complete <- d1.long %>% 
  group_by(id) %>% 
  mutate(eli = f.eligible(A))


# define eligible id
eli0 <- d1.long.complete %>% filter(time == 0 & eli == 1) %>% distinct(id)
eli1 <- d1.long.complete %>% filter(time == 1 & eli == 1) %>% distinct(id)
eli2 <- d1.long.complete %>% filter(time == 2 & eli == 1) %>% distinct(id)
eli3 <- d1.long.complete %>% filter(time == 3 & eli == 1) %>% distinct(id)
eli4 <- d1.long.complete %>% filter(time == 4 & eli == 1) %>% distinct(id)
eli5 <- d1.long.complete %>% filter(time == 5 & eli == 1) %>% distinct(id)
eli6 <- d1.long.complete %>% filter(time == 6 & eli == 1) %>% distinct(id)


# trial data
# censor a patient at time t if they start a new transfusion at a different threshold. 
# in case of no transfusion, continue follow-up

trial0 <- seq.trial.data.gen2(id.list = eli0$id, time.s = 0, data = d1.long.complete)
trial1 <- seq.trial.data.gen2(id.list = eli1$id, time.s = 1, data = d1.long.complete)
trial2 <- seq.trial.data.gen2(id.list = eli2$id, time.s = 2, data = d1.long.complete)
trial3 <- seq.trial.data.gen2(id.list = eli3$id, time.s = 3, data = d1.long.complete)
trial4 <- seq.trial.data.gen2(id.list = eli4$id, time.s = 4, data = d1.long.complete)
trial5 <- seq.trial.data.gen2(id.list = eli5$id, time.s = 5, data = d1.long.complete) # no A0s = {1,2}
trial6 <- seq.trial.data.gen2(id.list = eli6$id, time.s = 6, data = d1.long.complete) # no A0s = {1,2}


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
  group_by(trial, id) %>% 
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
# t.hor = seq(0,4,0.01)
d.time0 <- d1.long.complete %>% filter(time == 0)

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

# Some plots
par(mfrow = c(1,2))
plot(t.hor, orginal.weight.results$surv0, "l", ylim = c(0.5, 1), main = "no trunc")
lines(t.hor, orginal.weight.results$surv1, col = "red")
lines(t.hor, orginal.weight.results$surv2, col = "green")
lines(t.hor, orginal.weight.results$surv3, col = "orange")
legend("bottomright", 
       col = c("black", "red", "green", "orange"), 
       legend = c("No transfusions", "Low", "High", "Middle"),
       lty = 1,
       lwd = 2, box.lwd = "o")

plot(t.hor, trunc95.weight.results$surv0, "l", ylim = c(0.5, 1), main = "trunc 95")
lines(t.hor, trunc95.weight.results$surv1, col = "red")
lines(t.hor, trunc95.weight.results$surv2, col = "green")
lines(t.hor, trunc95.weight.results$surv3, col = "orange")

# Figure 2b:Estimated survival difference between the liberal and restrictive transfusion strategies, 
# with bootstrap confidence interval

# load bootstrap result
load("bootstrap_result.Rdata")

boot.data <- bootstrap.results %>% 
  pivot_longer(-c(time, boot.id))
surv.differce.sw.boot <- bootstrap.results %>% 
  group_by(time) %>% 
  summarise(diff21.5 = quantile(surv2.org - surv1.org, 0.05),
            diff21.95 = quantile(surv2.org - surv1.org, 0.95),
            diff21.truc95.5 = quantile(surv2.t95 - surv1.t95, 0.05),
            diff21.truc95.95 = quantile(surv2.t95 - surv1.t95, 0.95))

plot.d <- surv.differce.sw.boot %>% 
  mutate(diff = orginal.weight.results$surv2 - orginal.weight.results$surv1,
         diff.t95 = trunc95.weight.results$surv2 - trunc95.weight.results$surv1)

plot.d %>% 
  ggplot(aes(time + 1, diff)) +
  geom_ribbon(aes(ymin = diff21.5, ymax = diff21.95), fill = "grey70")+
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Day", 
       y = "Difference in survival probability") +
  scale_x_continuous(breaks = 1:7) +
  ylim(-0.4, 0.4)
