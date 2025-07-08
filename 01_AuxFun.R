# Aux functions for oblex analysis
# author: Thao Le

###########################################
# Function to define eligible patient #####
##########################################
f.eligible <- function(x){
  eli <- NULL
  len.x = length(x)
  
  if(len.x>1){
    eli[1] = 1
    for(i in 2:len.x){
      eli[i] = ifelse(any(x[1:i-1]!=0),0,1)
    }
  }else{eli = 1}
  
  eli
}

##########################################
# Function to create trial at time s #####
##########################################
# L: any mech support, rrt, severe barc
# L0: age, surgery, ecpr, rrt, diagnostic group
seq.trial.data.gen2 <- function(id.list, time.s, data){
  data %>% 
    filter(id %in% id.list & time >= time.s) %>% 
    group_by(id) %>% 
    mutate(
      # treatment at the start of trial s
      A0.s = first(A),
      # Baseline cov at the start of trial s
      L0s.rrt = first(RRT),
      L0s.anyMS = first(anyMS),
      L0s.barc = first(SevereBARC),
      
      # Time varying covariates
      # L.rrt = RRT,
      # L.anyMS = anyMS,
      # L.barc = SevereBARC
      L.lead.rrt = lead(RRT),
      L.lead.anyms = lead(anyMS),
      L.lead.barc = lead(SevereBARC),
      
      # Indicator: not censored at the end of the current follow-up
      # i.e. next A is equal to A0.s
      Cens = ifelse(lead(A) == A0.s, 1,
                    ifelse(lead(A) == 0 & A0.s !=0, 1,0)),
      
      # Others
      rownum = row_number(),
      trial = time.s,
      new.time = time - time.s, # first ecmo is time 0
      new.timestop = time.stop - time.s,
      no.meet.cond = ifelse(A0.s == 0, !(A == A0.s), !(A == A0.s | A == 0)),
      cum.no.meet = cumsum(no.meet.cond)) %>% 
    filter(cum.no.meet<1) %>% 
    select(-c(cum.no.meet, no.meet.cond))
}

#########################################################
# Function to generate data for bootstrap procedure #####
#########################################################
seq.trial.data.gen_boot_2 <- function(id.list, time.s, data){
  data %>% 
    filter(newid %in% id.list & time >= time.s) %>% 
    group_by(newid) %>% 
    mutate(
      # treatment at the start of trial s
      A0.s = first(A),
      # Baseline cov at the start of trial s
      L0s.rrt = first(RRT),
      L0s.anyMS = first(anyMS),
      L0s.barc = first(SevereBARC),
      
      # Time varying covariates
      # L.rrt = RRT,
      # L.anyMS = anyMS,
      # L.barc = SevereBARC
      L.lead.rrt = lead(RRT),
      L.lead.anyms = lead(anyMS),
      L.lead.barc = lead(SevereBARC),
      
      # Indicator: not censored at the end of the current follow-up
      # i.e. next A is equal to A0.s
      Cens = ifelse(lead(A) == A0.s, 1,
                    ifelse(lead(A) == 0 & A0.s !=0, 1,0)),
      
      # Others
      rownum = row_number(),
      trial = time.s,
      new.time = time - time.s, # first ecmo is time 0
      new.timestop = time.stop - time.s,
      no.meet.cond = ifelse(A0.s == 0, !(A == A0.s), !(A == A0.s | A == 0)),
      cum.no.meet = cumsum(no.meet.cond)) %>% 
    filter(cum.no.meet<1) %>% 
    select(-c(cum.no.meet, no.meet.cond))
}


#########################################################
# Function to fit aalen additive hazard models #####
#########################################################
surv.aalen.func <- function(data = mega.d,
                            weight,
                            Z = Z,
                            t.hor = t.hor){
  ah.sw = aalen(Surv(new.time, new.timestop, event) ~ A0.s + 
                  age + gender + ECPR + rrt_pre + Surgical + diag.group + 
                  L0s.rrt + L0s.anyMS + L0s.barc, 
                data = data, 
                weights = weight)
  
  
  # Obtain survival probability for "always treated" and "never treated", at several time point thor
  # obtain stepfunction for A, L and intercept from the aalen model
  
  ah.sw.stepfun.int = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"(Intercept)"])) # intercept
  ah.sw.stepfun.A1 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"A0.s1"])) # treatment A1
  ah.sw.stepfun.A2 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"A0.s2"])) # treatment A2
  ah.sw.stepfun.A3 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"A0.s3"])) # treatment A3
  ah.sw.stepfun.age = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"age"])) # age
  ah.sw.stepfun.gender = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"gender"])) 
  ah.sw.stepfun.ecpr = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"ECPR"])) 
  ah.sw.stepfun.rrtpre = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"rrt_pre"])) 
  ah.sw.stepfun.surgical = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"Surgical"])) 
  ah.sw.stepfun.diag1 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"diag.groupchronic_cardio"])) 
  ah.sw.stepfun.diag2 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"diag.groupmyocarditis"]))
  ah.sw.stepfun.diag3 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"diag.groupOther_acute_cardio"])) 
  # ah.sw.stepfun.diag4 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"diag.groupperi_support"])) 
  ah.sw.stepfun.diag4 = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"diag.grouppulmonary"])) 
  ah.sw.stepfun.l0rrt = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"L0s.rrt"])) 
  ah.sw.stepfun.l0anyms = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"L0s.anyMS"])) 
  ah.sw.stepfun.l0barc = stepfun(ah.sw$cum[,1], c(0,ah.sw$cum[,"L0s.barc"])) 
  
  
  # surv probabilities
  
  ## A = 0 (no transfusion)
  surv0.seq.sw <- sapply(t.hor, FUN = function(h){
    mean(exp(-(ah.sw.stepfun.int(h) + 
                 ah.sw.stepfun.age(h)*Z[,"age"] + 
                 ah.sw.stepfun.gender(h)*Z[,"gender"] + 
                 ah.sw.stepfun.ecpr(h)*Z[,"ECPR"] + 
                 ah.sw.stepfun.rrtpre(h)*Z[,"rrt_pre"] + 
                 ah.sw.stepfun.surgical(h)*Z[,"Surgical"] + 
                 ah.sw.stepfun.diag1(h)*Z[,"diag.groupchronic_cardio"] + 
                 ah.sw.stepfun.diag2(h)*Z[,"diag.groupmyocarditis"] + 
                 ah.sw.stepfun.diag3(h)*Z[,"diag.groupOther_acute_cardio"] + 
                 # ah.sw.stepfun.diag4(h)*Z[,"diag.groupperi_support"] + 
                 ah.sw.stepfun.diag4(h)*Z[,"diag.grouppulmonary"] + 
                 ah.sw.stepfun.l0rrt(h)*Z[,"RRT"] + 
                 ah.sw.stepfun.l0anyms(h)*Z[,"anyMS"] + 
                 ah.sw.stepfun.l0barc(h)*Z[,"SevereBARC"]
    )
    )
    )
  }
  )
  ## A = 1 (low)
  surv1.seq.sw <- sapply(t.hor, FUN = function(h){
    mean(exp(-(ah.sw.stepfun.int(h) + ah.sw.stepfun.A1(h) + 
                 ah.sw.stepfun.age(h)*Z[,"age"] + 
                 ah.sw.stepfun.gender(h)*Z[,"gender"] + 
                 ah.sw.stepfun.ecpr(h)*Z[,"ECPR"] + 
                 ah.sw.stepfun.rrtpre(h)*Z[,"rrt_pre"] + 
                 ah.sw.stepfun.surgical(h)*Z[,"Surgical"] + 
                 ah.sw.stepfun.diag1(h)*Z[,"diag.groupchronic_cardio"] + 
                 ah.sw.stepfun.diag2(h)*Z[,"diag.groupmyocarditis"] + 
                 ah.sw.stepfun.diag3(h)*Z[,"diag.groupOther_acute_cardio"] + 
                 # ah.sw.stepfun.diag5(h)*Z[,"diag.groupperi_support"] + 
                 ah.sw.stepfun.diag4(h)*Z[,"diag.grouppulmonary"] + 
                 ah.sw.stepfun.l0rrt(h)*Z[,"RRT"] + 
                 ah.sw.stepfun.l0anyms(h)*Z[,"anyMS"] + 
                 ah.sw.stepfun.l0barc(h)*Z[,"SevereBARC"]
    )
    )
    )
  }
  )
  
  ## A = 2 (high)
  surv2.seq.sw <- sapply(t.hor, FUN = function(h){
    mean(exp(-(ah.sw.stepfun.int(h) + ah.sw.stepfun.A2(h) + 
                 ah.sw.stepfun.age(h)*Z[,"age"] + 
                 ah.sw.stepfun.gender(h)*Z[,"gender"] + 
                 ah.sw.stepfun.ecpr(h)*Z[,"ECPR"] + 
                 ah.sw.stepfun.rrtpre(h)*Z[,"rrt_pre"] + 
                 ah.sw.stepfun.surgical(h)*Z[,"Surgical"] + 
                 ah.sw.stepfun.diag1(h)*Z[,"diag.groupchronic_cardio"] + 
                 ah.sw.stepfun.diag2(h)*Z[,"diag.groupmyocarditis"] + 
                 ah.sw.stepfun.diag3(h)*Z[,"diag.groupOther_acute_cardio"] + 
                 # ah.sw.stepfun.diag4(h)*Z[,"diag.groupperi_support"] + 
                 ah.sw.stepfun.diag4(h)*Z[,"diag.grouppulmonary"] + 
                 ah.sw.stepfun.l0rrt(h)*Z[,"RRT"] + 
                 ah.sw.stepfun.l0anyms(h)*Z[,"anyMS"] + 
                 ah.sw.stepfun.l0barc(h)*Z[,"SevereBARC"]
    )
    )
    )
  }
  )
  
  ## A = 3 (middle)
  surv3.seq.sw <- sapply(t.hor, FUN = function(h){
    mean(exp(-(ah.sw.stepfun.int(h) + ah.sw.stepfun.A3(h) + 
                 ah.sw.stepfun.age(h)*Z[,"age"] + 
                 ah.sw.stepfun.gender(h)*Z[,"gender"] + 
                 ah.sw.stepfun.ecpr(h)*Z[,"ECPR"] + 
                 ah.sw.stepfun.rrtpre(h)*Z[,"rrt_pre"] + 
                 ah.sw.stepfun.surgical(h)*Z[,"Surgical"] + 
                 ah.sw.stepfun.diag1(h)*Z[,"diag.groupchronic_cardio"] + 
                 ah.sw.stepfun.diag2(h)*Z[,"diag.groupmyocarditis"] + 
                 ah.sw.stepfun.diag3(h)*Z[,"diag.groupOther_acute_cardio"] + 
                 # ah.sw.stepfun.diag4(h)*Z[,"diag.groupperi_support"] + 
                 ah.sw.stepfun.diag4(h)*Z[,"diag.grouppulmonary"] + 
                 ah.sw.stepfun.l0rrt(h)*Z[,"RRT"] + 
                 ah.sw.stepfun.l0anyms(h)*Z[,"anyMS"] + 
                 ah.sw.stepfun.l0barc(h)*Z[,"SevereBARC"]
    )
    )
    )
  }
  )
  return(list(fit = ah.sw,
              surv0 = surv0.seq.sw,
              surv1 = surv1.seq.sw,
              surv2 = surv2.seq.sw,
              surv3 = surv3.seq.sw))
}


