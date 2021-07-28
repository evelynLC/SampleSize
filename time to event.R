library(tidyr)
library(dplyr)
library("survival")

# Example: DREAMM 15
# Design: Du-endpoints(mPFS + MDR rate)
# mPFS = 52m vs 36m, HR=0.692, alpha=0.02 one side
# MDR% = 40 vs 20, diff=0.2, alpha=0.005 one side
# Global sample size 750, 1:1, 31pts/m, drop out = 8%/year
# 

## mPFS
## POWER for China

subject.sim <- function(n, enrollment, follow, dropout){
  # population
  enroll.time <- runif(n)*enrollment
  trt.group <- c(rep(c(1,2), n/2))
  pop <- data.frame(cbind(enroll.time, trt.group)) %>% arrange(trt.group)
  
  # survival
  time.t <- rexp(n/2, rate=log(2)/mt.t)
  time.p <- rexp(n/2, rate=log(2)/mt.p)
  survtime<- c(time.t, time.p)
  
  # drop
  time.drop <- rexp(n, rate=dropout)
  
  # combine
  ana.time <- c(NA)
  event <- c(NA)
  ana.time2 <- c(NA)
  event2 <- c(NA)
  for (i in 1:n){
    ana.time[i] <- min(time.drop[i], survtime[i], enrollment+follow - enroll.time[i]) 
    event[i] <- ifelse(ana.time[i]==survtime[i], 1, 0)
    # if accrual start delay for 6 month
    ana.time2[i] <- min(time.drop[i], survtime[i], enrollment+(follow-6) - enroll.time[i]) 
    event2[i]    <- ifelse(ana.time2[i]==survtime[i], 1, 0)
  }
  
  return(cbind(pop, ana.time, event, ana.time2, event2))
}

mt.p <- 12
mt.t <- 20

# observed events by trt
############################
mydata <- subject.sim(100000, 100/6, 25-100/6, 0.05/12)
mysum <- mydata %>%
  group_by(trt.group) %>%
  summarise(sum=sum(event)/1000)
############################

# power calculation
########################################
sumfit <- c()
event.num <- c()
upper <- c()

for (i in 1:2000){
  
  temp <- subject.sim(750, 750/31, 23, 0.08/12) %>% 
    mutate(trial=i)
  
  fit <- coxph(Surv(ana.time, event) ~ as.factor(trt.group), temp, 
               method="breslow")
  sumfit[i]<- 1/exp(fit$coefficients) 
  upper[i] <- 1/exp(confint(fit, level = 0.96)[1])
  event.num[i] <- fit$nevent
  
}

length(upper[upper<=1])/2 #power = ~85%
########################################


# simulation
########################################
gensurvival <- function(trial.n, study.pop, cn.pop, enrollment, follow, drop) {
  hr.cn <- c(NA)
  hr.all<- c(NA)
  event.cn <- c(NA)
  event.all<- c(NA)
  upper.all <- c(NA)  
  
  
  for (i in 1:trial.n) {
    
    sample.all <- subject.sim(study.pop, enrollment, follow, drop) 
    fit <- coxph(Surv(ana.time, event) ~ as.factor(trt.group), sample.all, 
                 method="breslow")
    hr.all[i]<- 1/exp(fit$coefficients) 
    upper.all[i] <- 1/exp(confint(fit, level = 0.96)[1])
    event.all[i] <- fit$nevent
    
    sample.cn <- sample.all[sample(1:study.pop, cn.pop),]
    fit <- coxph(Surv(ana.time2, event2) ~ as.factor(trt.group), sample.cn, 
                 method="breslow")
    hr.cn[i]<- 1/exp(fit$coefficients) 
    event.cn[i] <- fit$nevent
  }      
  
  return(as.data.frame(cbind(hr.all, upper.all, event.all, hr.cn, event.cn, 
                             all.success=ifelse(upper.all<1,1,0),
                             cn_trend=ifelse(hr.cn<1 & upper.all<1 ,1,0),
                             cn_half_eff=ifelse(log(hr.cn)<log(hr.all)/2 & upper.all<1,1,0))
  ))
}


# 15%: 750*0.15=112.5
s1 <- gensurvival(2000, 750, 112, 750/31, 23, 0.08/12)

sum(s1$cn_trend)/sum(s1$all.success)
sum(s1$cn_half_eff)/sum(s1$all.success)

# 10%: 750*10% = 75
s2 <- gensurvival(2000, 750, 76, 750/31, 23, 0.08/12)

sum(s2$cn_trend)/sum(s2$all.success)
sum(s2$cn_half_eff)/sum(s2$all.success)