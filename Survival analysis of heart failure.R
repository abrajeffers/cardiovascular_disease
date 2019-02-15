library(survival)
library(randomForestSRC)
library(survivalsvm)
library(Hmisc)

setwd("/Users/abrajeffers/programming_projects/cardiovascular_disease/")
df = read.csv("./S1Data.csv")
View(df)
event_horizon = 60
df$ReachedEvent <- rep(0,nrow(df))
df$ReachedEvent[which(df$Event == 1 & df$TIME < event_horizon)] <- 1
mean(df$ReachedEvent)

df$index <- rep(NA,nrow(df))
for(i in 1:nrow(df)){
  df$index[i] = i
}

survfit_KM <- survfit(Surv(df$TIME,df$Event)~1)
plot(survfit_KM, conf.int = FALSE, mark.time = FALSE, main = "KM of Congestive Heart Failure Patients",xlab = "Time in days",ylab = "Survival function", lwd = 1)


summary(df)

hist(df$CPK,freq = FALSE)




h.CPK <- hist(df$CPK, breaks = 100, plot=FALSE)
h.CPK$counts=h.CPK$counts/sum(h.CPK$counts)
plot(h.CPK, ylab = "Probability",xlab = "CPK",ylim = c(0,1), main ='Histogram of CPK')

h.age <- hist(df$Age,breaks = seq(40,100,5),plot =FALSE)
h.age$counts = h.age$counts/sum(h.age$counts)
plot(h.age, ylab = "Probability",xlab = "Age",ylim = c(0,1),xlim = c(40,100), main ='Histogram of age')

h.ef <- hist(df$Ejection.Fraction,breaks = seq(10,80,5), plot = FALSE)
h.ef$counts = h.ef$counts/sum(h.ef$counts)
plot(h.ef, ylab = "Probability",xlab = "Ejection fraction",ylim = c(0,1), main ='Histogram of ejection fraction')

h.sodium <- hist(df$Sodium, breaks = seq(110,150,2),plot = FALSE)
h.sodium$counts = h.sodium$counts/sum(h.sodium$counts)
plot(h.sodium, ylab = "Probability",xlab = "Sodium",ylim = c(0,1), main ='Histogram of sodium')


rf.df = rfsrc(Surv(TIME, Event) ~ Gender + Smoking + Diabetes + BP + Anaemia + Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK, data = df, seed = 123, importance = "random")

rf.df

sort(rf.df$importance)

df$chf.oob = rep(0,nrow(df))
for(i in 1:nrow(df)){
  df$chf.oob[i] <- rf.df$chf.oob[i,which(rf.df$time.interest == event_horizon)]
}

for(i in 1:nrow(df)){
  df$chf.inbag[i] <- rf.df$chf[i,which(rf.df$time.interest == event_horizon)]
}

df$r <- rep(0,nrow(df))
for(i in 1:nrow(df)){
  df$r[i] = 1-exp(-rf.df$chf.oob[i,which(rf.df$time.interest == event_horizon)])
}
h.r <- hist(df$r,breaks = seq(0,1,.05),plot = FALSE)
h.r$counts = h.r$counts/sum(h.r$counts)
plot(h.r, ylab = "Probability",xlab = "Assigned risk",ylim = c(0,1),xlim = c(0,1), main ='Histogram of random forest assigned risk')


#quartiles of assigned risk
r_quartile <- stats::quantile(df$r,seq(0,1,by=1/4))
df$k_r <- (df$r > r_quartile[2]) + (df$r > r_quartile[3]) + (df$r > r_quartile[4]) +1

#survival object and survfit for each quartile's dataframe

survfit_r1 <- survfit(Surv(df$TIME[which(df$k_r == 1)],df$Event[which(df$k_r == 1)])~1)
survfit_r2 <- survfit(Surv(df$TIME[which(df$k_r == 2)],df$Event[which(df$k_r == 2)])~1)
survfit_r3 <- survfit(Surv(df$TIME[which(df$k_r == 3)],df$Event[which(df$k_r == 3)])~1)
survfit_r4 <- survfit(Surv(df$TIME[which(df$k_r == 4)],df$Event[which(df$k_r == 4)])~1)

survfit_all <- survfit(Surv(time = df$TIME, event = df$Event == 1) ~ 1, type = "kaplan-meier")

#function which takes survfit object and endpt (t=60) and returns survival probability
#at the greatest failure time <= endpt
FFpi <- function(data, endpt) {
  idx <- which(summary(data)$time <= endpt)
  return(summary(data)$surv[idx[length(idx)]])
}


#function that returns the upper confidence limit of the survival probability
FFupper <- function(data, endpt) {
  idx <- which(summary(data)$time <= endpt)
  return(summary(data)$upper[idx[length(idx)]])
}

#function that returns the lower confidnece limit of the survival probability
FFlower <- function(data, endpt) {
  idx <- which(summary(data)$time <= endpt)
  return(summary(data)$lower[idx[length(idx)]])
}

#function that returns the standard error of the survival probability at the greatest
#failure time elss than or equal to the endpt
FFse <- function(data, endpt) {
  idx <- which(summary(data)$time <= endpt)
  return(summary(data)$std.err[idx[length(idx)]])
}

#pi_i is risk of event (death) in group i
pi1 <- 1- FFpi(survfit_r1,event_horizon)
pi2 <- 1 - FFpi(survfit_r2,event_horizon)
pi3 <- 1 - FFpi(survfit_r3,event_horizon)
pi4 <- 1 - FFpi(survfit_r4,event_horizon)

pi_all <- 1- FFpi(survfit_all, event_horizon)
se_all <- FFse(survfit_all,event_horizon)

#robserved risk of death
pi <- c(pi1,pi2,pi3,pi4)
se <- c(FFse(survfit_r1,event_horizon),FFse(survfit_r2,event_horizon),FFse(survfit_r3,event_horizon),FFse(survfit_r4,event_horizon))

x.RSF <- c(mean(df[which(df$k_r ==1),]$r),mean(df[which(df$k_r ==2),]$r),mean(df[which(df$k_r ==3),]$r),mean(df[which(df$k_r ==4),]$r))
y.RSF <- c(pi1,pi2,pi3, pi4) 

#Attribute Diagram for RSF
errbar(x.RSF,y.RSF,y.RSF-1.96*se, y.RSF+1.96*se, xlim =c(0,1), ylim=c(0,1), main = 'Attribute diagram for RSF',xlab="Assigned Risk", ylab="Observed Risk")
abline(0,1)
title('Attribute diagram for RSF risks')

# set.seed(123)
# random_splits <- runif(nrow(df))
# train_df_official <- df[random_splits < .5,]
# dim(train_df_official)
# 
# validate_df_official <- df[random_splits >= .5,]
# dim(validate_df_official)
# 
# train_df_classification  <- train_df_official 
# validate_df_classification  <- validate_df_official

cases <- df[which(df$ReachedEvent == 1),]
controls <- df[which(df$TIME > event_horizon),]
auc <- rbind(cases,controls)

my.roc.rsf <- roc(auc$Event,auc$r,ci = TRUE)
my.auc.rsf <- auc(my.roc.rsf)
my.auc.rsf
ci.auc(my.auc.rsf)

############################Standard Cox PH model with bootstrapping##########
#Number of bootstraps
B = 1000
set.seed(123)
inBagIndex <- data.frame(matrix(NA, nrow = nrow(df), ncol = B))
for(i in 1:B){
  inBagIndex[,i] = sample(seq(1,nrow(df),1),nrow(df), replace = TRUE)
}
inBag<-  rep(list(data.frame(matrix(NA,nrow = nrow(df),ncol = ncol(df)))),B)
for(i in 1:B){
  inBag[[i]] = df[inBagIndex[,i],]
}
outBag = rep(list(data.frame()),B)
for(i in 1:B){
  outBag[[i]] = df[-inBagIndex[,i],]
}
Betas.df <- data.frame(matrix(NA, nrow = B, ncol = 11))
for(i in 1:B){
  cox.fit <- coxph(Surv(TIME, Event) ~ Gender + Smoking + Diabetes + BP + Anaemia + Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK, data = inBag[[i]])
  bh=basehaz(cox.fit)
  bh.eventhorizon = bh$hazard[max(which(bh$time <= event_horizon))]
  
  outBag[[i]]$r.cox = rep(0,nrow(outBag[[i]]))
  inBag[[i]]$r.cox = rep(0,nrow(inBag[[i]]))
  for(j in 1:11){
    Betas.df[i,j] = cox.fit$coefficients[j]
  }
  for(j in 1:nrow(outBag[[i]])){
    innerproduct= 0
    for(k in 1:11){
      innerproduct = innerproduct +  exp(cox.fit$coefficients[k]*outBag[[i]][j,k+2])
    }
    outBag[[i]]$r.cox[j] =exp(-bh.eventhorizon)^innerproduct  
  }
  # for(j in 1:nrow(inBag[[i]])){
  #   innerproduct= 0
  #   for(k in 1:11){
  #     iBinnerproduct = iBinnerproduct +  exp(cox.fit$coefficients[k]*inBag[[i]][j,k+2])
  #   }
  #   inBag[[i]]$r.cox[j] = exp(-bh.eventhorizon)^IBinnerproduct
  # }
}
Betas = colMeans(Betas.df)
ci = data.frame(matrix(NA,nrow = 2, ncol = length(Betas)))
for(i in 1:length(Betas)){
  ci[1,i] = quantile(Betas.df[,i],0.05,na.rm = TRUE)
  ci[2,i] = quantile(Betas.df[,i],0.95,na.rm = TRUE)
}
rlist = vector("list", length = nrow(df))
for(r in 1:nrow(df)){
  rlist[[r]] = r
}
for(i in 1:B){
  for(j in 1:nrow(df)){
    for(k in 1:nrow(outBag[[i]])){
      if(outBag[[i]][k,]$index == rlist[[j]][1]){
        rlist[[j]] = list.append(rlist[[j]],outBag[[i]][k,]$r.cox)
      }
    }
  }
}
#############


df$r.cox = rep(NA,nrow(df))
for(i in 1:length(rlist)){
  rmean = mean(rlist[[i]][-1])
  df[i,]$r.cox = rmean
}
par(mfrow=c(1,2)) 
plot(h.r, ylab = "Probability",xlab = "Assigned risk",ylim = c(0,1),xlim = c(0,1), main ='Random Forest')
h.r.cox <- hist(df$r.cox,breaks = seq(0,1,.05),plot = FALSE)
h.r.cox$counts = h.r.cox$counts/sum(h.r.cox$counts)
plot(h.r.cox, ylab = "Probability",xlab = "Assigned risk",ylim = c(0,1),xlim = c(0,1), main ='Cox PH')

#quartiles of assigned risk
r_quartile.cox <- stats::quantile(df$r.cox,seq(0,1,by=1/4))
df$k_r.cox <- (df$r.cox > r_quartile.cox[2]) + (df$r.cox > r_quartile.cox[3]) + (df$r.cox > r_quartile.cox[4]) +1

#survival object and survfit for each quartile's dataframe

survfit_r1.cox <- survfit(Surv(df$TIME[which(df$k_r.cox == 1)],df$Event[which(df$k_r.cox == 1)])~1)
survfit_r2.cox <- survfit(Surv(df$TIME[which(df$k_r.cox == 2)],df$Event[which(df$k_r.cox == 2)])~1)
survfit_r3.cox <- survfit(Surv(df$TIME[which(df$k_r.cox == 3)],df$Event[which(df$k_r.cox == 3)])~1)
survfit_r4.cox <- survfit(Surv(df$TIME[which(df$k_r.cox == 4)],df$Event[which(df$k_r.cox == 4)])~1)

#pi_i is risk of event (death) in group i
pi1.cox <- 1- FFpi(survfit_r1.cox,event_horizon)
pi2.cox <- 1 - FFpi(survfit_r2.cox,event_horizon)
pi3.cox <- 1 - FFpi(survfit_r3.cox,event_horizon)
pi4.cox <- 1 - FFpi(survfit_r4.cox,event_horizon)

#robserved risk of death
pi <- c(pi1.cox,pi2.cox,pi3.cox,pi4.cox)
se.cox <- c(FFse(survfit_r1.cox,event_horizon),FFse(survfit_r2.cox,event_horizon),FFse(survfit_r3.cox,event_horizon),FFse(survfit_r4.cox,event_horizon))

x.cox <- c(mean(df[which(df$k_r.cox ==1),]$r.cox),mean(df[which(df$k_r.cox ==2),]$r.cox),mean(df[which(df$k_r.cox ==3),]$r.cox),mean(df[which(df$k_r.cox ==4),]$r.cox))
y.cox <- c(pi1.cox,pi2.cox,pi3.cox, pi4.cox) 

par(mfrow=c(1,1)) 
#Attribute Diagram for Cox PH
errbar(x.cox,y.cox,y.cox-1.96*se.cox, y.cox+1.96*se.cox, xlim =c(0,1), ylim=c(0,1), main = 'Cox PH',xlab="Assigned Risk", ylab="Observed Risk")
abline(0,1)
title('Attribute diagram for Cox PH risks')

par(mfrow=c(1,1)) 
errbar(x.RSF,y.RSF,y.RSF-1.96*se, y.RSF+1.96*se, xlim =c(0,1), ylim=c(0,1), main = 'Attribute diagram for RSF',xlab="Assigned Risk", ylab="Observed Risk")
abline(0,1)
title('Attribute diagram for RSF risks')


cases <- df[which(df$ReachedEvent == 1),]
controls <- df[which(df$TIME > event_horizon),]
auc <- rbind(cases,controls)

##AUC and CL for AUC
my.roc.cox <- roc(auc$Event, auc$r.cox,ci = TRUE)
my.auc.cox <- auc(my.roc.cox)
print("Cox PH:")
my.auc.cox
ci.auc.cox = ci.auc(my.auc.cox)


print("Random Survival Forest:")
my.auc.rsf
ci.auc.rsf = ci.auc(my.auc.rsf)


par(mfrow=c(1,1)) 
plot.roc(my.roc.rsf,title="Receiver operator curve")
lines.roc(my.roc.cox, print.auc = TRUE, col = "red")
legend("bottomright",legend = c(paste("RSF. AUC=",round(my.auc.rsf,3),", 95%CI (",round(ci.auc.rsf[1],3),",",round(ci.auc.rsf[2],3),")",sep=""),paste("Cox PH. AUC=",round(my.auc.cox,3)," 95%CI (",round(ci.auc.cox[1],3),", ",round(ci.auc.cox[2],3),")",sep = "")),lty = 1,cex = .7,col = c("black","red"))
title("Receiver operator curve",line = 2.5)

