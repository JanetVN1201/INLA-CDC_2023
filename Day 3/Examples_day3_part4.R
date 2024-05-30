### Examples CDC training - day 3 - part 4

library(INLA)
library(INLAjoint)
# Example: Joint longitudinal - Survival

set.seed(1)
# data generation - one longitudinal marker
nsujet=500 # number of individuals
# L1 (longitudinal continuous)
b_0=2 # intercept
b_1=-0.3 # slope
b_e=0.1 # residual error
# S1 (survival)
phi_b0=1 # random intercept association with survival
gap=1 # gap between measurements
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # total number of measurements
id<-rep(1:nsujet, each=nmesindiv) # individual id
# random effect
b_i <- rnorm(nsujet, 0, 0.4)
b_i0 <- rep(b_i, each=nmesindiv) # random intercept
# linear predictor
linPredL1 <- b_i0 + b_0 + b_1 * time
# continuous outcome L1
L1 <- rnorm(nmesy, linPredL1, b_e)
lon <- data.frame(L1, id, time)
## generation of exponential death times
u <- runif(nsujet) # uniform distribution for survival times generation
baseScale=0.1
deathTimes <- -(log(u) / (baseScale * exp(b_i*phi_b0)))
d <- as.numeric(deathTimes<followup) # deathtimes indicator
## censoring individuals at end of follow-up (not at random)
deathTimes[deathTimes>=followup]=followup
surv <- data.frame(id=1:nsujet,deathTimes, d) # survival times dataset
## removing longi measurements after death
ind <- rep(NA, nsujet*length(mestime))
for (i in 1:nsujet){
  for(j in 1:length(mestime)){
    if(lon[(i-1)*length(mestime)+j, "time"]<=surv[i,"deathTimes"]) ind[(i-1)*length(mestime)+j]=1
  }
}
lon <- lon[!is.na(ind),]
summary(surv)
summary(lon)
head(surv)
head(lon, 20)

#
#
# # with INLA:
# NL <- dim(lon)[1]
# NS <- dim(surv)[1]
#
# data <- list(
#   IntL1 = c(rep(1, NL), rep(NA, NS)),
#   TimeL1 = c(lon$time, rep(NA, NS)),
#   IntS1 = c(rep(NA, NL), rep(1, NS)),
#   b_i0L1 = c(lon$id, rep(NA, NS)),
#   b_i0S1 = c(rep(NA, NL), surv$id),
#   Yjoint = list(L1 = c(lon$L1, rep(NA, NS)),
#                 S1 = inla.surv(time = c(rep(NA, NL), surv$deathTimes),
#                                event = c(rep(NA, NL), surv$d)))
# )
# formula <- Yjoint ~ -1 + IntL1 + TimeL1 + IntS1 +
#   f(b_i0L1, model="iid") + f(b_i0S1, copy="b_i0L1", fixed=FALSE)
#
# JointM <- inla(formula = formula, data = data,
#                family=c("gaussian", "exponentialsurv"),
#                control.inla = list(int.strategy="eb"))
# summary(JointM)
#
# lambda <- inla.tmarginal(function(x) exp(x), JointM$marginals.fixed$IntS1)
# inla.zmarginal(lambda)
# ResErr <- inla.tmarginal(function(x) sqrt(1/x),
#                          JointM$marginals.hyperpar$`Precision for the Gaussian observations`)
# inla.zmarginal(ResErr)
# SD_b <- inla.tmarginal(function(x) sqrt(1/x),
#                        JointM$marginals.hyperpar$`Precision for b_i0L1`)
# inla.zmarginal(SD_b)


library(INLAjoint)
JointM_2 <- joint(formLong = L1 ~ time + (1|id),
                  formSurv = inla.surv(deathTimes, d) ~ 1, basRisk = "rw1",
                  id="id", timeVar = "time", assoc = "SRE",
                  dataLong = lon, dataSurv = surv)
summary(JointM_2, sdcor=TRUE)




# different parameterizations
JointM_SRE_ind <- joint(formLong = L1 ~ time + (1+time|id),
                        formSurv = inla.surv(deathTimes, d) ~ 1, basRisk = "rw1",
                        id="id", timeVar = "time", assoc = "SRE_ind",
                        dataLong = lon, dataSurv = surv)
summary(JointM_SRE_ind, sdcor=TRUE)

JointM_SRE <- joint(formLong = L1 ~ time + (1+time|id),
                    formSurv = inla.surv(deathTimes, d) ~ 1, basRisk = "rw1",
                    id="id", timeVar = "time", assoc = "SRE",
                    dataLong = lon, dataSurv = surv)
summary(JointM_SRE, sdcor=TRUE)

JointM_CV <- joint(formLong = L1 ~ time + (1+time|id),
                   formSurv = inla.surv(deathTimes, d) ~ 1, basRisk = "rw1",
                   id="id", timeVar = "time", assoc = "CV",
                   dataLong = lon, dataSurv = surv)
summary(JointM_CV, sdcor=TRUE)

JointM_CS <- joint(formLong = L1 ~ time + (1+time|id),
                   formSurv = inla.surv(deathTimes, d) ~ 1, basRisk = "rw1",
                   id="id", timeVar = "time", assoc = "CS",
                   dataLong = lon, dataSurv = surv)
summary(JointM_CS, sdcor=TRUE)

JointM_CV_CS <- joint(formLong = L1 ~ time + (1+time|id),
                      formSurv = inla.surv(deathTimes, d) ~ 1, basRisk = "rw1",
                      id="id", timeVar = "time", assoc = "CV_CS",
                      dataLong = lon, dataSurv = surv)
summary(JointM_CV_CS, sdcor=TRUE)




# predictions
# create 2 new individuals with different profiles
DN <- function(x) dnorm(x, mean = 0, sd = 0.4)
curve(DN, from=-2, to=2)
abline(v=-1, col="red")
abline(v=1, col="red")

ID1 <- 1 + b_0 + b_1 * seq(0, 1, len=5) + rnorm(5, 0, 0.1)
ID2 <- -1 + b_0 + b_1 * seq(0, 1, len=5) + rnorm(5, 0, 0.1)
NewData <- data.frame("L1"=c(ID1, ID2), "id"=c(rep(1, length(ID1)), rep(2, length(ID2))), "time"=rep(seq(0, 1, len=5), 2))

P <- predict(JointM_2, NewData, horizon=10, inv.link=T, survival=T)

library(ggplot2)
theme_set(theme_minimal())
ggplot(P$PredL) +
  geom_line(aes(x=time, y=quant0.5, group=id, color=id)) +
  geom_line(aes(x=time, y=quant0.025, group=id, color=id), linetype="dashed")+
  geom_line(aes(x=time, y=quant0.975, group=id, color=id), linetype="dashed")+
  theme(legend.position = "none")+
  geom_point(data = NewData, aes(x = time, y = L1, group=id, color=id))

TrueRisk <- data.frame("Risk" = c(baseScale*exp(1), baseScale*exp(-1)), id=1:2)

ggplot(P$PredS) + facet_wrap(~id, ncol=2) +
  geom_line(aes(x=time, y=Haz_quant0.5)) +
  geom_line(aes(x=time, y=Haz_quant0.025), linetype="dashed")+
  geom_line(aes(x=time, y=Haz_quant0.975), linetype="dashed")+
  ylab("Risk") + theme(legend.position = "none")+
  geom_hline(data=TrueRisk, aes(yintercept=Risk, group=id, color="red"))


ggplot(P$PredS) + facet_wrap(~id, ncol=2) +
  geom_line(aes(x=time, y=Surv_quant0.5)) +
  geom_line(aes(x=time, y=Surv_quant0.025), linetype="dashed")+
  geom_line(aes(x=time, y=Surv_quant0.975), linetype="dashed")+
  ylab("Survival probability") + theme(legend.position = "none")



















