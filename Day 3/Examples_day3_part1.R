### Examples CDC training - day 3 - part 1 

library(INLA)
library(INLAjoint)
library(JM) # This package contains the dataset for Models 1, 3
library(smcure) # This package contains the dataset for Model 2
library(frailtypack) # This package contains the dataset for Models 5, 6

data(pbc2) # primary biliary cholangitis dataset
pbc2$death <- ifelse(pbc2$status=="dead",1,0)
pbc2$trans <- ifelse(pbc2$status=="transplanted",1,0) # competing event
# extracting time to event data
SurvData <- pbc2[c(which(diff(as.numeric(pbc2[,which(colnames(pbc2)=="id")]))==1),
                   length(pbc2[,which(colnames(pbc2)=="id")])),c("id", "years", "death","drug", "sex")]

head(SurvData)


## Model 1 - Proportional hazards model with parametric baseline hazards
M1.exp <- inla(inla.surv(time = years, event = death) ~ drug + sex, data = SurvData, family="exponentialsurv")
summary(M1.exp)

M1.wei <- inla(inla.surv(years, death) ~ drug + sex, data = SurvData, family="weibullsurv")
summary(M1.wei)


# same model with INLAjoint
M1.exp_2 <- joint(formSurv = inla.surv(years, death) ~ drug + sex, dataSurv = SurvData, basRisk="exponentialsurv")
summary(M1.exp_2)
summary(M1.exp_2, hr=TRUE) # hazard ratios instead of betas

# Weibull baseline risk
M1.wei_2 <- joint(formSurv = inla.surv(years, death) ~ drug + sex, dataSurv = SurvData, basRisk="weibullsurv")
summary(M1.wei_2) # shape = 1 => exponential
# for weibull variants, see: inla.doc("weibullsurv")
# variant=0 is the commonly used weibull survival model
# variant=1 is a more stable, slightly differently parametrized but equivalent version
# => Can be formulated as AFT!

# smooth semiparametric baseline risk (random walk order 2)
M1_3 <- joint(formSurv = inla.surv(years, death) ~ drug + sex, dataSurv = SurvData, basRisk="rw2", NbasRisk = 30)
summary(M1_3)
plot(M1_3, hr=T) # posterior marginals and baseline risk

NewData <- SurvData[c(1,3),]
NewData$years=0
NewData$death=0
P <- predict(M1_3, NewData, horizon=14, survival=TRUE, id="id")
P1 <- P$PredS[P$PredS$id==1,]
P2 <- P$PredS[P$PredS$id==3,]
plot(P1$time, P1$Haz_quant0.5, type="l")
plot(P1$time, P1$Surv_quant0.5, type="l", ylim=c(0,1), main="Survival probabilities")
lines(P1$time, P1$Surv_quant0.025, type="l", ylim=c(0,1), lty=2)
lines(P1$time, P1$Surv_quant0.975, type="l", ylim=c(0,1), lty=2)
lines(P2$time, P2$Surv_quant0.5, col=2)
lines(P2$time, P2$Surv_quant0.025, col=2, lty=2)
lines(P2$time, P2$Surv_quant0.975, col=2, lty=2)
legend("topright", c("Female", "Male"), lty=c(1,1), col=c(1,2))
summary(SurvData) # more data om females => less uncertainty


## Model 2: Competing risks
SurvData2 <- pbc2[c(which(diff(as.numeric(pbc2[,which(colnames(pbc2)=="id")]))==1),
                   length(pbc2[,which(colnames(pbc2)=="id")])),c("id", "years", "death", "trans","drug", "sex")]

M2 <- joint(formSurv = list(inla.surv(years, death) ~ drug + sex,
                            inla.surv(years, trans) ~ drug + sex),
            basRisk = c("rw2", "rw2"), dataSurv = SurvData2)
summary(M2)
summary(M2, hr=T)
plot(M2)
NewData2 <- SurvData2[c(2,3,5,14),]
NewData2$years <- 0
NewData2$death <- NewData2$trans <- 0
NewData2$id <- 1:4

P2 <- predict(M2, NewData2, horizon=14, CIF=TRUE, id="id")
# CIF is preferred because it accounts for other competing events
library(ggplot2)
ggplot(P2$PredS, aes(x=time, y=CIF_quant0.5, group=id)) +
  geom_line(aes(color=id)) +
  geom_line(aes(x=time, y=CIF_quant0.025, color=id), linetype="dashed")+
  geom_line(aes(x=time, y=CIF_quant0.975, color=id), linetype="dashed")+
  facet_wrap(~Outcome+id, ncol=4)




## Model 3: Multi-state model (e.g., illness-death)
data(SurvMS) # package INLAjoint
E12 <- inla.surv(time = SurvMS[[1]]$Tstop, event = SurvMS[[1]]$status) # transition 1->2
E13 <- inla.surv(time = SurvMS[[2]]$Tstop, event = SurvMS[[2]]$status) # transition 1->3
E23 <- inla.surv(time = SurvMS[[3]]$Tstop, truncation=SurvMS[[3]]$Tstart,
                 event =SurvMS[[3]]$status) # transition 2->3

M3 <- joint(formSurv=list(E12 ~ X, E13 ~ X, E23 ~ X),
            basRisk = c("rw2", "rw1", "exponentialsurv"),
            dataSurv = SurvMS)
summary(M3)

## Model 4: Mixture cure model
data(bmt) # package smcure
surv.obj.cure <- inla.surv(time = bmt$Time, event = bmt$Status,
                           cure = cbind("Int"=1, "TRT"=bmt$TRT))
M4 <- joint(formSurv = surv.obj.cure ~ TRT, basRisk = "weibullsurv",
            dataSurv = bmt, control = list(variant=0))
summary(M4)


## Model 5: shared frailty model for recurrent events
data(readmission) # package frailtypack
?readmission
M5 <- joint(formSurv=inla.surv(t.stop, event) ~ chemo + (1|id), id="id", dataSurv = readmission)
summary(M5)

## Model 6: shared frailty model for recurrent and terminal event
terminalData <- readmission[readmission$event==0,]
M6 <- joint(formSurv=list(inla.surv(t.stop, event) ~ chemo + (1|id),
                           inla.surv(t.stop, death) ~ chemo), id="id",
             basRisk=c("rw2", "exponentialsurv"), assocSurv=TRUE,
             dataSurv = list(readmission,terminalData))
summary(M6)











