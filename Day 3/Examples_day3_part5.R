### Examples CDC training - day 3 - part 5

library(INLA)
library(INLAjoint)
library(splines)

setwd("/home/dr/Documents/CDC24")
load("Data/Data_day3.RData") # contains 3 datasets (pbc2 from package JM, bmt from package smcure and readmission from package frailtypack)
pbc2 <- pbc2[which(pbc2$id %in% c(1:5)),]
# extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age",
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])
Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]
Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2


# empirical Bayes: same frequentist properties / cheaper computational cost

# Joint model with functions of time and current value + current slope association
Nsplines <- ns(Longi$year, knots=1)
f1 <- function(x) predict(Nsplines, x)[,1]
f2 <- function(x) predict(Nsplines, x)[,2]
curve(f1, from=0, to=max(Longi$year), ylim=c(-1,1))
curve(f2, from=0, to=max(Longi$year), add=T)
# f1 <- function(x) x^2
# f2 <- function(x) x^3

M1 <- joint(formSurv =  inla.surv(time = years, event = death)  ~ drug,
            formLong = serBilir ~ (1 + f1(year) + f2(year))*drug +
              (1 + f1(year) + f2(year)|id), family = "lognormal",
            dataLong = Longi, dataSurv = Surv, id = "id", timeVar = "year", assoc = "CV",
            basRisk = "rw2", control=list(int.strategy="eb"))
summary(M1)
plot(M1, sdcor=T)

NewData <- Longi[Longi$id %in%c(2, 13),]
P <- predict(M1, NewData, horizon=14, inv.link=T, survival=T)#, Csurv=0 (enter at risk)

library(ggplot2)
theme_set(theme_minimal())
ggplot(P$PredL) +
  geom_line(aes(x=year, y=quant0.5, group=id, color=id)) +
  geom_line(aes(x=year, y=quant0.025, group=id, color=id), linetype="dashed")+
  geom_line(aes(x=year, y=quant0.975, group=id, color=id), linetype="dashed")+
  theme(legend.position = "none")+
  geom_point(data = NewData, aes(x = year, y = serBilir, group=id, color=id))

ggplot(P$PredS) +
  geom_line(aes(x=year, y=Surv_quant0.5, group=id, color=id)) +
  geom_line(aes(x=year, y=Surv_quant0.025, group=id, color=id), linetype="dashed")+
  geom_line(aes(x=year, y=Surv_quant0.975, group=id, color=id), linetype="dashed")+
  ylab("Survival probability") + theme(legend.position = "none")








# id 1:5
M1_ <- joint(formSurv =  inla.surv(time = years, event = death)  ~ drug,
            formLong = serBilir ~ (1 + f1(year) + f2(year))*drug +
              (1 + f1(year) + f2(year)|id), family = "lognormal",
            dataLong = Longi, dataSurv = Surv, id = "id", timeVar = "year", assoc = "SRE_ind",
            basRisk = "rw2", control=list(int.strategy="eb"))
summary(M1_)
M1_$.args$data






















# Joint model for longitudinal and competing risk
M6 <- joint(formSurv = list(inla.surv(time = years, event = death)  ~ sex + drug,
                            inla.surv(time = years, event = trans) ~ edema * sex),
            formLong = serBilir ~ year * (drug + sex) + (1+year|id),
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", family = "lognormal",
            basRisk = c("rw1", "rw1"), assoc = c("SRE", "SRE_ind"), control=list(int.strategy="eb", safemode=F))
summary(M6)





# Joint model for longitudinal and multi-state
data(SurvMS) # package INLAjoint
data(LongMS) # package INLAjoint
E12 <- inla.surv(time = SurvMS[[1]]$Tstop, event = SurvMS[[1]]$status) # transition 1->2
E13 <- inla.surv(time = SurvMS[[2]]$Tstop, event = SurvMS[[2]]$status) # transition 1->3
E23 <- inla.surv(time = SurvMS[[3]]$Tstop, truncation=SurvMS[[3]]$Tstart,
                 event =SurvMS[[3]]$status) # transition 2->3
M7 <- joint(formSurv=list(E12 ~ X, E13 ~ X, E23 ~ X),
            formLong=list(y ~ time + X + (1+time|id)),
            basRisk = c("rw2", "rw1", "exponentialsurv"), timeVar = "time",
            assoc = list(c("CV", "CV", "CV")), id="id",
            dataSurv = SurvMS, dataLong = LongMS,
            cutpoints=seq(0, max(SurvMS[[3]]$Tstop), len=15))
summary(M7)



# Joint model for 3 longitudinal and competing risk
M8 <- joint(formLong = list(serBilir ~ year * drug + sex + (1|id),
                            platelets ~ year + f1(year) + drug + sex + (1|id),
                            albumin ~ year + f1(year) + f2(year) + drug + (1|id)),
            formSurv = list(inla.surv(time = years, event = death) ~ drug,
                            inla.surv(time = years, event = trans) ~ drug),
            dataLong = Longi, dataSurv=Surv, id = "id", corLong=TRUE, timeVar = "year",
            family = c("lognormal", "poisson", "gaussian"), basRisk = c("rw1", "rw1"),
            assoc = list(c("CV", "CV"), c("SRE", ""), c("CV_CS", "CS")),
            control=list(int.strategy="eb"))
summary(M8)




# Multivariate joint model
# From: https://doi.org/10.1093/biostatistics/kxad019
Nsplines <- ns(Longi$year, knots=c(1,4))
f1 <- function(x) predict(Nsplines, x)[,1]
f2 <- function(x) predict(Nsplines, x)[,2]
f3 <- function(x) predict(Nsplines, x)[,3]
# inla.setOption(num.threads="1:1") # in case of limited random access memory! => slower as it does not uses parallel computations

M16 <-joint(formSurv = list(inla.surv(time = years, event = death) ~ drug,
                            inla.surv(time = years, event = trans) ~ drug),
            formLong = list(serBilir ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                              (1 + f1(year) + f2(year) + f3(year) | id),
                            SGOT ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                              (1 + f1(year) + f2(year) + f3(year) | id),
                            albumin ~ (1 + year) * drug + (1 + year | id),
                            platelets ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                              (1 + f1(year) + f2(year) + f3(year) | id),
                            spiders ~ (1 + year) * drug + (1 + year | id)),
            dataLong = Longi, dataSurv = Surv, id = "id", timeVar = "year", corLong = F,
            family = c("lognormal", "lognormal", "gaussian", "poisson", "binomial"),
            basRisk = c("rw2", "rw1"), NbasRisk = 15, assoc = list(c("CV_CS", "CV"),
             c("CV", ""), c("CV", "CV"), c("CV", "CV"), c("CV", "")),
            control=list(int.strategy="eb"))
summary(M16)

# 2000 seconds on laptop with 1 thread

browseVignette("INLAjoint")













