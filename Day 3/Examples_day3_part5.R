### Examples CDC training - day 3 - part 5 

library(INLA)
library(INLAjoint)
library(JM) # This package contains the dataset

data(pbc2) # dataset
#pbc2 <- pbc2[which(pbc2$id %in% c(1:5)),]
# extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age",
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])
Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]
Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2


# explain empirical Bayes // paper gives same results with eb and default

# Joint model with functions of time and current value + current slope association
f1 <- function(x) x^2
f2 <- function(x) x^3

M4 <- joint(formSurv =  inla.surv(time = years, event = death)  ~ drug,
            formLong = serBilir ~ (1 + year + f1(year) + f2(year))*drug +
              (f1(year) + f2(year) |id), family = "lognormal",
            dataLong = Longi, dataSurv = Surv, id = "id", timeVar = "year", assoc = "CV_CS",
            basRisk = "rw2", NbasRisk=25, control=list(int.strategy="eb", safemode=F))
summary(M4)
plot(M4, sdcor=T)


# Joint model for longitudinal and competing risk
M6 <- joint(formSurv = list(inla.surv(time = years, event = death)  ~ sex + drug,
                            inla.surv(time = years, event = trans) ~ edema * sex),
            formLong = serBilir ~ year * (drug + sex) + (1+year|id),
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", family = "lognormal",
            basRisk = c("rw1", "rw1"), assoc = c("SRE", "SRE_ind"), control=list(int.strategy="eb", safemode=F))
summary(M6)

# Joint model for 3 longitudinal and competing risk

M7 <- joint(formLong = list(serBilir ~ year * drug + sex + (1|id),
                            platelets ~ year + f1(year) + drug + sex + (1|id),
                            albumin ~ year + f1(year) + f2(year) + drug + (1|id)),
            formSurv = list(inla.surv(time = years, event = death) ~ drug,
                            inla.surv(time = years, event = trans) ~ drug),
            dataLong = Longi, dataSurv=Surv, id = "id", corLong=TRUE, timeVar = "year",
            family = c("lognormal", "poisson", "gaussian"), basRisk = c("rw1", "rw1"),
            assoc = list(c("CV", "CV"), c("SRE", ""), c("CV_CS", "CS")),
            control=list(int.strategy="eb"))
summary(M7)


# Joint model for longitudinal and multi-state
data(SurvMS) # package INLAjoint
data(LongMS) # package INLAjoint
E12 <- inla.surv(time = SurvMS[[1]]$Tstop, event = SurvMS[[1]]$status) # transition 1->2
E13 <- inla.surv(time = SurvMS[[2]]$Tstop, event = SurvMS[[2]]$status) # transition 1->3
E23 <- inla.surv(time = SurvMS[[3]]$Tstop, truncation=SurvMS[[3]]$Tstart,
                 event =SurvMS[[3]]$status) # transition 2->3
M9 <- joint(formSurv=list(E12 ~ X, E13 ~ X, E23 ~ X),
            formLong=list(y ~ time + X + (1+time|id)),
            basRisk = c("rw2", "rw1", "exponentialsurv"), timeVar = "time",
            assoc = list(c("CV", "CV", "CV")), id="id",
            dataSurv = SurvMS, dataLong = LongMS)
summary(M9)


# Multivariate joint model
Nsplines <- ns(Longi$year, knots=c(1,4))
f1 <- function(x) predict(Nsplines, x)[,1]
f2 <- function(x) predict(Nsplines, x)[,2]
f3 <- function(x) predict(Nsplines, x)[,3]
inla.setOption(num.threads="1:1") # in case of limited random access memory! => slower as it does not uses parallel computations

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















