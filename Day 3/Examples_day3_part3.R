### Examples CDC training - day 3 - part 3 

library(INLA)
library(INLAjoint)
# Example: Multivariate longitudinal

set.seed(1)
# data generation - two longitudinal markers
nsujet=500 # number of individuals

# Y1 (continuous)
b1_0=0.2 # intercept
b1_1=-0.1 # slope
b1_e=0.1 # residual error

# Y2 (counts)
b2_0=3 # intercept
b2_1=-0.1 # slope

gap=1 # gap between measurements
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # total number of measurements
id<-rep(1:nsujet, each=nmesindiv) # individual id

# random effects variance-covariance matrix
Sigma <- matrix(c(0.16, 0.03, 0.02, 0.04,
                  0.03, 0.09, 0.03, 0.00,
                  0.02, 0.03, 0.25, 0.08,
                  0.04, 0.00, 0.08, 0.16),ncol=4,nrow=4)

MVnorm <- mvtnorm::rmvnorm(nsujet, rep(0, 4), Sigma)
b1_i0 <- rep(MVnorm[,1], each=nmesindiv) # random intercept Y1
b1_i1 <- rep(MVnorm[,2], each=nmesindiv) # random slope Y1
b2_i0 <- rep(MVnorm[,3], each=nmesindiv) # random intercept Y2
b2_i1 <- rep(MVnorm[,4], each=nmesindiv) # random slope Y2

# linear predictor
linPredY1 <- (b1_i0+b1_0) + (b1_i1+b1_1)*time
linPredY2 <- (b2_i0+b2_0) + (b2_i1+b2_1)*time

# observed outcomes:
# continuous outcome Y1
Y1 <- rnorm(nmesy, linPredY1, b1_e)
# count outcome Y2
Y2 <- rpois(nmesy, exp(linPredY2))
lon <- data.frame(Y1, Y2, id, time)
summary(lon) # dataset


# first fit a mixed effects model for the first marker
ids <- length(unique(lon$id))+lon$id

formula <- Y1 ~ time + f(id, model="iid2d", n=length(unique(lon$id))*2)+
  f(ids, time, copy="id")
M1 <- inla(formula, data=lon)
summary(M1)


formula <- Y1 ~ time + f(id, model="iidkd", order=2, n=length(unique(lon$id))*2)+
  f(ids, time, copy="id")
M1 <- inla(formula, data=lon)
summary(M1)

MC_samples <- inla.iidkd.sample(10^4, M1, "id", return.cov = T)
VarCov <- matrix(unlist(MC_samples), nrow=2*2)
VarCovMeans <- matrix(rowMeans(VarCov), 2, 2); round(VarCovMeans, 3)
VarCovSD <- matrix(apply(VarCov, 1, sd), 2, 2);round(VarCovSD, 3)
VarCov025 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.025)), 2, 2) ; round(VarCov025, 3)
VarCov05 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.5)), 2, 2) ; round(VarCov05, 3)
VarCov975 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.975)), 2, 2) ; round(VarCov975, 3)


formula <- Y1 ~ time + f(id, model="iidkd", order=2, n=length(unique(lon$id))*2,
                         hyper=list(theta1=list(param=c(10, 1, 1, 0))))+
  f(ids, time, copy="id")
M1 <- inla(formula, data=lon)
summary(M1)

x <- seq(-10, 10, by=0.1)
plot(x, dnorm(x, 0, 1), type='l')
abline(v=0.2, col='red')

x <- seq(-100, 1000, by=0.1)
plot(x, dnorm(x, 500, 1), type='l')
abline(v=0.2, col='red')

inla.priors.used(M1)

M1 <- inla(formula, data=lon,
           control.fixed=list(mean.intercept=500, prec.intercept=1,
                              mean=0, prec=1))
summary(M1)



N1 <- length(lon$Y1)
N2 <- length(lon$Y2)
NS <- length(unique(lon$id))

data <- list(
  IntY1 <- c(rep(1, N1), rep(NA, N2)),
  IntY2 <- c(rep(NA, N1), rep(1, N2)),
  TimeY1 <- c(lon$time, rep(NA, N2)),
  TimeY2 <- c(rep(NA, N1), lon$time),
  b1_i0Y1 <- c(lon$id, rep(NA, N2)),
  b1_i0Y2 <- c(rep(NA, N1), NS+lon$id),
  b1_i1Y1 <- c(NS+NS+lon$id, rep(NA, N2)),
  b1_i1Y2 <- c(rep(NA, N1), NS+NS+NS+lon$id),
  Yjoint <- list(Y1 = c(lon$Y1, rep(NA, N2)),
                 Y2 = c(rep(NA, N1), lon$Y2))
)

formula <- Yjoint ~ -1 + IntY1 + IntY2 + TimeY1 + TimeY2 +
  f(b1_i0Y1, model="iidkd", order=4, n=NS*4) +
  f(b1_i0Y2, copy="b1_i0Y1") +
  f(b1_i1Y1, TimeY1, copy="b1_i0Y1") +
  f(b1_i1Y2, TimeY2, copy="b1_i0Y1")


M2 <- inla(formula, data=as.data.frame(data), family=c("gaussian", "poisson"),
           control.inla=list(int.strategy="eb"))
summary(M2)

MC_samples <- inla.iidkd.sample(10^4, M2, "b1_i0Y1", return.cov = T)
VarCov <- matrix(unlist(MC_samples), nrow=4*4)
VarCovMeans <- matrix(rowMeans(VarCov), 4, 4); round(VarCovMeans, 3)
VarCovSD <- matrix(apply(VarCov, 1, sd), 4, 4);round(VarCovSD, 3)
VarCov025 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.025)), 4, 4) ; round(VarCov025, 3)
VarCov05 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.5)), 4, 4) ; round(VarCov05, 3)
VarCov975 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.975)), 4, 4) ; round(VarCov975, 3)


# INLAjoint instead of INLA
head(lon)
M3 <- joint(formLong=list(Y1 ~ time + (1+time|id),
                          Y2 ~ time + (1+time|id)), dataLong=lon,
            family=c("gaussian", "poisson"), id="id", timeVar="time", corLong=T)
summary(M3)
summary(M3, sdcor=T) # standard deviation and correlation instead of variance-covariance


# independent markers: corLong=F

# make predictions with observations from 0, 1 and 2 markers


NewData <- lon[lon$id==1,]
P <- predict(M3, NewData, horizon=10, inv.link=T)

library(ggplot2)
theme_set(theme_minimal())
ggplot(P$PredL) + facet_wrap(~Outcome, ncol=2, scales="free") +
  geom_line(aes(x=time, y=quant0.5)) +
  geom_line(aes(x=time, y=quant0.025), linetype="dashed")+
  geom_line(aes(x=time, y=quant0.975), linetype="dashed")+
  theme(legend.position = "none")+
geom_point(data = data.frame(x = NewData$time, y = NewData$Y1, Outcome = "Y1"),
           aes(x = NewData$time, y = NewData$Y1))+
  geom_point(data = data.frame(x = NewData$time, y = NewData$Y2, Outcome = "Y2"),
             aes(x = NewData$time, y = NewData$Y2))



































