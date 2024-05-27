### Examples CDC training - day 3 - part 3

library(INLA)
library(INLAjoint)
# Example: Multivariate longitudinal

set.seed(1)
# data generation - two longitudinal markers
nsujet=300 # number of individuals

# Y1 (continuous)
b1_0=2 # intercept
b1_1=-0.5 # slope
b1_e=0.1 # residual error

# Y2 (counts)
b2_0=4 # intercept
b2_1=-0.5 # slope

gap=1 # gap between measurements
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # total number of measurements
id<-rep(1:nsujet, each=nmesindiv) # individual id

s1 <- 0.4 # random intercept Y1 SD
s2 <- 0.3 # random slope Y1 SD
s3 <- 0.5 # random intercept Y2 SD
s4 <- 0.4 # random slope Y1 SD
c12 <- 0.8 # correlations
c13 <- 0.8
c14 <- 0.7
c23 <- 0.8
c24 <- 0.8
c34 <- 0.4

cov_12 <- s1*s2*c12 # covariances
cov_13 <- s1*s3*c13
cov_14 <- s1*s4*c14
cov_23 <- s2*s3*c23
cov_24 <- s2*s4*c24
cov_34 <- s3*s4*c34

# Random effects variance-covariance matrix
Sigma=matrix(c(s1^2,cov_12,cov_13,cov_14,
               cov_12,s2^2,cov_23,cov_24,
               cov_13,cov_23,s3^2,cov_34,
               cov_14,cov_24,cov_34,s4^2),ncol=4,nrow=4)

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
head(lon, 20)













# first fit a mixed effects model for the first marker
ids <- length(unique(lon$id))+lon$id

formula <- Y1 ~ time + f(id, model="iid2d", n=length(unique(lon$id))*2)+
  f(ids, time, copy="id")
M1 <- inla(formula, data=lon)
summary(M1)

# alternative parametrization
formula <- Y1 ~ time + f(id, model="iidkd", order=2, n=length(unique(lon$id))*2,
                         hyper = list(theta1 = list(param = c(10, 1, 1, 0))))+
  f(ids, time, copy="id")
M1 <- inla(formula, data=lon)
summary(M1)

MC_samples <- inla.iidkd.sample(10^4, M1, "id", return.cov = T)
VarCov <- matrix(unlist(MC_samples), nrow=2*2)
VarCovMeans <- matrix(rowMeans(VarCov), 2, 2); round(VarCovMeans, 3)
# VarCovSD <- matrix(apply(VarCov, 1, sd), 2, 2);round(VarCovSD, 3)
# VarCov025 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.025)), 2, 2) ; round(VarCov025, 3)
# VarCov05 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.5)), 2, 2) ; round(VarCov05, 3)
# VarCov975 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.975)), 2, 2) ; round(VarCov975, 3)











# bivariate model
# structure for multi-outcome models:
# L1 .
# L1 .
# L1 .
# .  L2
# .  L2
# .  L2
# .  L2
# => size = N1 + N2

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
  f(b1_i0Y1, model="iidkd", order=4, n=NS*4,
    hyper = list(theta1 = list(param = c(10, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)))) +
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

Sigma[c(1,3,2,4), c(1,3,2,4)]











# INLAjoint instead of INLA
head(lon)
M3 <- joint(formLong=list(Y1 ~ time + (1+time|id),
                          Y2 ~ time + (1+time|id)), dataLong=lon,
            family=c("gaussian", "poisson"), id="id", timeVar="time", corLong=T)
summary(M3)
summary(M3, sdcor=T) # standard deviation and correlation instead of variance-covariance


# independent markers: corLong=F












# predictions
NewData0 <- lon[1,]
NewData0$Y1 <- NewData0$Y2 <- NA
NewData0
P0 <- predict(M3, NewData0, horizon=10, inv.link=T)

library(ggplot2)
theme_set(theme_minimal())
ggplot(P0$PredL) + facet_wrap(~Outcome, ncol=2, scales="free") +
  geom_line(aes(x=time, y=quant0.5)) +
  geom_line(aes(x=time, y=quant0.025), linetype="dashed")+
  geom_line(aes(x=time, y=quant0.975), linetype="dashed")+
  theme(legend.position = "none")

# 1 marker observed
NewData1 <- lon[lon$id==1,]
NewData1$Y1 <- NA
NewData1
P1 <- predict(M3, NewData1, horizon=10, inv.link=T)

library(ggplot2)
theme_set(theme_minimal())
ggplot(P1$PredL) + facet_wrap(~Outcome, ncol=2, scales="free") +
  geom_line(aes(x=time, y=quant0.5)) +
  geom_line(aes(x=time, y=quant0.025), linetype="dashed")+
  geom_line(aes(x=time, y=quant0.975), linetype="dashed")+
  theme(legend.position = "none")+
  geom_point(data = data.frame(x = NewData1$time, y = NewData1$Y2, Outcome = "Y2"),
             aes(x = NewData1$time, y = NewData1$Y2))


# 2 markers observed
NewData2 <- lon[lon$id==1,]
NewData2
P2 <- predict(M3, NewData2, horizon=10, inv.link=T)

library(ggplot2)
theme_set(theme_minimal())
ggplot(P2$PredL) + facet_wrap(~Outcome, ncol=2, scales="free") +
  geom_line(aes(x=time, y=quant0.5)) +
  geom_line(aes(x=time, y=quant0.025), linetype="dashed")+
  geom_line(aes(x=time, y=quant0.975), linetype="dashed")+
  theme(legend.position = "none")+
geom_point(data = data.frame(x = NewData2$time, y = NewData2$Y1, Outcome = "Y1"),
           aes(x = NewData2$time, y = NewData2$Y1))+
  geom_point(data = data.frame(x = NewData2$time, y = NewData2$Y2, Outcome = "Y2"),
             aes(x = NewData2$time, y = NewData2$Y2))

















