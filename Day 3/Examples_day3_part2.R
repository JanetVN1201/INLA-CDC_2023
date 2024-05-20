library(INLAjoint) 
library(JM)

data(pbc2) # dataset

pbc2time05 <- pbc2[c(2,which(diff(as.integer(pbc2$id))==1)+2),]
pbc2time05 <- pbc2time05[-which(diff(as.integer(pbc2time05$id))==0),]
pbc2time05$year <- 0.5

pbc2_ni <- rbind(pbc2[c(1,which(diff(as.integer(pbc2$id))==1)+1),], pbc2time05)
pbc2_ni <- pbc2_ni[order(pbc2_ni$id),]
pbc2_2 <- pbc2[, c("id", "year", "serBilir", "drug")]
print(head(pbc2_2, 10), row.names=F)

M1 <- joint(formLong = serBilir ~ year + drug  +  (1 + year|id),
            dataLong = pbc2, id = "id", timeVar = "year",
            family = "lognormal")
summary(M1)
summary(M1, sdcor=TRUE)
inla.priors.used(M1)
plot(M1)
plot(M1, priors=TRUE, sdcor=TRUE)$Covariances



M2 <- joint(formLong = serBilir ~ year + drug  +  (1 + year|id),
            dataLong = pbc2, id = "id", timeVar = "year",
            family = "lognormal", control=list(priorRandom=list(r=100, R=1)))

M3 <- joint(formLong = serBilir ~ year + drug  +  (1 + year|id),
            dataLong = pbc2_ni, id = "id", timeVar = "year",
            family = "lognormal", control=list(priorRandom=list(r=10, R=1)))

M4 <- joint(formLong = serBilir ~ year + drug  +  (1 + year|id),
            dataLong = pbc2_ni, id = "id", timeVar = "year",
            family = "lognormal", control=list(priorRandom=list(r=100, R=1)))

plot(M1, priors=TRUE, sdcor=TRUE)$Covariances
plot(M2, priors=TRUE, sdcor=TRUE)$Covariances
plot(M3, priors=TRUE, sdcor=TRUE)$Covariances
plot(M4, priors=TRUE, sdcor=TRUE)$Covariances



# This chunk makes the plot for the correlation parameter only:
# A <- plot(M1, priors=TRUE, sdcor=TRUE)$Covariances
# B <- plot(M2, priors=TRUE, sdcor=TRUE)$Covariances
# C <- plot(M3, priors=TRUE, sdcor=TRUE)$Covariances
# D <- plot(M4, priors=TRUE, sdcor=TRUE)$Covariances
#
# A$L1$data <- A$L1$data[which(A$L1$data$Effect=="year_L1:Intercept_L1"), ]
# A$L1$data$Effect <- "Identifiable model, prior 1"
# B$L1$data <- B$L1$data[which(B$L1$data$Effect=="year_L1:Intercept_L1"), ]
# B$L1$data$Effect <- "Identifiable model, prior 2"
# C$L1$data <- C$L1$data[which(C$L1$data$Effect=="year_L1:Intercept_L1"), ]
# C$L1$data$Effect <- "Non-identifiable model, prior 1"
# D$L1$data <- D$L1$data[which(D$L1$data$Effect=="year_L1:Intercept_L1"), ]
# D$L1$data$Effect <- "Non-identifiable model, prior 2"
#
# E <- rbind(A$L1$data, B$L1$data, C$L1$data, D$L1$data)
#
# CorPlots <- ggplot(E, aes(x=x,y=y,group=group)) +
#   xlab('') +
#   ylab('Density') +
#   geom_line(aes(color=group, linetype=group)) +
#   facet_wrap(~Effect, scales='free') + xlim(-1,1) + ylim(0,6.1)
#
# CorPlots







