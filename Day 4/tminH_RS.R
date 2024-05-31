
library(sp)
library(sf)
library(ggplot2)
library(INLA)
library(inlabru)
library(INLAspacetime)

saveRDS(infodata, "infodataRS.rds")
infodata <- readRDS("infodataRS.rds")

## visualize
ggplot() +
    theme_minimal() + 
    geom_sf(
        data = st_as_sf(infodata),
        aes(color = precApr),
        cex = 7
    ) +
    scale_color_distiller(
        na.value = 'transparent',
        palette = 'OrRd',
        direction = 1
    )

## spacetime data: 
stdata <- readRDS(
    file = "../Data/tempMinH_RS202404.rds"
)

dim(stdata)
stdata[1:2, 1:5, 1:2]

(n.h <- dim(stdata)[1])
(n.d <- dim(stdata)[2])
(n.s <- dim(stdata)[3])

plot(infodata)
stlines(stdata[1, , ], infodata, pch = 8)

plot(infodata)
stlines(stdata[, 1, ], infodata, pch = 8)

dim(stdata[, 1, ])

plot(infodata)
stlines(rbind(stdata[, 1, ],
              stdata[, 2, ],
              stdata[, 3, ]),
        infodata, pch = 8)

n.t <- n.h ## 1 day

ldata <- data.frame(
    xloc = rep(coordinates(infodata)[, 1], each = n.t),
    yloc = rep(coordinates(infodata)[, 2], each = n.t),
    time = rep(1:n.t, n.s),
    elev = rep(infodata$elev, each = n.t),
    tminH = as.vector(stdata[, 30, ])
)
head(ldata)

### prior for sigma_e
ctrlf <- list(
    hyper = list(
        prec = list(
            prior = "pcprec",
            param = c(1, 0.1)
        )
    )
)

### likelihood object
datalike <- like(
  formula = tminH ~ ., 
  family = "gaussian",
  control.family = ctrlf, 
  data = ldata
)

## temporal mesh
ht <- 2
tmesh <- inla.mesh.1d(
    loc = seq(1, n.t+ht, ht)
)
tmesh$n
tmesh$mid
n.t

## spatial mesh
smesh <- inla.mesh.2d(
  loc = coordinates(infodata),
  max.edge = c(0.15, 0.5) * 5,
  offset = c(0.5, 1.5)
)
smesh$n

## visualize
ggplot() +
    theme_minimal() + 
    gg(smesh) +
    geom_sf(
        data = st_as_sf(infodata),
        aes(color = elev),
        cex = 7
    ) +
    scale_color_distiller(
        palette = 'RdBu'
    )

sspde <- inla.spde2.pcmatern(
    mesh = smesh,
    prior.range = c(5, NA), ## range fixed to 5
    prior.sigma = c(1, 0.05)
)

stmodel1 <- stModel.define(
    smesh = smesh, ## spatial mesh
    tmesh = tmesh, ## temporal mesh
    model = '102', ## model, see the paper
    control.priors = list(
        prs = c(0.5, 0.1), ## P(spatial range < 1) = 0.1
        prt = c(10, 0.0), ## temporal range fixed to 10
        psigma = c(1, 0.1) ## P(sigma > 1) = 0.1
    )
)

stmodel2 <- stModel.define(
    smesh = smesh, ## spatial mesh
    tmesh = tmesh, ## temporal mesh
    model = '121', ## model, see the paper
    control.priors = list(
        prs = c(0.5, 0.1), ## P(spatial range < 1) = 0.1
        prt = c(20, 0.0), ## temporal range fixed to 20
        psigma = c(1, 0.1) ## P(sigma > 1) = 0.1
    )
)

linpred1 <- ~ 1 + elev +
    spatial(cbind(xloc, yloc),
            model = sspde) + 
    field(list(space = cbind(xloc, yloc), 
               time = time),
          model = stmodel1)

fit1 <- bru(
    components = linpred1,
    datalike,
    options = list(
        control.inla = list(
            int.strategy = "eb"
        ),
        control.predictor = list(link = 1),
        verbose = !TRUE)
)

linpred2 <- ~ 1 + elev +
    spatial(cbind(xloc, yloc),
            model = sspde) + 
    field(list(space = cbind(xloc, yloc), 
               time = time),
          model = stmodel2)

fit2 <- bru(
    components = linpred2,
    datalike,
    options = list(
        control.inla = list(
                int.strategy = "eb"
        ),
        control.predictor = list(link = 1),
        verbose = !TRUE)
)

fit1$misc$nfunc
fit2$misc$nfunc

grep("nnz", fit1$logfile, value = TRUE)
grep("nnz", fit2$logfile, value = TRUE)

rbind(fit1$cpu.used,
      fit2$cpu.used)

rbind(
    fit1$summary.fixed, 
    fit2$summary.fixed
)

fit1$summary.hyperpar[, 1:2]
exp(fit1$summary.hyperpar$mean[3:4])

fit2$summary.hyperpar[, 1:2]
exp(fit2$summary.hyperpar$mean[3:4])

r1 <- fit1$summary.fitted.values[1:(n.s*n.t), ]
r2 <- fit2$summary.fitted.values[1:(n.s*n.t), ]

summary(r1$mean)
summary(r2$mean)

plot(r1$mean, ldata$tminH, asp = 1, bty = "n")
points(r2$mean, ldata$tminH, col = 2, pch = 8)
abline(0:1, lty = 2)

d.all <- data.frame(
    station = rep(1:n.s, each = n.t),
    ldata,
    r1 = r1,
    r2 = r2
)

ggplot(d.all, aes(x = time)) +
    theme_minimal() +
    geom_point(
        aes(y = tminH)
    ) +
    geom_ribbon(
        aes(ymin = r2.0.025quant,
            ymax = r2.0.975quant),
            fill = gray(0.9, 0.5)
    ) +
    geom_line(
        aes(y = r2.mean)
    ) +
    facet_wrap(~station) 


r2mat <- matrix(d.all$r2.mean, n.t)

par(mfrow = c(5,9), mar = c(1.5, 1.5, 0.1, 0.1), mgp = c(1,0.5,0))
for(i in 1:n.s) {
    plot(stdata[, 30, i], bty = "n", 
         ylim = range(r2mat[, i], na.rm = TRUE),
         xlab = "", ylab = ""
         )
    lines(r2mat[, i])
}

par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
plot(infodata)
stlines(r2mat, infodata)

par(mfrow = c(1, 2))
hist(fit1$summary.random$field$mean) # do not worry!
hist(fit2$summary.random$field$mean) # do not worry!
