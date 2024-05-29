library(spdep)
library(INLA)
library(ggplot2)
library(gridExtra)

load("../Data/brazil.RData")
ls()


table(substr(map_df$code, 1, 1))

### select regions 3 and 4
idselect <- which(substr(map_df$code, 1, 1) %in% 3:4)

map_df <- map_df[idselect, ]
dim(map_df)

n.areas <- nrow(map_df)

## neighbour list
nb <- poly2nb(map_df, queen = FALSE)
nnb <- card(nb)

## as graph
graph <- sparseMatrix(
    i = rep(1:n.areas, nnb),
    j = unlist(nb[nnb>0])
)

### get the data
tbdata <- read.table(
    file = "../Data/tbBrazil.txt",
    header = TRUE)
head(tbdata)

range(tbdata$area[tbdata$time==1][idselect] )
length(unique(tbdata$area[tbdata$time==1][idselect]))

tbdata <- tbdata[tbdata$area %in% tbdata$area[tbdata$time==1][idselect], ]
tbdata$area <- tbdata$area - min(tbdata$area) + 1

range(tbdata$area)
nrow(map_df)

(n.times <- length(unique(tbdata$time)))

tx0 <- 1e5 * sum(tbdata$tb) / sum(tbdata$pop)
tx0

m1 <- tb ~ f(area, model = 'bym2', graph = graph)

r1 <- inla(
    formula = m1,
    family = 'poisson',
    E = pop,
    data = tbdata)

r1$cpu.used

r1$summary.hy[, 1:2]

plot(r1$summary.fitted.values$mean,
     tbdata$tb / tbdata$pop,
     asp = 1,
     bty = 'n')
abline(0:1)

ggplot(map_df) + theme_minimal() +
    geom_sf(aes(fill = r1$summary.random$area$mean[1:n.areas])) +
    scale_fill_distiller(
        palette = "RdBu")

tbdata$area.st <- tbdata$area

m2 <- update(
    m1, . ~ . +
            f(area.st, model = 'besag',
              graph = graph, scale.model = TRUE,
              group = time,
              control.group = list(model = 'ar1'))
)

r2 <- inla(
    formula = m2,
    family = 'poisson',
    E = pop,
    verbose = !TRUE,
    data = tbdata)

r1$cpu.used
r2$cpu.used

grep("Number of constraints", r2$logfile, value = TRUE)
1 + n.times

plot(r2$summary.fitted.values$mean,
     tbdata$tb / tbdata$pop,
     asp = 1,
     bty = 'n')
abline(0:1)

r2$summary.hy

## spatial effect
ggplot(map_df) + theme_minimal() +
    geom_sf(aes(fill = r2$summary.ran$area$mean[1:n.areas])) + 
    scale_fill_distiller(
        palette = "RdBu")

r.st <- matrix(
    r2$summary.random$area.st$mean,
    n.areas,
    n.times)

## map of the 1st year
ggplot(map_df) + theme_minimal() +
    geom_sf(aes(fill = r.st[, 1])) +
    scale_fill_distiller(
        palette = "RdBu")

range(r2$summary.random$area.st$mean)

gg0 <- ggplot(map_df) + theme_minimal() +
    scale_fill_gradient2(
        low = "#49B8F1",
        mid = "white",
        high = "red",
        limits = c(-3, 3))

tsel <- seq(1, n.times, 3)
tsel

lgg <- lapply(tsel, function(i)
    gg0 + 
    geom_sf(aes(fill = r.st[, i]), lwd = 0.1, show.legend = !FALSE) +
    ggtitle(paste("year", 2000+i)))

do.call("grid.arrange", c(lgg, nrow = 2))

