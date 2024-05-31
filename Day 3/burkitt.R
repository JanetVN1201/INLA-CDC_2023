library(sf)
library(splancs)
library(fmesher)
library(INLA)
library(inlabru)
library(ggplot2)
library(gridExtra)

## load Burkitt data from the splancs package
data('burkitt')

ls()

head(burkitt)

par(mar = c(0,0,0,0))
plot(burbdy, bty = 'n', type = "l", asp = 1)
with(burkitt, points(x, y, cex = age/10))

## data as sf
bkt <- st_as_sf(burkitt, coords = c("x", "y"))
bnd <- st_sf(st_sfc(st_polygon(list(burbdy))))

## mesh
mesh <- fm_mesh_2d(
    boundary = bnd, 
    max.edge = c(5, 15),
    cutoff = 2)

## visualize
ggplot() + theme_minimal() + 
    geom_sf(data = bnd) +
    gg(mesh) +
    geom_sf(data = bkt) 

## The spatial model
spde <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(5, 0.01),
    prior.sigma = c(0.3, 0.01)
)

## model components 
components <- ~  
    Intercept(1) +        ## "inlabru" way of doing it
    spatial(geometry,     ## index on the geometry
            model = spde) ## the actual model definition

## likelihood object
lhood <- like(
    geometry ~ .,         ## locations are the outcome
    family = "cp",
    data = bkt,
    domain = list(geometry = mesh),
    samplers = bnd)

## fit the model
fit <- bru(
    components,
    lhood)

fit$cpu.used
fit$bru_timings

## summaries
fit$summary.hyperpar

## prepare for visualization
post.range <- spde.posterior(fit, name = "spatial", what = "range")
post.sigma <- spde.posterior(fit, name = "spatial", what = "variance")
post.corr <- spde.posterior(fit, name = "spatial", what = "matern.correlation")

## visualize
grid.arrange(
    plot(post.range),
    plot(post.sigma),
    plot(post.corr),
    nrow = 1)

fit$summary.fixed

nrow(burkitt) / st_area(bnd)
exp(-6.61)

## the bounding box
st_bbox(bnd)
bb <- matrix(st_bbox(bnd), 2)
bb

apply(bb, 1, diff)


## setup a grid
grid <- fm_pixels(
    mesh = mesh,
    dims = c(95, 182),
    mask = bnd,
    format = "sf")

str(grid)

## inlabru::predict()
##  drawn (Monte Carlo: independent) samples
## from the model parameters posterior fitted by INLA
## and compute functions from these
pred <- predict(
    fit, 
    grid, 
    ~ data.frame(
    lambda = exp(Intercept + spatial),
    loglambda = Intercept + spatial
    )
)

str(pred)

## visualize
ggplot() + theme_minimal() +
    geom_sf(data = pred$loglambda,            
            aes(color = mean)) +
    scale_color_distiller(
        palette = 'RdBu'
    ) +
    geom_sf(data = bkt)

## some interpretation 

hist(pred$lambda$mean)

summary(pred$lambda$mean)

nrow(burkitt) / st_area(bnd)

exp(-6.61)

## extra

## Explore \lambda() 

## setup integration points
ipts <- fm_int(
    domain = mesh,
    samplers = bnd)

str(ipts)

ggplot(ipts) +
    theme_minimal() +
    geom_sf(data=bnd, fill = gray(0.9)) + 
    geom_sf(aes(color = weight)) 

## as sf
ipts.sf <- st_sf(
    data.frame(
        st_drop_geometry(ipts), 
        geometry = st_geometry(ipts)
    )
)

str(ipts.sf)

## Lambda samples
iLambda <- predict(
  fit,
  ipts.sf,
  ~ weight * exp(spatial + Intercept)
)

str(iLambda)

c(sum(iLambda$mean),
  sum(iLambda$q0.025),
  sum(iLambda$q0.975))

### setup sub-areas
st_bbox(bnd)
bb
apply(bb, 1, diff)

hx <- 46

subAreas.grid <- GridTopology(
    bb[, 1] + hx/2,
    cellsize = c(hx, hx),
    cells.dim = c(2, 4)
)
subAreas.grid

Gridb <- SpatialGrid(subAreas.grid)

plot(Gridb)
points(st_coordinates(st_geometry(ipts)), pch = 4, cex = 0.5)

gridSP <- as.SpatialPolygons.GridTopology(subAreas.grid)
ipts.b <- fm_int(
    domain = mesh,
    samplers = gridSP
)

str(ipts.b)
table(ipts.b$.block)

plot(ipts.b, col = ipts.b$.block)

## Lambda samples with the new setup
bLambda <- predict(
  fit,
  ipts.b,
  ~ weight * exp(spatial + Intercept)
)

str(bLambda)

tapply(bLambda$mean, bLambda$.block, sum)
tapply(bLambda$q0.025, bLambda$.block, sum)
tapply(bLambda$q0.975, bLambda$.block, sum)

obs.counts <- st_within(
    bkt,
    st_as_sf(gridSP)
)

table(unlist(obs.counts))
tapply(bLambda$mean, bLambda$.block, sum)

## see a complete example at
## https://inlabru-org.github.io/inlabru/articles/2d_lgcp.html
