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
    prior.range = c(10, 0.01),
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

## setup a grid
grid <- st_make_grid(
    x = bnd, 
    cellsize = 2
)

ggplot() +
    geom_sf(data = grid)

pred <- predict(
    fit,
    newdata = grid,
    ~ exp(Intercpt + spatial),
    n.samples = 1000
)

bb <- matrix(st_bbox(bnd), 2)
bb

apply(bb, 1, diff)

grid <- fm_pixels(
    mesh = mesh,
    dims = c(95, 182),
    mask = bnd,
    format = "sf")

pred <- predict(
    fit, 
    grid, 
    ~ data.frame(
    lambda = exp(Intercept + spatial),
    loglambda = Intercept + spatial
    )
)

str(pred)

ggplot() + theme_minimal() +
    geom_sf(data = pred$loglambda,            
            aes(color = mean)) +
    scale_color_distiller(
        palette = 'RdBu'
    )

hist(pred$lambda$mean)
summary(pred$lambda$mean)
nrow(burkitt) / st_area(bnd)
exp(-6.61)
