## --------------------------------------------------------------------------------------------------------------
#| label: packages
options("rgdal_show_exportToProj4_warnings"="none")
options(warn = -1)
# Robert Hijmans raster and vector data
library(terra, warn.conflicts=FALSE, quiet = TRUE) # replaces `raster`
# still needed to convert to `sp`
library(raster, warn.conflicts=FALSE, quiet = TRUE) 
# Pebesma et al. spation-temporal data
# Simple Features
library(sf, warn.conflicts=FALSE, quiet = TRUE)          
# `sp` spatial classes -- still needed for conversions
library(sp, warn.conflicts=FALSE, quiet = TRUE)
# variogram modelling
library(gstat, warn.conflicts=FALSE, quiet = TRUE)
# Co-occurrence vectors
library(motif, warn.conflicts=FALSE, quiet = TRUE) 
# multivariate distance metrics
library(philentropy, warn.conflicts=FALSE, quiet = TRUE) 
# compare polygon map spatial structures, V measure
library(sabre, warn.conflicts=FALSE, quiet = TRUE)   
# FRAGSTATS-style metrics
library(landscapemetrics, warn.conflicts=FALSE, quiet = TRUE) 
# utility functions for raster* landscape objects)
library(landscapetools, warn.conflicts=FALSE, quiet = TRUE)
# ggplot graphics
library(ggplot2, warn.conflicts=FALSE, quiet = TRUE)
# multiple graphics in one plot
library(gridExtra, warn.conflicts=FALSE, quiet = TRUE)
# ggplot with terra SpatRaster
library(tidyterra, warn.conflicts=FALSE, quiet = TRUE)
# data wrangling
library(dplyr, warn.conflicts=FALSE, quiet = TRUE)
# colour palettes for graphics
library(RColorBrewer, warn.conflicts=FALSE, quiet = TRUE)
# to access NRCS databases
library(soilDB, warn.conflicts=FALSE, quiet = TRUE)
# supercells
# install.packages("supercells", repos = "https://nowosad.r-universe.dev")
library(supercells)
# compare two rasters directly -- in development
# devtools::install_github("Nowosad/spquery")
# library(spquery)
# analyze cross-classification matrices
library(diffeR) 


## ----dirs------------------------------------------------------------------------------------------------------
#| label: set-base-dir
file.dir <- path.expand("~/ds_reference/Compare_DSM/")


## ----import.maps-----------------------------------------------------------------------------------------------
#| label: import
file.dir <- path.expand("~/ds_reference/Compare_DSM/")
(gn <- rast(paste0(file.dir, "gNATSGO/lat4243_lon-77-76/ph1to1h2o_r_05_250.tif")))
(sg <- rast(paste0(file.dir, "SoilGrids250/lat4243_lon-77-76/phh2o_0-5cm_mean.tif")))


## ----mult.gn---------------------------------------------------------------------------------------------------
values(sg) <- values(sg)/10


## --------------------------------------------------------------------------------------------------------------
#| label: fig-sg-gn
#| fig-cap: "pH, 0-5 cm"
#| warning: false
range.sg.gn <- range(range(values(sg), na.rm = TRUE), 
                      range(values(gn), na.rm = TRUE))
par(mfrow=c(1,2))
terra::plot(sg, main = "SoilGrids v2.0", 
     range = range.sg.gn, col=(sp::bpy.colors(50)))
terra::plot(gn, main = "gNATSGO", 
     range = range.sg.gn, col=(sp::bpy.colors(50)))
par(mfrow=c(1,1))


## ----crop.to.quarter-------------------------------------------------------------------------------------------
#| label: crop
test.tile.size <- 0.25  # degrees
test.tile.x.offset <- 0.2 # lrc west from right edge
test.tile.y.offset <- 0.3  # lrc north from bottom edge
ext.crop <- round(as.vector(ext(sg)),2) # line up to .00 decimal degrees
ext.crop["xmax"] <- ext.crop["xmax"] - test.tile.x.offset
ext.crop["xmin"] <- ext.crop["xmax"] - test.tile.size
ext.crop["ymin"] <- ext.crop["ymin"] + test.tile.y.offset
ext.crop["ymax"] <- ext.crop["ymin"] + test.tile.size
ext(ext.crop)
gn <- crop(gn, ext(ext.crop));sg <- crop(sg, ext(ext.crop))
ext(gn)
ext(sg)


## --------------------------------------------------------------------------------------------------------------
#| label: crop.2
(range.sg <- range(values(sg), na.rm = TRUE))
(range.gn <- range(values(gn), na.rm = TRUE))
range.sg.gn <- range(range(values(sg), na.rm = TRUE), 
                     range(values(gn), na.rm = TRUE))
par(mfrow=c(1,2))
plot(sg, main = "SoilGrids v2.0", 
     range = range.sg.gn, col=(sp::bpy.colors(50)))
plot(gn, main = "gNATSGO", 
     range = range.sg.gn, col=(sp::bpy.colors(50)))
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: get.utm
# a function to find the correct UTM zome
long2UTM <- function(long) { (floor((long + 180)/6) %% 60) + 1 }
# find the zone from the central meridian
utm.zone <- long2UTM(st_bbox(sg)$xmin + 
                       0.5*(st_bbox(sg)$xmax - st_bbox(sg)$xmin))
cat(paste("UTM Zone", utm.zone))
epsg.utm <- paste0("epsg:326", utm.zone)
cat(paste("CRS code:", epsg.utm))


## --------------------------------------------------------------------------------------------------------------
#| label: resample
st_bbox(gn)
res(gn)
gn.utm <- terra::project(gn, epsg.utm, 
                         res = c(250, 250), method = "bilinear")
st_bbox(gn.utm)
res(gn.utm)
st_bbox(sg)
sg.utm <- terra::project(sg, epsg.utm, 
                         res = c(250, 250), method = "bilinear")
st_bbox(sg.utm)


## ----make.extents.identical------------------------------------------------------------------------------------
ext(gn.utm)
ext(sg.utm)
sg.utm <- resample(sg.utm, gn.utm)
ext(sg.utm)


## --------------------------------------------------------------------------------------------------------------
#| label: mask
sg.utm <- mask(sg.utm, gn.utm)
# SoilGrids now has some `NA` added from gNATSGO
gn.utm <- mask(gn.utm, sg.utm)
# The added `NA` are already in gNATSGO, now it gets `NA` originally on SoilGrids


## --------------------------------------------------------------------------------------------------------------
#| label: fig-sg-gn-utm
#| fig-cap: "pH, 0-5 cm"
#| warning: false
par(mfrow=c(1,2))
plot(sg.utm, main = "SoilGrids v2.0", 
     range = range.sg.gn, col=(sp::bpy.colors(50)))
plot(gn.utm, main = "gNATSGO", 
     range = range.sg.gn, col=(sp::bpy.colors(50)))
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: download-ssurgo
# get the polygons with their key
system.time(
  mu.poly <- SDA_spatialQuery(gn, 
                            what = "mupolygon", 
                            db = "SSURGO", 
                            geomIntersection = TRUE)
)
class(mu.poly)
st_crs(mu.poly)$proj4string
summary(mu.poly)
head(mu.poly)
# plot with area in acres
plot(mu.poly, y = "area_ac",
     type = "continuous",
     main = "SSURGO map units, area in acres")


## --------------------------------------------------------------------------------------------------------------
mu.key <- SDA_spatialQuery(gn.utm, 
                           what = "mukey", 
                           db = "SSURGO", 
                           geomIntersection = TRUE)
head(mu.key)


## --------------------------------------------------------------------------------------------------------------
# format the list of map units for SQL
IS <- soilDB::format_SQL_in_statement(mu.poly$mukey)
# query string -- all components
ws <- sprintf("mukey IN %s", IS)
# format the SQL query
query <- paste("SELECT * FROM mapunit WHERE", ws)
# and run it
mu.info <- SDA_query(query)
dim(mu.info)
names(mu.info)
series.name <- "Ovid" # look for a map unit by name
length(ix <- which(substr(mu.info$muname, 1, 
                          nchar(series.name)) == series.name))
mu.info[ix, ]


## --------------------------------------------------------------------------------------------------------------
st_crs(mu.poly)$proj4string
mu.poly <- terra::project(mu.poly, epsg.utm)
st_crs(mu.poly)$proj4string


## --------------------------------------------------------------------------------------------------------------
gn.sp <- as(raster(gn.utm), "SpatialPointsDataFrame")
gn.sf <- st_as_sf(gn.sp)
names(gn.sf)
sg.sp <- as(raster(sg.utm), "SpatialPointsDataFrame")
sg.sf <- st_as_sf(sg.sp)
names(sg.sf)


## --------------------------------------------------------------------------------------------------------------
#| label: variogram_parameters
range.init <- 1000  # estimated range, m 
cutoff.init <- range.init*5  # cutoff for empirical variogram, m
width.init <- 250   # bin width


## --------------------------------------------------------------------------------------------------------------
#| label: variogram_empirical
#| fig-width: 8
#| fig-height: 4
v.sg <- variogram(phh2o_0.5cm_mean ~ 1, loc = sg.sf, 
                  cutoff=cutoff.init, width=width.init)
#
v.gn <- gstat::variogram(ph1to1h2o_r ~ 1, loc = gn.sf, 
                         cutoff=cutoff.init, width=width.init)
ylim.v <- max(v.gn$gamma, v.sg$gamma)
p1 <- plot(v.sg, ylim = c(0, ylim.v), main = "SoilGrids v2.0")
p2 <- plot(v.gn, ylim = c(0, ylim.v), main = "gNATSGO")
grid.arrange(p1, p2, nrow = 1)


## --------------------------------------------------------------------------------------------------------------
#| label: variogram_model
vm.gn <- vgm(psill = 0.8*max(v.gn$gamma), 
             model = "Exp", 
             range = range.init, 
             nugget = 0)
print(vmf.gn <- fit.variogram(v.gn, model=vm.gn))
vm.sg <- vgm(0.8*max(v.sg$gamma), "Exp", range.init, 0)
print(vmf.sg <- fit.variogram(v.sg, model=vm.sg))


## --------------------------------------------------------------------------------------------------------------
#| label: fitted_variogram_model
#| fig-width: 8
#| fig-height: 4
p1 <- plot(v.sg, model=vmf.sg, ylim = c(0, ylim.v), main = "SoilGrids v2.0", 
     xlab = "separation m", ylab = expression(paste(Delta, plain(pH)^2)))
p2 <- plot(v.gn, model=vmf.gn, ylim = c(0, ylim.v), main = "gNATSGO", 
     xlab = "separation m", ylab = expression(paste(Delta, plain(pH)^2)))
grid.arrange(p1, p2, nrow = 1)


## --------------------------------------------------------------------------------------------------------------
#| label: wts-matrix
# for a 5x5 matrix
# there must be a more elegant way to do this!
(vl <- variogramLine(vm.sg,   
                     dist_vector = c(0,
                                     250, 250*sqrt(2), 
                                     500, 250*sqrt(5), 
                                     500*sqrt(2))))
(w.r <- 1- vl$gamma)  # relative weights
(w.m <- matrix(c(w.r[6], w.r[5], w.r[4], w.r[5], w.r[6],
                w.r[5], w.r[3], w.r[2], w.r[3], w.r[5],
                w.r[4], w.r[2], w.r[1], w.r[2], w.r[4],
                w.r[5], w.r[3], w.r[2], w.r[3], w.r[5],
                w.r[6], w.r[5], w.r[4], w.r[5], w.r[6]), 
                nrow = 5, ncol = 5))


## --------------------------------------------------------------------------------------------------------------
#| label: moving-window-5
sg.utm.autocor <- terra::autocor(sg.utm, w=w.m, 
                                 method="moran", global = FALSE)
gn.utm.autocor <- terra::autocor(gn.utm, w=w.m, 
                                 method="moran", global = FALSE)
(range.sg.autocor <- range(values(sg.utm.autocor), na.rm = TRUE))
(range.gn.autocor <- range(values(gn.utm.autocor), na.rm = TRUE))
range.autocor <- range(range.sg.autocor, range.gn.autocor)
# hcl.pals(type = "diverging")
par(mfrow=c(1,2))
terra::plot(sg.utm.autocor, main = "SG250, Moran's I, 5x5", 
            range = range.autocor, col = rev(hcl.colors(32, palette = "RdYlGn")))
terra::plot(gn.utm.autocor, main = "gNATSGO, Moran's I, 5x5", 
            range = range.autocor, col = rev(hcl.colors(32, palette = "RdYlGn")))
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: global-autocor
terra::autocor(sg.utm, w=w.m, method="moran", global = TRUE)
terra::autocor(gn.utm, w=w.m, method="moran", global = TRUE)


## ----cut.matrix------------------------------------------------------------------------------------------------
range(values(sg.utm, na.rm = TRUE))
sg.quant <- cut(values(sg.utm), breaks = 16, labels = 0:15, include.lowest = TRUE)
table(sg.quant)
# show the breakpoints
levels(cut(values(sg.utm), breaks = 16, include.lowest = TRUE))
sg.utm.quant <- sg.utm; values(sg.utm.quant) <- sg.quant
plot(sg.utm.quant, col = rainbow(16), main = "pH, 16 levels")


## ----fig.caption = "Example GLCM"------------------------------------------------------------------------------
require(GLCMTextures)  # v0.4.1
test.quant <- cut(values(sg.utm), breaks = 16, 
                  labels = 0:15, include.lowest = TRUE)
table(test.quant)
(l.16 <- levels(cut(values(sg.utm), breaks = 16, include.lowest = TRUE)))
test.rast <- sg.utm; values(test.rast) <- test.quant
dim(test.rast); ext(test.rast)
(xy <- xyFromCell(test.rast, cellFromRowCol(test.rast, 50:55, 50:55)))
w.sg <- crop(sg.utm, xy)
w.sg.16 <- crop(test.rast, xy)
par(mfrow = c(1,2))
plot(w.sg, main = "pH")
plot(w.sg.16, main = "grey level", col = grey.colors(16))
par(mfrow = c(1,1))
test.matrix <- as.matrix(w.sg.16, wide = TRUE)
make_glcm(test.matrix,  
          n_levels = 8, shift = c(1, 0), # shift one cell to the right
          normalize = FALSE )


## ----glcm------------------------------------------------------------------------------------------------------
require(glcm)
# convert to the older `raster` format
sg.utm.raster <- raster(sg.utm)
gn.utm.raster <- raster(gn.utm)


## ----measures--------------------------------------------------------------------------------------------------
stat.list <- c("mean","variance","homogeneity","contrast",
               "entropy","dissimilarity","second_moment",
               "correlation")
glcm.sg <- rast(glcm(sg.utm.raster,
                   window = c(5, 5),
                   n_grey = 32, # number of levels in the GLCM
                   shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)), # all directions
                   na_opt = "ignore",
                   statistics = stat.list))
# gNATSGO is not perfectly square
glcm.gn <- rast(glcm(gn.utm.raster,
                   window = c(5, 5),
                   n_grey = 32, # number of levels in the GLCM
                   shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)), # all directions
                   na_opt = "ignore",
                   statistics = stat.list))
class(glcm.sg)
summary(glcm.sg)
summary(glcm.gn)
summary(glcm.sg-glcm.gn)
plot(glcm.sg)
plot(glcm.gn)


## ----fig.caption = "GLCM mean, 5 x 5 grid cells"---------------------------------------------------------------
zlim <- range(range(values(glcm.sg[["glcm_mean"]]), na.rm = TRUE), 
                     range(values(glcm.gn[["glcm_mean"]]), na.rm = TRUE))
par(mfrow=c(1,2))
plot(glcm.sg[["glcm_mean"]], main = "SoilGrids v2.0", 
     range = zlim, col=(sp::bpy.colors(32)))
plot(glcm.gn[["glcm_mean"]], main = "gNATSGO", 
     range = zlim, col=(sp::bpy.colors(32)))
par(mfrow=c(1,1))


## ----fig.caption = "GLCM variance, 5 x 5 grid cells"-----------------------------------------------------------
zlim <- range(range(values(glcm.sg[["glcm_variance"]]), na.rm = TRUE), 
                     range(values(glcm.gn[["glcm_variance"]]), na.rm = TRUE))
par(mfrow=c(1,2))
plot(glcm.sg[["glcm_variance"]], main = "SoilGrids v2.0", 
     range = zlim, col=(cm.colors(32)))
plot(glcm.gn[["glcm_variance"]], main = "gNATSGO", 
     range = zlim, col=(cm.colors(32)))
par(mfrow=c(1,1))


## ----fig.caption = "GLCM contrast, 5 x 5 grid cells"-----------------------------------------------------------
zlim <- range(range(values(glcm.sg[["glcm_contrast"]]), na.rm = TRUE), 
                     range(values(glcm.gn[["glcm_contrast"]]), na.rm = TRUE))
par(mfrow=c(1,2))
plot(glcm.sg[["glcm_contrast"]], main = "SoilGrids v2.0", 
     range = zlim, col=(topo.colors(32)))
plot(glcm.gn[["glcm_contrast"]], main = "gNATSGO", 
     range = zlim, col=(topo.colors(32)))
par(mfrow=c(1,1))


## ----fig.caption = "GLCM dissimilarity, 5 x 5 grid cells"------------------------------------------------------
zlim <- range(range(values(glcm.sg[["glcm_dissimilarity"]]), na.rm = TRUE), 
                     range(values(glcm.gn[["glcm_dissimilarity"]]), na.rm = TRUE))
par(mfrow=c(1,2))
plot(glcm.sg[["glcm_dissimilarity"]], main = "SoilGrids v2.0", 
     range = zlim, col=(topo.colors(32)))
plot(glcm.gn[["glcm_dissimilarity"]], main = "gNATSGO", 
     range = zlim, col=(topo.colors(32)))
par(mfrow=c(1,1))


## ----fig.caption = "GLCM entropy, 5 x 5 grid cells"------------------------------------------------------------
zlim <- range(range(values(glcm.sg[["glcm_entropy"]]), na.rm = TRUE), 
                     range(values(glcm.gn[["glcm_entropy"]]), na.rm = TRUE))
par(mfrow=c(1,2))
plot(glcm.sg[["glcm_entropy"]], main = "SoilGrids v2.0", 
     range = zlim, col=(topo.colors(32)))
plot(glcm.gn[["glcm_entropy"]], main = "gNATSGO", 
     range = zlim, col=(topo.colors(32)))
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: histo-equal
#| fig-width: 8
n.class <- 8
# combined values
values.all <- c(values(gn.utm),
                values(sg.utm))
values.all.sort <- sort(values.all)
#
n <- length(values.all) - sum(is.na(values.all))
(cut.positions <- round(n/n.class))
(cuts <- values.all.sort[cut.positions * 1:(n.class-1)])
hist(values.all, breaks=36, main="Histogram equalization",
     xlab = "pH")
abline(v=cuts, col="blue", lwd=2)


## --------------------------------------------------------------------------------------------------------------
#| label: histo-equal-one-by-one
#| fig-width: 14
par(mfrow=c(1,2))
hist(values(gn.utm), breaks=36, main="gNATSGO",
     xlab = "pH")
abline(v=cuts, col="blue", lwd=2)
hist(values(sg.utm), breaks=36, main="SoilGrids v2.0",
     xlab = "pH")
abline(v=cuts, col="blue", lwd=2)
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: classify.setup
(zlim <- c(min(values.all, na.rm = TRUE),
                max(values.all, na.rm=TRUE)))
(cut.names <- cut(zlim, breaks=c(zlim[1], cuts, zlim[2]),
                  ordered_result=TRUE, include.lowest = TRUE)) 
# make sure lowest value is included
#
# common colour ramp
color.ramp <- bpy.colors(n.class+1)
#
(cuts <- round(c(zlim[1], cuts, zlim[2]),2))


## --------------------------------------------------------------------------------------------------------------
#| label: classify-raster-hist
gn.class <- terra::classify(gn.utm, rcl= cuts)
# gn.class <- as.factor(gn.class)
table(values(gn.class))
names(gn.class) <- "class"
sg.class <- terra::classify(sg.utm, rcl= cuts)
table(values(sg.class))
names(sg.class) <- "class"


## --------------------------------------------------------------------------------------------------------------
#| label: show.classified-hist
par(mfrow=c(1, 2))
.l <- range(values(gn.class), na.rm=TRUE)
terra::plot(gn.class,
            col=color.ramp[.l[1]:.l[2]+1], type="classes",
            main="gNATSGO")
.l <- range(values(sg.class), na.rm=TRUE)
terra::plot(sg.class,
            col=color.ramp[.l[1]:.l[2]+1], type="classes",
            main="SG2")
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: ph.classes
range.all <- range(values(gn.utm),
                   values(sg.utm),
                   na.rm = TRUE)
lim.low <- floor(10*range.all[1])/10
lim.low <- ifelse((lim.low %% .2) != 0, lim.low - 0.1, lim.low)
lim.high <- ceiling(10*range.all[2])/10
lim.high <- ifelse((lim.high %% .2) != 0, lim.high + 0.1, lim.high)
(cuts <- seq(lim.low, lim.high, by = 0.2))


## --------------------------------------------------------------------------------------------------------------
#| label: classify-raster-cuts
gn.class <- terra::classify(gn.utm, rcl= cuts)
# gn.class <- as.factor(gn.class)
table(values(gn.class))
names(gn.class) <- "class"
sg.class <- terra::classify(sg.utm, rcl= cuts)
table(values(sg.class))
names(sg.class) <- "class"


## --------------------------------------------------------------------------------------------------------------
#| label: show.classified-cuts
par(mfrow=c(1, 2))
color.ramp <- bpy.colors(length(cuts))
.l <- range(values(gn.class), na.rm=TRUE)
terra::plot(gn.class,
            col=color.ramp[.l[1]:.l[2]+1], type="classes",
            main="gNATSGO pH 0.2 units")
.l <- range(values(sg.class), na.rm=TRUE)
terra::plot(sg.class,
            col=color.ramp[.l[1]:.l[2]+1], type="classes",
            main="SG2 pH 0.2 units")
par(mfrow=c(1,1))


## --------------------------------------------------------------------------------------------------------------
#| label: diffeR
dim(ccm <- diffeR::crosstabm(sg.class, gn.class))
sum(ccm)
prod(dim(gn.class))  # total pixels, includes some NA
ccm[1:5, 1:5]


## --------------------------------------------------------------------------------------------------------------
#| label:  diffeR.p
dim(ccm.p <- diffeR::crosstabm(sg.class, gn.class, percent = TRUE))
sum(ccm.p)


## --------------------------------------------------------------------------------------------------------------
ccm.p[1:5, 1:5]


## --------------------------------------------------------------------------------------------------------------
#| label: diffeR.table
(dt <- diffeR::diffTablej(ccm))


## --------------------------------------------------------------------------------------------------------------
#| label: components-plot
#| fig-width: 4
#| fig-height: 6
#| out-width: "60%"
overallComponentsPlot(ctmatrix = ccm, units = "pixels")


## --------------------------------------------------------------------------------------------------------------
#| label: crosstab-continuous
# whole pH units
print(ccm.c <- crosstabm(sg.utm, gn.utm))
diffTablej(ccm.c)
# to one decimal pH unit
# show the first four reference classes
print((ccm.c10 <- crosstabm(sg.utm*10, gn.utm*10))[1:5, ])
diffTablej(ccm.c10)[1:5, ]


## --------------------------------------------------------------------------------------------------------------
#| label: coma
coma.gn <- lsp_signature(gn.class, type="coma", neighbourhood = 8)
print(coma.gn.matrix <- as.matrix(coma.gn$signature)[[1]])
sum(diag(coma.gn.matrix))/sum(coma.gn.matrix)
coma.sg <- lsp_signature(sg.class, type="coma", neighbourhood = 8)
print(coma.sg.matrix <- as.matrix(coma.sg$signature)[[1]])
sum(diag(coma.sg.matrix))/sum(coma.sg.matrix)


## --------------------------------------------------------------------------------------------------------------
#| label:  metrics.cove
# normalized co-occurence vector 8 x 8
cove.gn <- lsp_signature(gn.class, type="cove", neighbourhood = 8)
cove.sg <- lsp_signature(sg.class, type="cove", neighbourhood = 8)


## --------------------------------------------------------------------------------------------------------------
#| label: silt-map
#| fig-width: 8
(sg.silt <- rast(paste0(file.dir,
                        "SoilGrids250/lat4243_lon-77-76/silt_0-5cm_mean.tif")))
sg.silt <- crop(sg.silt, ext(ext.crop))
values(sg.silt) <- values(sg.silt)/10 # convert from ppt to %
sg.silt.utm <- terra::project(sg.silt, epsg.utm, 
                         res = c(250, 250), method = "bilinear")
sg.silt.utm <- resample(sg.silt.utm, sg.utm) # make extents identical
cuts <- seq(10, 90, by = 5) 
sg.silt.class <- terra::classify(sg.silt.utm, rcl= cuts)
table(values(sg.silt.class))
names(sg.silt.class) <- "class"
plot(sg.silt.class, col = topo.colors(11))


## --------------------------------------------------------------------------------------------------------------
#|.label: coma-cove
coma.sg.silt <- lsp_signature(sg.silt.class, type="coma", neighbourhood = 8)
print(coma.sg.silt.matrix <- as.matrix(coma.sg.silt$signature)[[1]])
sum(diag(coma.sg.silt.matrix))/sum(coma.sg.silt.matrix)
cove.sg.silt <- lsp_signature(sg.silt.class, type="cove", neighbourhood = 8)


## --------------------------------------------------------------------------------------------------------------
#| label: distance-ph-silt
cove.df <- data.frame(cove.sg)$signature[[1]][1,]
cove.df <- rbind(cove.df, cove.sg.silt$signature[[1]][1,])
cove.dists <- round(
  philentropy::distance(cove.df, method = "jensen-shannon", 
                        use.row.names =TRUE, 
                        as.dist.obj = FALSE,
                        diag = FALSE) ,4)
print(cove.dists)


## --------------------------------------------------------------------------------------------------------------
#| label: incove.dist
sg.ph.silt.class <- c(sg.class, sg.silt.class)
incove.sg <- lsp_signature(sg.ph.silt.class,
                              type = "incove",
                              neighbourhood = 8,
                              ordered = TRUE,  # the pH classes are ordered
                              window = 16,
                              normalization = "pdf")  #sum to one
summary(incove.sg.dist <- lsp_to_dist(incove.sg,
                                 dist_fun = "jensen-shannon"))
dim(incove.sg.dist)


## --------------------------------------------------------------------------------------------------------------
#| label: incove.cluster
sg.hclust <- hclust(incove.sg.dist, method = "ward.D2")
plot(sg.hclust, main = "clusters of distance between `incove`")


## --------------------------------------------------------------------------------------------------------------
#| label: incove.cluster.cut
#| fig-width: 8
sg.clusters <- as.factor(cutree(sg.hclust, h = 0.5))  # cutpoint by visual inspection
levels(sg.clusters)
sg.grid.sf = lsp_add_clusters(incove.sg, sg.clusters)
sg.grid.sf$clust <- as.factor(sg.grid.sf$clust)
my.pal <- colorRampPalette(brewer.pal(8, "Accent"))(length(levels(sg.grid.sf$clust)))
ggplot(data = sg.grid.sf) + 
  geom_sf(aes(fill = clust), alpha = 0.7) +
  scale_fill_discrete(type = my.pal) +
  labs(title = "Clusters: distance between integrated co-occurrence vectors",
       fill = "cluster")


## --------------------------------------------------------------------------------------------------------------
#| label: incove-clusters-map-3
p1 <- ggplot(data = sg.grid.sf) + 
  geom_sf(aes(fill = clust), alpha = 0.7) +
  scale_fill_discrete(type = my.pal) +
  labs(fill = "cluster")  +
  theme(legend.position="none")
p2 <- ggplot() +
  tidyterra::geom_spatraster(data = sg.class, aes(fill = class)) +
   theme(legend.position="none")
p3 <- ggplot() +
  tidyterra::geom_spatraster(data = sg.silt.class, aes(fill = class)) +
   theme(legend.position="none")
gridExtra::grid.arrange(p1, p2, p3, nrow=1)


## --------------------------------------------------------------------------------------------------------------
#| label: list_metrics_patch
landscapemetrics::list_lsm(level="patch") %>% print(n=Inf)


## --------------------------------------------------------------------------------------------------------------
#| label: list_metrics_class
landscapemetrics::list_lsm(level="class") %>% print(n=Inf)


## --------------------------------------------------------------------------------------------------------------
#| label: list_metrics_landscape
landscapemetrics::list_lsm(level="landscape") %>% print(n=Inf)


## --------------------------------------------------------------------------------------------------------------
#| label: compute.metrics.check
check_landscape(gn.class)
check_landscape(sg.class)


## --------------------------------------------------------------------------------------------------------------
#| label:  show.patches.global
#| fig-width: 8
show_patches(gn.class, class = "global")
show_patches(sg.class, class = "global")


## --------------------------------------------------------------------------------------------------------------
#| label:  show.patches.all
#| fig-width: 14
show_patches(sg.class, class = "all", nrow = 3)
show_patches(gn.class, class = "all", nrow = 3)


## --------------------------------------------------------------------------------------------------------------
lst <- paste0("lsm_l_", c("shdi", "shei", "lsi", "ai",  "frac_mn"))
ls.metrics.gn <- calculate_lsm(gn.class, what=lst)
ls.metrics.sg <- calculate_lsm(sg.class, what=lst)
metrics.table <- data.frame(product=c("gNATSGO", "SG2"),
                            rbind(round(ls.metrics.gn$value, 3),
                                  round(ls.metrics.sg$value, 3)))
names(metrics.table)[2:6] <- ls.metrics.gn$metric
metrics.table


## --------------------------------------------------------------------------------------------------------------
#| label: difference.cove
names(cove.gn)
cove.df <- data.frame(cove.gn)$signature[[1]][1,]
cove.df <- rbind(cove.df, cove.sg$signature[[1]][1,])
cove.dists <- round(
  philentropy::distance(cove.df, method = "jensen-shannon", 
                        use.row.names =TRUE, 
                        as.dist.obj = FALSE,
                        diag = FALSE) ,4)
print(cove.dists)


## --------------------------------------------------------------------------------------------------------------
#| label: lsp.compare
lsp_compare(gn.class, sg.class, 
            type = "cove", dist_fun = "jensen-shannon",
            neighbourhood = 8, # queen's case
            output = "sf")


## --------------------------------------------------------------------------------------------------------------
#| label: lsp.compare.window
#| fig-width: 6
dim(sg.class)
x.dim <- diff(range(st_bbox(sg.class)[c(1,3)]))
y.dim <- diff(range(st_bbox(sg.class)[c(2,4)]))
(compare.16 <- lsp_compare(gn.class, sg.class, 
                           type = "cove", dist_fun = "jensen-shannon",
                           neighbourhood = 8, # queen's case
                           window = 16,
                           output = "sf"))
ggplot(data = compare.16) +
  geom_sf(aes(fill = dist)) +
  labs(title = "Distance between co-occurrence vectors, pH class")


## --------------------------------------------------------------------------------------------------------------
#| label: cove-patterns-plot
p1 <- ggplot(data = compare.16) + 
    geom_sf(aes(fill = dist))  +
   theme(legend.position="none")
p2 <- ggplot() +
  tidyterra::geom_spatraster(data = sg.class, aes(fill = class)) +
   theme(legend.position="none")
p3 <- ggplot() +
  tidyterra::geom_spatraster(data = gn.class, aes(fill = class)) +
   theme(legend.position="none")
gridExtra::grid.arrange(p1, p2, p3, nrow=1)


## --------------------------------------------------------------------------------------------------------------
#| label: polygonize
gn.poly <- terra::as.polygons(gn.class,
                              aggregate= TRUE,
                              values = TRUE,
                              dissolve = TRUE)
sg.poly <- terra::as.polygons(sg.class,
                              aggregate= TRUE,
                              values = TRUE, 
                              dissolve=TRUE)


## --------------------------------------------------------------------------------------------------------------
#| label: convert.polygons.sf
gn.sf <- st_as_sf(gn.poly)
gn.sf <- st_cast(gn.sf, "MULTIPOLYGON")
#
sg.sf <- st_as_sf(sg.poly)
sg.sf <- st_cast(sg.sf, "MULTIPOLYGON")


## --------------------------------------------------------------------------------------------------------------
#| label: make.valid
st_is_valid(gn.sf, reason=TRUE)
gn.sf.v <- sf::st_make_valid(gn.sf)  |> st_cast("MULTIPOLYGON")
#
st_is_valid(sg.sf, reason=TRUE)
sg.sf.v <- sf::st_make_valid(sg.sf)  |> st_cast("MULTIPOLYGON")


## --------------------------------------------------------------------------------------------------------------
#| label: make_class_maps
classes.both <- union(values(sg.poly)$class, values(gn.poly)$class)
my.pal <- colorRampPalette(brewer.pal(8, "PuBu"))(length(classes.both))
g0 <- ggplot(data=gn.sf.v) +
  geom_sf(aes(fill = class)) +
  coord_sf(crs = st_crs(gn.sf)) +
  labs(title = "gNATSGO")  +
  scale_fill_manual(values = my.pal, drop=TRUE,
                    limits = levels(gn.sf.v$class)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
g1 <- ggplot(data=sg.sf.v) +
  geom_sf(aes(fill = class)) +
  coord_sf(crs = st_crs(sg.sf)) +
  labs(title = "SG2")  +
  scale_fill_manual(values = my.pal, drop=TRUE,
                    limits = levels(gn.sf.v$class)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
grid.arrange(g0,g1, nrow=1, ncol=2)


## ----gNATSGO8.sg-----------------------------------------------------------------------------------------------
#| label: gnsg.regions
regions.sg.gn <- sabre::vmeasure_calc(x = gn.sf.v, 
                                 y = sg.sf.v, 
                                 x_name = class, y_name = class)
print(regions.sg.gn)
names(regions.sg.gn)
names(regions.sg.gn$map1)
attr(regions.sg.gn, "precision")  # NULL, means a system default


## --------------------------------------------------------------------------------------------------------------
#| label: vmaps.gnatsgo.sg.homogeneity
terra::plot(regions.sg.gn$map1["rih"], main = "Inhomogeneity -- SG2 vs. gNATSGO")


## --------------------------------------------------------------------------------------------------------------
#| label: vmaps.gnatsgo.sg.completeness
terra::plot(regions.sg.gn$map2["rih"], main = "Incompleteness -- SG2 vs. gNATSGO")


## --------------------------------------------------------------------------------------------------------------
#| label: SLIC-source
#| fig-width: 6
ggplot() +
  geom_spatraster(data=sg.utm) +
  labs(fill = "pH")


## --------------------------------------------------------------------------------------------------------------
#| label: supercells-not-compact
#| fig-width: 6
sg.utm.50 = supercells(sg.utm, k = 50, compactness = 0.5)
ggplot(data=sg.utm.50) +
  geom_sf(aes(fill = phh2o_0.5cm_mean)) +
  labs(fill = "mean pH")


## --------------------------------------------------------------------------------------------------------------
#| label: supercells-compact
#| fig-width: 6
sg.utm.50 = supercells(sg.utm, k = 50, compactness = 5)
ggplot(data=sg.utm.50) +
  geom_sf(aes(fill = phh2o_0.5cm_mean)) +
  labs(fill = "mean pH")


## --------------------------------------------------------------------------------------------------------------
#| label: supercells-multiple
#| fig-width: 6
r <- c(sg.utm, sg.silt.utm)
r.50 = supercells(r, k = 50, compactness = 0.5)
ggplot(data=r.50) +
  geom_sf(aes(fill = phh2o_0.5cm_mean)) +
  labs(fill = "mean pH") +
  scale_fill_continuous(type = "viridis")
ggplot(data=r.50) +
  geom_sf(aes(fill = silt_0.5cm_mean)) +
  labs(fill = "silt ppt") +
  scale_fill_continuous(type = "viridis")


## --------------------------------------------------------------------------------------------------------------
mu.template <- rast(mu.poly, res=c(20,20))
dim(mu.template)
mu.raster <- rasterize(mu.poly, mu.template, field="mukey")
summary(mu.raster)
check_landscape(mu.raster)


## --------------------------------------------------------------------------------------------------------------
#| label:  show.landscape.gn
#| fig-width: 14
#| fig-height: 10
head(unique(mu.raster$mukey))
show_patches(mu.raster, class = "global")
 show_patches(mu.raster, class = 295575)
(ix <- grep("Mardin", mu.key$muname, fixed = TRUE))
print(mu.key[ix, ])
show_patches(mu.raster, 
             class = c(mu.key[ix, "mukey"]), nrow = 3)


## --------------------------------------------------------------------------------------------------------------
gn.20 <- mu.raster


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.30.landscape
lst <- paste0("lsm_l_", c("shdi", "shei", "lsi", "ai",  "frac_mn"))
print(ls.metrics.gn20 <- calculate_lsm(gn.20, what=lst)[, c("metric", "value")])


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.20.class
c_pland.20 <- calculate_lsm(gn.20, what="lsm_c_pland")
head(sort(c_pland.20$value, decreasing = TRUE), 24)
ix <- order(c_pland.20$value, decreasing = TRUE)
head(c_pland.20$class[ix], 24)


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.20.patch.area
# each patch
head(sort(area.20 <-  
            calculate_lsm(gn.20, what="lsm_p_area")$value, decreasing = TRUE))
quantile(area.20 , seq(0,1,by=.12))


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.20.patch
# each patch
head(ncore.20 <- calculate_lsm(gn.20, what="lsm_p_ncore"))
# summarize by map unit
ncore.20.summary <- ncore.20 %>% group_by(class) %>%
  summarize(max_cores = max(value)) %>%
  arrange(class)
print(sort(ncore.20.summary$max_cores, decreasing = TRUE))
ix <- order(ncore.20.summary$max_cores, decreasing = TRUE)
cbind(mu.key[ix[1:8], ], ncore = ncore.20.summary$max_cores[ix[1:8]])


## --------------------------------------------------------------------------------------------------------------
#| label: show-patches-Erie-20
show_patches(gn.20, class = mu.key[ix[1], "mukey"])


## --------------------------------------------------------------------------------------------------------------
#| label: show-ncore-20
show_lsm(gn.20, "lsm_p_ncore")


## --------------------------------------------------------------------------------------------------------------
#| label: scale.change.100
gn.100 <- terra::aggregate(gn.20, fact=5, fun="modal")


## --------------------------------------------------------------------------------------------------------------
print(ls.metrics.gn20[, c("metric", "value")])
print(ls.metrics.gn100 <- calculate_lsm(gn.100, what=lst)[, c("metric", "value")])


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.100.class
c_pland.100 <- calculate_lsm(gn.100, what="lsm_c_pland")
head(sort(c_pland.20$value, decreasing = TRUE), 24)
head(sort(c_pland.100$value, decreasing = TRUE), 24)


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.100.patch.area
# each patch
head(sort(area.100 <-  
            calculate_lsm(gn.100, what="lsm_p_area")$value, decreasing = TRUE))
quantile(area.20, seq(0,1,by=.12))
quantile(area.100, seq(0,1,by=.12))


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.100.patch
# each patch
head(ncore.100 <- calculate_lsm(gn.100, what="lsm_p_ncore"))
# summarize by map unit
ncore.100.summary <- ncore.100 %>% group_by(class) %>%
  summarize(max_cores = max(value)) %>%
  arrange(class)
print(sort(ncore.20.summary$max_cores, decreasing = TRUE))
print(sort(ncore.100.summary$max_cores, decreasing = TRUE))
ix <- order(ncore.100.summary$max_cores, decreasing = TRUE)
cbind(mu.key[ix[1:8], ], ncore = ncore.100.summary$max_cores[ix[1:8]])


## --------------------------------------------------------------------------------------------------------------
#| label: show-patches-Erie-100
show_patches(gn.100, class = mu.key[ix[1], "mukey"])


## --------------------------------------------------------------------------------------------------------------
#| label: show-ncore-100
show_lsm(gn.100, "lsm_p_ncore")


## --------------------------------------------------------------------------------------------------------------
#| label: scale.change.300
gn.300 <- terra::aggregate(gn.100, fact=3, fun="modal")


## --------------------------------------------------------------------------------------------------------------
print(ls.metrics.gn20[, c("metric", "value")])
print(ls.metrics.gn100[, c("metric", "value")])
print(ls.metrics.gn300 <- calculate_lsm(gn.300, what=lst)[, c("metric", "value")])


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.300.class
c_pland.300 <- calculate_lsm(gn.300, what="lsm_c_pland")
head(sort(c_pland.20$value, decreasing = TRUE), 24)
head(sort(c_pland.100$value, decreasing = TRUE), 24)
head(sort(c_pland.300$value, decreasing = TRUE), 24)


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.300.patch.area
# each patch
head(sort(area.300 <-  
            calculate_lsm(gn.300, what="lsm_p_area")$value, decreasing = TRUE))
quantile(area.20, seq(0,1,by=.12))
quantile(area.100, seq(0,1,by=.12))
quantile(area.300, seq(0,1,by=.12))


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.300.patch
# each patch
head(ncore.300 <- calculate_lsm(gn.300, what="lsm_p_ncore"))
# summarize by map unit
ncore.300.summary <- ncore.300 %>% group_by(class) %>%
  summarize(max_cores = max(value)) %>%
  arrange(class)
print(ncore.20.summary)
print(ncore.100.summary)
print(ncore.300.summary)


## --------------------------------------------------------------------------------------------------------------
#| label: show-patches-Erie-300
show_patches(gn.300, class = mu.key[ix[1], "mukey"])


## --------------------------------------------------------------------------------------------------------------
#| label: show-ncore-300
show_lsm(gn.300, "lsm_p_ncore")


## ----echo=FALSE------------------------------------------------------------------------------------------------
ix <- grep("Mardin", mu.key$muname, fixed = TRUE)
print(mu.key[ix, ])


## --------------------------------------------------------------------------------------------------------------
#| label: dominant-series
length(names <- mu.key$muname)
# first name by spaces
names <- strsplit(names, " ")
head(names)
names.unique <- unique(unlist(lapply(names, function(x) x[1])))
print(names.unique)


## --------------------------------------------------------------------------------------------------------------
#| label: group
mu.poly.general <- mu.poly
names(values(mu.poly.general))
(l.all <- length(unique(values(mu.poly)$mukey)))
names.all <- unlist(lapply(names, function(x) x[1]))
(l.general <- length(names.unique <- unique(names.all)))
for (name in names.unique) {
  ix <- which(name == names.all)
  # map units in this group
  keys <- mu.key$mukey[ix]
  # polygons of these map units
  ix.polys <- which(values(mu.poly.general)$mukey %in% keys)
  # rename the key
  values(mu.poly.general)$mukey[ix.polys] <- keys[1]
}
names(values(mu.poly.general))


## --------------------------------------------------------------------------------------------------------------
#| label: dissolve-polygons
mu.poly.general <- terra::aggregate(mu.poly.general,
                                    by = "mukey",
                                    fun = "sum")
plot(mu.poly.general, y=2,
     type = "continuous",
     main = "SSURGO map units, area in acres")


## --------------------------------------------------------------------------------------------------------------
#| label: rasterize-general
mu.raster.general <- rasterize(mu.poly.general, mu.template, field="mukey")
check_landscape(mu.raster.general)
gn.20.g <- mu.raster.general


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.20.landscape.genereal
ls.metrics <- ls.metrics.gn.20.g <- 
  calculate_lsm(gn.20.g, what=lst)[, c("metric", "value")]
ls.metrics$detailed <- ls.metrics.gn20$value
names(ls.metrics) <- c("metric", "generalized", "detailed")
print(ls.metrics)


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.20.class.general
c_pland.20.g <- calculate_lsm(gn.20.g, what="lsm_c_pland")
head(sort(c_pland.20.g$value, decreasing = TRUE), 24)
head(sort(c_pland.20$value, decreasing = TRUE), 24)
ix <- order(c_pland.20.g$value, decreasing = TRUE)
head(c_pland.20.g$class[ix], 24)
head(c_pland.20$class[ix], 24)


## --------------------------------------------------------------------------------------------------------------
#| label: metrics.20.patch.area.general
# each patch
head(sort(area.20.g <-  
            calculate_lsm(gn.20.g, 
                          what="lsm_p_area")$value, 
          decreasing = TRUE))
quantile(area.20.g , seq(0,1,by=.12))
quantile(area.20 , seq(0,1,by=.12))


## --------------------------------------------------------------------------------------------------------------
#| label: cores-general
ncore.g <- calculate_lsm(gn.20.g, 
                         what="lsm_p_ncore")
ncore.g.summary <- ncore.g %>% group_by(class) %>%
  summarize(max_cores = max(value)) %>%
  arrange(class)
print(ncore.g.summary)
g1 <- show_lsm(gn.20.g, "lsm_p_ncore")
g2 <- show_lsm(gn.20, "lsm_p_ncore")
par(mfrow = c(1,2))
g1; g2
par(mfrow = c(1,2))


## ----eval=FALSE------------------------------------------------------------------------------------------------
## require(rgeopat2)

