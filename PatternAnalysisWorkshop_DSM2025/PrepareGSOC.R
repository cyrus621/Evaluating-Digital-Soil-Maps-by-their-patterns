# Prepare inputs from GSOC map for GeoPAT
# D G Rossiter 16-Aug-2024

# Crop, and export layers from GLOSIS, as inputs to GeoPAT
# Global SOC stocks 0-30 cm t ha^{-1}

# packages
require(terra)

## Manually download from FAO https://data.apps.fao.org/glosis/?lang=en
## and unzip into the `dsm_path`.

##

# input from FAO
dsm_path <- "/Users/rossiter/ds_reference/GSOCmap"
# output -- for GeoPAT input
ds_path <- "/Users/rossiter/tmp/GeoPAT/GloSIS/ds"

# source map
input_TIFF <- file.path(dsm_path, "GSOCmap1.6.1.tif")
(r <- rast(input_TIFF))


# select a study area
# 1200 x 1200 cells; approx. 1 km resolution = 0.008333333 deg.
(delta.ll <- res(r)*1200) # 10 x 10 degrees

# study area extent
# 1 -- Sonora/New Mexico/Texas

# target crop
ext.crop <- ext(c(-110, -100, 27, 37))
# larger area so we can get a square
ext.crop.extended <- ext(c(-110.5, -99.5, 26.5, 37.5))
(r.nm.ext <- crop(r, ext.crop.extended))
plot(r.nm.ext, range = c(0,100))

# project
nm_albers = crs("+proj=aea +lat_1=29 +lat_2=35 +lat_0=32 +lon_0=-105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
(r.nm.ext.p <- project(r.nm.ext, nm_albers, method = "bilinear", res = 1000))
plot(r.nm.ext.p, range = c(0,100))

# make square -- 800 x 800 km
ext.crop.p <- ext(c(-450000, +450000, -450000, +450000))
# corners in Long/Lat
v <- vect(ext.crop.p, crs=nm_albers)
project(v, "EPSG:4326")

(r.nm.p <- crop(r.nm.ext.p, ext(ext.crop.p)))
plot(r.nm.p, range = c(0,100))

hist(values(r.nm.p), main = "SOC stock, t/ha")
quantile(values(r.nm.p), c(0, .25, .5, .8, .9, .95, .99, .995, .998, .999, .9995, 1), na.rm = TRUE)
# limit upper outliers to 0.9995 quantile.
(r.nm.p <- clamp(r.nm.p, lower = -Inf, upper = quantile(values(r.nm.p), 0.9995, na.rm = TRUE)))
plot(r.nm.p)

# classify for GeoPAT -- SOC 
(r.nm.cl <- classify(r.nm.p, seq(0, ceiling(max(values(r.nm.p), na.rm = TRUE)), by = 4)))
plot(r.nm.cl)
table(values(r.nm.cl))


# export
# output cropped and projected source map (continuous)
output_TIFF <- file.path(ds_path, "gsoc_nm_p.tif")
writeRaster(r.nm.p, output_TIFF, overwrite = TRUE)
# output for input to GeoPAT (classes)
output_TIFF <- file.path(ds_path, "gsoc_nm.tif")
writeRaster(r.nm.cl, output_TIFF, overwrite = TRUE)

### --------------------------------------------------------

# 2 -- central Europe
ext.crop <- ext(c(11, 24, 46, 53))
r.ce <- crop(r, ext.crop)
r.ce
plot(r.ce, range = c(0,200))

# custom Albers equal-area for this region
ce_albers = crs("+proj=aea +lat_1=48 +lat_2=51 +lat_0=46 +lon_0=17.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
r.ce.p <- project(r.ce, ce_albers, method = "bilinear", res = 1000)
r.ce.p
plot(r.ce.p, range = c(0,200))

hist(r.ce.p, breaks = 24)
# square root transformation
hist(sqrt(r.ce.p))
r.ce.p <- sqrt(r.ce.p)
r.ce.p
plot(r.ce.p)

# classify for GeoPAT -- SOC 
r.ce.cl <- classify(r.ce.p, seq(0, ceiling(max(values(r.ce.p), na.rm = TRUE)), by = 1))
table(values(r.ce.cl))

# export
# output cropped and projected source map (continuous)
output_TIFF <- file.path(ds_path, "gsoc_ce_sqrt_p.tif")
writeRaster(r.ce.p, output_TIFF, overwrite = TRUE)
# output for input to GeoPAT (classes)
output_TIFF <- file.path(ds_path, "gsoc_ce_sqrt.tif")
writeRaster(r.ce.cl, output_TIFF, overwrite = TRUE)

