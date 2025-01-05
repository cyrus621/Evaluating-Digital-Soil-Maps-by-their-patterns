# Get a set of SG250 tiles (properties, depth slices)
# run the sections of interest
# can follow this with SoilGrids250_MakeRasterStack.Rmd
require(rmarkdown)
# set lower-right corner of the area of interest, the tile size,
lrc.long <- 78; lrc.lat <-  10; tile.size <- 1
# set the quantile
# quantile.list <- c("Q0.05", "Q0.5", "Q0.95", "mean")
quantile <- 4  # mean



# Example: one property, all depth slices
# voi.list.sg <- c("clay", "silt", "sand", "phh2o", "cec", "soc", "bdod", "cfvo", "nitrogen", "ocd")
# depth.list <- paste0(c("0-5", "5-15", "15-30", "30-60", "60-100", "100-200"),"cm")
voi <- 5 # example: cec
tmp <- lapply(1:6, function (i) {
  rmarkdown::render("SoilGrids250_WCS_import.Rmd", output_format = NULL, 
                    output_file = tempfile(), # don't write any output
                    params = list(lrc_long = lrc.long, lrc_lat = lrc.lat, 
                                  size = tile.size, 
                                  quantile.n = quantile, 
                                  voi.n = voi,
                                  depth.n = i))
})

# Example: all properties, one depth slice
# voi.list.sg <- c("clay", "silt", "sand", "phh2o", "cec", "soc", "bdod", "cfvo", "nitrogen", "ocd")
# depth.list <- paste0(c("0-5", "5-15", "15-30", "30-60", "60-100", "100-200"),"cm")
depth <- 1 # example: surface layer
tmp <- lapply(1:10, function (i) {
  rmarkdown::render("SoilGrids250_WCS_import.Rmd", output_format = NULL, 
                    output_file = tempfile(), # don't write any output
                    params = list(lrc_long = lrc.long, lrc_lat = lrc.lat, 
                                  size = tile.size, 
                                  quantile.n = quantile, 
                                  voi.n = i,
                                  depth.n = depth))
})


# Example: selected properties, all depth slices
# voi.list.sg <- c("clay", "silt", "sand", "phh2o", "cec", "soc", "bdod", "cfvo", "nitrogen", "ocd")
# depth.list <- paste0(c("0-5", "5-15", "15-30", "30-60", "60-100", "100-200"),"cm")
tmp <- lapply(c(2,4,5,8), function (i) {
  lapply((1:6), function (j) {
    rmarkdown::render("SoilGrids250_WCS_import.Rmd", output_format = NULL, 
                      output_file = tempfile(), # don't write any output
                      params = list(lrc_long = lrc.long, lrc_lat = lrc.lat, 
                                    size = tile.size, 
                                    quantile.n = quantile, 
                                    voi.n = i,
                                    depth.n = j))
  })
})


