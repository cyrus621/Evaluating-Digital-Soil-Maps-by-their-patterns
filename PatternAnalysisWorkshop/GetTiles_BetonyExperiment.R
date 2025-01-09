## test from Betony
library(XML)
library(gdalUtilities)

voi = "bdod" # variable of interest
depth = "0-5cm"
quantile = "mean"

voi_layer = paste(voi,depth,quantile, sep="_") # layer of interest 


wcs_path = paste0("https://maps.isric.org/mapserv?map=/map/",voi,".map") # Path to the WCS. See maps.isric.org
wcs_service = "SERVICE=WCS"
wcs_version = "VERSION=2.0.1" # This works for gdal >=2.3; "VERSION=1.1.1" works with gdal < 2.3.

bb=c(8397053, 2449029, 8535870, 2337709) # Example bounding box (homolosine)
igh="+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs" # proj string for Homolosine projection

wcs = paste(wcs_path,wcs_service,wcs_version,sep="&") # This works for gdal >= 2.3

l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", voi_layer, parent=l1)

# Save to local disk
xml.out = "./sg_david_bdod.xml"
saveXML(l1, file = xml.out)

# Download raster as GeoTIFF (Warning: it can be large!)
file.out <- './test_bdod_david.tif'

gdal_translate(xml.out, file.out,
               tr=c(250,250), projwin=bb,
               projwin_srs =igh, co=c("TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","BIGTIFF=YES")
)

