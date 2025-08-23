### Get DEM ###

library(elevatr)
library(sf)
library(terra)

# Define AOI
bbox <- st_bbox(c(xmin = -114.9, ymin = 49.6, xmax = -114.6, ymax = 49.95), crs = st_crs(4326))
aoi <- st_as_sfc(bbox)
aoi <- st_sf(geometry = aoi)

# Plot AOI
plot(aoi, main = "Area of Interest")

# Get DEM (zoom level 11)
dem <- get_elev_raster(aoi, z = 11, clip = "bbox")

# Plot DEM
plot(dem, main = "Digital Elevation Model")

terra::writeRaster(dem, "dat/demz11.tiff", overwrite = T)
