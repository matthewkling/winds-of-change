
library(raster)
library(tidyverse)



# prepare temperature, elevation, and coastal distance data


##### prep 1km CHELSA data #####

# future -- mean across 10 gcms
f <- data.frame(path = list.files("F:/chelsa/cmip5", full.names=T),
                stringsAsFactors = F) %>%
      filter(grepl("rcp85", path),
             grepl("2080", path),
             grepl("_tas_", path)) %>%
      mutate(model = sub("F:/chelsa/cmip5/CHELSA_tas_mon_", "", path),
             model = substr(model, 1, 10)) %>%
      split(.$model)

mm <- lapply(f, function(x){
      message(x$path[1])
      mean(stack(x$path))
})

m <- stack(mm) %>% mean()


# historic
h <- raster("F:/chelsa/bio19/CHELSA_bio10_1.tif")
m <- mask(m, h)

clim <- stack(h, m) %>% "/"(10) %>% writeRaster("data/geographic/raw/chelsa_tmean_1km.tif", overwrite=T)



##### upscale to CFSR grid #####

r <- brick("data/geographic/raw/chelsa_tmean_1km.tif")
template <- raster("f:/cfsr/land.tif") %>% rotate()

# aggregate a bit (for computation reasons) and convert to points
p <- r %>% aggregate(4) %>% rasterToPoints() %>% as.data.frame()
coordinates(p) <- c("x", "y")
crs(p) <- crs(r)
u <- rasterize(p, template, field=c("tmean_1km.1", "tmean_1km.2"), fun=mean, na.rm=T)
writeRaster(u, "data/geographic/processed/temperature.tif")





######## elevation ########
r <- raster("F:/chelsa/elevation/mn30_grd/w001001.adf") 
p <- r %>% aggregate(4) %>% rasterToPoints() %>% as.data.frame()
coordinates(p) <- c("x", "y")
crs(p) <- crs(r)
u <- rasterize(p, template, field="w001001", fun=mean, na.rm=T)
writeRaster(u, "data/geographic/processed/elevation.tif")


####### coastal distance ########
r <- raster("data/geographic/raw/GMT_intermediate_coast_distance_01d/GMT_intermediate_coast_distance_01d.tif") 
p <- r %>% reclassify(c(0, Inf, NA)) %>%
      aggregate(10) %>% rasterToPoints() %>% as.data.frame()
coordinates(p) <- c("x", "y")
crs(p) <- crs(r)
u <- rasterize(p, template, field="layer", fun=mean, na.rm=T)
u <- u * -1
writeRaster(u, "data/geographic/processed/coastal_distance.tif", overwrite=T)
