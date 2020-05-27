
# prepare temperature, elevation, and coastal distance data

library(raster)
library(tidyverse)


##### prep 1km CHELSA data #####

for(var in c("tmean", "bio5", "bio6", "ppt", "bio13", "bio14")){
      
      # the bits that differ by variable
      v <- var %>% switch(tmean = "_tas_", bio5 = "_tasmax_", bio6 = "_tasmin_",
                          ppt = "_pr_", bio13 = "_pr_", bio14 = "_pr_")
      v1 <- var %>% switch(tmean = "1", bio5 = "5", bio6 = "6",
                           ppt = "12", bio13 = "13", bio14 = "14")
      v2 <- var %>% switch(tmean = "tmean", bio5 = "bio5", bio6 = "bio6",
                           ppt = "precip", bio13 = "bio13", bio14 = "bio14")
      v3 <- var %>% switch(tmean = "temperature", bio5 = "bio5", bio6 = "bio6",
                           ppt = "precipitation", bio13 = "bio13", bio14 = "bio14")
      fx <- var %>% switch(tmean = mean, bio5 = max, bio6 = min,
                           ppt = mean, bio13 = max, bio14 = min)
      
      # future -- ensemble means across 10 gcms
      f <- data.frame(path = list.files("F:/chelsa/cmip5", full.names=T),
                      stringsAsFactors = F) %>%
            filter(grepl("rcp85", path),
                   grepl("2080", path),
                   grepl(v, path)) %>%
            mutate(model = sub(paste0("F:/chelsa/cmip5/CHELSA", v, "mon_"), "", path),
                   model = substr(model, 1, 10)) %>%
            split(.$model)
      
      # average months into annual for each model
      mm <- lapply(f, function(x){
            message(x$path[1])
            fx(stack(x$path))
      })
      
      # average across models
      m <- stack(mm) %>% mean()
      
      # historic
      h <- paste0("F:/chelsa/bio19/CHELSA_bio10_", v1, ".tif") %>%
            raster()
      m <- mask(m, h)
      
      out <- paste0("data/geographic/raw/chelsa_", v2, "_1km.tif")
      clim <- stack(h, m) %>% "/"(10) %>% writeRaster(out, overwrite=T)
      
      
      ##### upscale to CFSR grid #####
      
      r <- brick(out)
      template <- raster("f:/cfsr/land.tif") %>% rotate()
      
      # aggregate a bit (for computation reasons) and convert to points
      p <- r %>% aggregate(4) %>% rasterToPoints() %>% as.data.frame()
      coordinates(p) <- c("x", "y")
      crs(p) <- crs(r)
      u <- rasterize(p, template, 
                     field=c(paste0("chelsa_", v2, "_1km.1"), paste0("chelsa_", v2, "_1km.2")), 
                     fun=mean, na.rm=T)
      
      if(var %in% c("ppt", "bio13", "bio14")){
            u[[2]] <- u[[2]] * 10
            u <- log10(u + 1)
      }
      
      out2 <- paste0("data/geographic/processed/", v3, ".tif")
      writeRaster(u, out2, overwrite = T)
}



######## elevation ########
r <- raster("F:/chelsa/elevation/mn30_grd/w001001.adf") 
p <- r %>% aggregate(4) %>% rasterToPoints() %>% as.data.frame()
coordinates(p) <- c("x", "y")
crs(p) <- crs(r)
u <- rasterize(p, template, field = "w001001", fun = mean, na.rm=T)
writeRaster(u, "data/geographic/processed/elevation.tif")


####### coastal distance ########
# data from https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/
r <- raster("data/geographic/raw/GMT_intermediate_coast_distance_01d/GMT_intermediate_coast_distance_01d.tif") 
p <- r %>% reclassify(c(0, Inf, NA)) %>%
      aggregate(10) %>% rasterToPoints() %>% as.data.frame()
coordinates(p) <- c("x", "y")
crs(p) <- crs(r)
u <- rasterize(p, template, field = "layer", fun = mean, na.rm=T)
u <- u * -1
writeRaster(u, "data/geographic/processed/coastal_distance.tif", overwrite=T)
