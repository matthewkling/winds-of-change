
library(raster)
library(tidyverse)

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

clim <- stack(h, m) %>% "/"(10) %>% writeRaster("data/chelsa/tmean_1km.tif")
