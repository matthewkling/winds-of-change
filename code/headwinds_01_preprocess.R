

# process CFSR data:

# 1. calculate monnhly means of U and V wind components
# 2. calculate long-term means
# 2. calculate wind direction


library(raster)
#library(ncdf4)
library(tidyverse)
library(rNOMADS)




# binary land-ocean layer
land <- raster("f:/cfsr/soilm1.gdas.197901.grb2") %>% 
      reclassify(c(-Inf, Inf, 1)) %>%
      writeRaster("f:/cfsr/land.tif")



modir <- "E:/wind/winds_of_change/tailwinds/data/cfsr_monthly"


#### calculate monthly means from raw hourly data

# wind
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
means <- function(file){
      #file <- f[1]
      message(basename(file))
      out <- paste0(modir, "/wnd10m/", sub("grb2", "tif", basename(file)))
      if(file.exists(sub("\\.tif", "_v.tif", out))) return("pass")
      w <- stack(file)
      even <- function(x) x %% 2 == 0
      u <- w[[which(!even(1:nlayers(w)))]] %>% calc(mean, filename=sub("\\.tif", "_u.tif", out))
      v <- w[[which(even(1:nlayers(w)))]] %>% calc(mean, filename=sub("\\.tif", "_v.tif", out))
}
if(T) lapply(f[241:360], means)

# temperature
f <- list.files("f:/CFSR/tmp2m", full.names=T)
f <- f[!grepl("inv", f)]
means <- function(file){
      #file <- f[1]
      message(basename(file))
      out <- paste0(modir, "/tmp2m/", sub("grb2", "tif", basename(file)))
      if(file.exists(out)) return("pass")
      w <- stack(file) %>% calc(mean, filename=out)
}
if(F) lapply(f, means)


