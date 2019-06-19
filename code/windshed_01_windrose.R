

library(windshed)
library(tidyverse)
library(raster)


cfsr_rose <- function(infile, outdir, ncores, weighting){
      message(infile)
      
      # file admin
      outfile <- paste0(outdir, "/", basename(infile))
      outfile <- sub("grb2", "tif", outfile)
      if(any(file.exists(c(outfile,
                           sub("\\.tif", paste0("_", ncores, ".tif"), outfile))))){
            return("skipping")}
      
      # load inventory file and identify layers to process
      # (the "analysis" and "6-hour forecast" timesteps overlap)
      inv <- paste0(infile, ".inv") %>%
            readLines()
      process <- which(!grepl("6 hour fcst", inv))
      
      # generate windrose
      x <- infile %>%
            brick() %>%
            subset(process) %>%
            windrose_rasters(ncores=ncores, 
                             weighting=weighting,
                             outfile=outfile)
      return("success")
}


# hourly input data
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]


# velocity
intdir <- "data/roses_velocity/cfsr_monthly"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, weighting="velocity")
f <- list.files(intdir, full.names=T) %>%
      lapply(stack) %>%
      Reduce("+", .) %>%
      writeRaster("data/roses_velocity/cfsr_climatology/roses_cfsr_1980s.tif", overwrite=T)


# force
intdir <- "data/roses_force/cfsr_monthly"
f <- f[grepl("gdas\\.20", f)]
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, weighting="force")

f <- list.files(intdir, full.names=T)[1:12] %>%
      lapply(stack) %>%
      Reduce("+", .) %>%
      "/"(24 * 365 * 1) %>%
      writeRaster("data/roses_force/cfsr_climatology/roses_cfsr_1980s.tif", overwrite=T)

