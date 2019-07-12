

library(windscape)
library(tidyverse)
library(raster)


cfsr_rose <- function(infile, outdir, ncores, p){
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
                             p=p,
                             outfile=outfile)
      return("success")
}


# hourly input data
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
f <- f[grepl("gdas\\.20", f)]



# drag
intdir <- "data/windrose/monthly_p2"
map(f[1:120], possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=2)
ny <- 10
s <- list.files(intdir, full.names=T)[1:(72*ny)] %>% 
      lapply(stack) %>%
      Reduce("+", .) %>%
      "/"(24 * 365 * ny) %>%
      writeRaster("data/windrose_p2_2000s.tif", overwrite=T)



stop("wootwoot")

# velocity
intdir <- "data/roses_velocity/cfsr_monthly"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=1)
ny <- 10
s <- list.files(intdir, full.names=T)[1:(72*ny)] %>% 
      lapply(stack) %>%
      Reduce("+", .) %>%
      "/"(24 * 365 * ny) %>%
      writeRaster("data/roses_velocity/cfsr_climatology/roses_cfsr_2000s.tif", overwrite=T)



# force
intdir <- "data/roses_force/cfsr_monthly"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=3)
ny <- 10
s <- list.files(intdir, full.names=T)[1:(72*ny)] %>% 
      lapply(stack) %>%
      Reduce("+", .) %>%
      "/"(24 * 365 * ny) %>%
      writeRaster("data/roses_force/cfsr_climatology/roses_cfsr_2000s.tif", overwrite=T)

