

library(windshed)
library(tidyverse)
library(raster)


# generate a windrose dataset for a single month of cfsr data
cfsr_rose <- function(infile, outdir){
      #infile <- f[2]
      message(paste(basename(infile), "--", Sys.time()))
      outfile <- paste0(outdir, "/", sub("grb2", "tif", basename(infile)))
      if(file.exists(outfile)) return("skipped")
      
      # load and collate data
      w <- stack(infile)
      even <- function(x) x %% 2 == 0
      u <- w[[which(!even(1:nlayers(w)))]]
      v <- w[[which(even(1:nlayers(w)))]]
      w <- stack(u, v)
      
      # some checks
      if(nlayers(u) != nlayers(v)) stop("error1")
      if(!all.equal(substr(names(u), 13, 18),
                    substr(names(v), 13, 18))) stop("error2")
      
      rose <- raster::calc(w, fun=windrose, forceapply=TRUE, filename=outfile)
      
      #beginCluster(7)
      #rosefun <- function(x) calc(x, fun=windrose, forceapply=TRUE)
      #rose <- clusterR(w, rosefun, filename=outfile)
      #rose <- clusterR(w, calc, args = list(fun = windrose, forceapply = TRUE))
      #writeRaster(rose, outfile)
      #endCluster()
}


f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
modir <- "data/roses/cfsr_monthly"
lapply(f, cfsr_rose, outdir = modir)


# sum the monthly roses
f <- list.files(modir, full.names=T) %>%
      lapply(stack) %>%
      Reduce("+", .) %>%
      writeRaster("data/roses/cfsr_climatology/windrose_cfsr.tif", overwrite=T)

