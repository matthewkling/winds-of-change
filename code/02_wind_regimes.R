

# derive wind regime stats for all grid cells 

library(windscape)
library(tidyverse)
library(raster)


regime <- function(x, ...){
      uv <- matrix(x, ncol=2, byrow=F)
      u <- mean(uv[,1])
      v <- mean(uv[,2])
      
      speeds <- sqrt(uv[,1]^2 + uv[,2]^2)
      speed <- mean(speeds)
      dir <- do.call(atan2, as.list(colMeans(uv))) / pi * 180
      
      dirs <- atan2(uv[,1], uv[,2])
      x <- sum(sin(dirs) * speeds) / sum(speeds)
      y <- sum(cos(dirs) * speeds) / sum(speeds)
      
      rbar <- sqrt(x^2 + y^2) # mean resultant vector
      iso <- sqrt(1-rbar) # circular standard deviation
      
      c(speed = speed, dir = dir, 
        rbar = rbar, aniso = 1 - iso,
        u = u, v = v, x = x, y = y)
}

windregime <- function(infile, outdir, ncores=1){
      
      # file admin
      message(infile)
      outfile <- paste0(outdir, "/", basename(infile))
      outfile <- sub("grb2", "tif", outfile)
      if(any(file.exists(c(outfile,
                           sub("\\.tif", paste0("_", ncores, ".tif"), outfile))))){
            return("skipping")}
      
      # load and organize hourly wind data
      inv <- paste0(infile, ".inv") %>%
            readLines()
      process <- which(!grepl("6 hour fcst", inv))
      w <- infile %>%
            brick() %>%
            subset(process)
      even <- function(x) x %% 2 == 0
      w <- w[[c(which(!even(1:nlayers(w))),
                which(even(1:nlayers(w))))]]
      
      if(ncores == 1){
            wr <- raster::calc(w, fun=regime, forceapply=TRUE, filename = outfile)
            return(wr)
      } else {
            
            # split dataset into batches
            nlayer <- nlayers(w)/2
            core <- rep(1:ncores, each=floor(nlayer/ncores))
            core <- c(core, rep(ncores, nlayer %% ncores))
            w <- lapply(1:ncores, function(x) list(i=x,
                                                   data=w[[c(which(core == x),
                                                             which(core == x) + nlayer)]]))
            
            # process batches in parallel
            require(doParallel)
            cl <- makeCluster(ncores)
            registerDoParallel(cl)
            wr <- foreach(x = w,
                          .export = "regime",
                          .packages=c("raster", "windscape")) %dopar% {
                                raster::calc(x$data, fun=regime, forceapply=TRUE, 
                                             filename=paste0(substr(outfile, 1, nchar(outfile)-4),
                                                             "_", x$i,
                                                             substr(outfile, nchar(outfile)-3, nchar(outfile))))
                          }
            stopCluster(cl)
      }
}

if(F){
      # hourly input data
      f <- list.files("f:/CFSR/wnd10m", full.names=T)
      f <- f[!grepl("inv", f)]
      f <- rev(f)
      map(f, windregime, outdir="data/wind_regime/monthly", ncores=6)
}


# compile long-term regimes from monthlies

mean_regime <- function(x){
      
      # average speed
      speed <- lapply(x, function(y) stack(y)[[1]]) %>%
            Reduce("+", .) %>%
            "/"(length(x))
      
      # prevailing direction
      u <- lapply(x, function(y) stack(y)[[5]]) %>%
            Reduce("+", .) %>%
            "/"(length(x))
      v <- lapply(x, function(y) stack(y)[[6]]) %>%
            Reduce("+", .) %>%
            "/"(length(x))
      dir <- u
      dir[] <- atan2(values(u), values(v)) / pi * 180
      
      # anisotropy (mean of partial mean resultant vectors, weighted by speed of each)
      ss <- lapply(x, function(y) stack(y)[[1]]) %>%
            Reduce("+", .)
      xs <- lapply(x, function(y) stack(y)[[7]] * stack(y)[[1]]) %>%
            Reduce("+", .) %>%
            "/"(ss)
      ys <- lapply(x, function(y) stack(y)[[8]] * stack(y)[[1]]) %>%
            Reduce("+", .) %>%
            "/"(ss)
      rbar <- sqrt(xs^2 + ys^2)
      aniso <- 1 - sqrt(1 - rbar)
      
      stack(speed, dir, aniso, u, v)
}

if(F){
      wr <- list.files("data/wind_regime/monthly", full.names=T) %>%
            mean_regime() %>%
            writeRaster("data/wind_regime/regimes.tif", overwrite=T)
}


