

# create a probability stack indicating the velocity-weighted frequency
# with which the wind blows in each semi-cardinal direction


# velocity from u&v components
velocity <- function(x) sqrt(sum(x^2))

# force is proportional to the cube of velocity
force <- function(x) velocity(x) ^3 

# direction from u&v components -- in degrees, clockwise from 12:00
direction <- function(x, y) atan2(x, y) * 180 / pi

spin90 <- function(x){
      x <- x - 90
      x[x<(-180)] <- x[x<(-180)] + 360
      x
}

# convert an angle [0-45] to a loading [0-1] in the 45-degree direction
loading <- function(theta){ # theta in [0,45]
      sin((4*(theta-22.5))*pi/180)/2+.5
}

#velocity-weighted frequency of wind in each quadrant
windrose <- function(x){
      #x0 <- s[277,349]
      
      # restructure: row=timestep, col=u&v components
      m <- matrix(x, ncol=2, byrow=F)
      
      # calcualte force and direction
      frc <- apply(m, 1, force) # fun should be force or velocity
      dir <- apply(m, 1, function(x) spin90(direction(x[2], -1*x[1])))
      
      # octant bounds for each vector 
      d1 <- plyr::round_any(dir, 45, floor)
      d1[d1 == -180] <- 180
      d2 <- plyr::round_any(dir, 45, ceiling)
      d2[d2 == -180] <- 180
      
      # loadings on each octant bound, based on angle
      theta <- dir %% 45
      l1 <- loading(theta)
      l2 <- 1-l1
      
      # add some dummy data to ensure all quadrants are represented
      frc <- c(frc, rep(0,8)) 
      l1 <- c(l1, rep(0,8)) 
      l2 <- c(l2, rep(0,8)) 
      d1 <- c(d1, c(-135, -90, -45, 0, 45, 90, 135, 180))
      d2 <- c(d2, c(-135, -90, -45, 0, 45, 90, 135, 180))
      
      # repackage
      d <- c(d1, d2)
      l <- c(l1, l2)
      frc <- c(frc, frc)
      
      # summarize
      d <- split(frc*l, d)
      d <- unlist(lapply(d, sum))
      return(as.numeric(d))
}


roses <- function(indir, outdir, months, cpus){
      
      require(doParallel)
      require(raster)
      
      uf <- list.files(indir, pattern="uwnd", full.names=T)
      vf <- list.files(indir, pattern="vwnd", full.names=T)
      
      #cl <- makeCluster(cpus)
      #registerDoParallel(cl)
      
      windrose <- windrose
      
      #r <- foreach(y = 1:length(uf),
      #             .packages=c("raster")) %dopar% {
      
      r <- c()
      for(y in 1:length(uf)){ 
            
            message(uf[y])
            
            u <- stack(uf[y])
            v <- stack(vf[y])
            if(!all.equal(names(u), names(v))) stop("houston, we have a problem")
            
            month <- substr(names(u), 7,8)
            mosi <- which(month %in% stringr::str_pad(months, 2, "left", 0))
            s <- stack(u[[mosi]], v[[mosi]])
            
            rose <- raster::calc(s, fun=windrose, forceapply=TRUE)
            outfile <- paste0(outdir, "/temp_", y, ".tif")
            writeRaster(rose, outfile, overwrite=T)
            #return(outfile)
            r <- c(r, outfile)
      }
      
      # sum years into final stack 
      s <- lapply(r, stack)
      s <- Reduce("+", s)
      writeRaster(s, paste0(outdir, "/windrose_", months[1], "_", months[length(months)], ".tif"), overwrite=T)
      file.remove(r)
      
      #stopCluster(cl)
}

roses(indir="E:/wind/windsurf/narr_data_raw",
      outdir="E:/wind/windsurf/windrose_data",
      months=2:4,
      cpus=5)


##################################

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
lapply(f, cfsr_rose, outdir = "data/roses/cfsr_monthly")


# sum the monthly roses

