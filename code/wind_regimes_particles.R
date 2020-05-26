
# load libraries
library(raster)
library(tidyverse)

# load raster dataset with 2 layers: zonal and meiridional windspeeds
w <- stack("data/regimes.tif") %>% stack(raster("data/land.tif")) %>% rotate()
names(w) <- c("speed", "direction", "anisotropy", "u", "v", "land")


# latitude-based adjustment factor to correct for geodesic distortion
distortion <- function(lat){
      xd <- geosphere::distGeo(c(0, lat), c(.0001, lat))
      yd <- geosphere::distGeo(c(0, lat), c(0, lat - sign(lat) * .0001))
      yd / xd
}



# to get winkel-tripel to work, borrowing heavily from 
# https://wilkelab.org/practicalgg/articles/Winkel_tripel.html

library(tidyverse)
library(sf)        # for manipulation of simple features objects
library(lwgeom)    # for st_transform_proj()
library(rworldmap) # for getMap()

world <- getMap(resolution = "low") %>% rgeos::gBuffer(width = 0) %>% rgeos::gUnaryUnion()
world_sf <- st_as_sf(world)

crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
world_wintri <- st_transform_proj(world_sf, crs = crs_wintri)

grat_wintri <- 
      st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
      st_transform_proj(crs = crs_wintri)

# vectors of latitudes and longitudes that go once around the 
# globe in 1-degree steps
lats <- c(90:-90, -90:90, 90)
longs <- c(rep(c(180, -180), each = 181), 180)

# turn into correctly projected sf collection
wintri_outline <- 
      list(cbind(longs, lats)) %>%
      st_polygon() %>%
      st_sfc( # create sf geometry list column
            crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      ) %>% 
      st_sf() %>%
      st_transform_proj(crs = crs_wintri) # transform to Winkel tripel


######## speed-based sampling ########

particle_trails <- function(d, # vector field rasters
                            np, # number of particles
                            nt, # number of timesteps
                            r, # smoothing neighborhood radius
                            scale){ # multiplier converting velocity to lat-lon
      #browser()
      x <- sort(unique(coordinates(d)[,1]))
      y <- desc(sort(unique(coordinates(d)[,2])))
      dx <- as.matrix(d[[1]])
      dy <- as.matrix(d[[2]])
      distort <- sapply(y, distortion)
      
      
      # sample starting locations, in proportion to area and 1/speed
      speed <- sqrt(sum(d^2))
      lat <- speed
      lat[] <- coordinates(lat)[,2]
      weight <- cos(lat/180*pi) / speed
      start <- coordinates(d)[sample(1:ncell(speed), np, replace=T, 
                                     prob = values(weight)),]
      
      # add some noise
      rs <- mean(res(d)) / 2
      start <- start +  runif(nrow(start)*2, -rs, rs)
      
      ptx <- pty <- matrix(NA, nt, np)
      for(p in 1:np){
            #px <- runif(1, min(x), max(x))
            #py <- runif(1, min(y), max(y))
            #py <- sample(Y, 1, replace=T, prob=Yw)
            px <- start[p,1]
            py <- start[p,2]
            
            for(t in 1:nt){
                  xi <- which(x > px - r & x < px + r)
                  yi <- which(y > py - r & y < py + r)
                  if(length(xi)==0 | length(yi)==0) browser()
                  px <- px + mean(dx[yi, xi]) * scale * mean(distort[yi])
                  py <- py + mean(dy[yi, xi]) * scale
                  ptx[t, p] <- px
                  pty[t, p] <- py
                  
                  if(py > 90 | py < -90) break()
                  if(px > 180) px <- px - 360
                  if(px < -180) px <- px + 360
                  if(px > 180 | px < -180) break()
            }
      }
      
      
      frame <- function(pd){
            colnames(pd) <- paste0("p", 1:ncol(pd))
            as.data.frame(pd) %>%
                  mutate(time = 1:nrow(.)) %>%
                  gather(id, z, -time)
      }
      
      left_join(frame(ptx) %>% rename(x = z),
                frame(pty) %>% rename(y = z))
}


#f <- particle_trails(d = w[[c("u", "v")]], np = 100000, nt = 50, r = 1, scale = .1)
#saveRDS(f, "data_speedsampling.rds")

