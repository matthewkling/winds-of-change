


library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)

select <- dplyr::select



geo_buffer <- function(pts, width) {
      # https://stackoverflow.com/questions/25411251/buffer-geospatial-points-in-r-with-gbuffer
      angles <- seq(from = 0, to = 360, by = 5)
      vertices <- geosphere::destPoint(p = pts, b = angles, d = width)
      poly <- Polygon(vertices, hole=F)
      poly <- Polygons(list(poly), ID=NA)
      poly <- SpatialPolygons(Srl = list(poly), 
                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
      poly
}

if(F){
      
      ws <- function(x, y, windrose, 
                     radius = 1000){
            
            start <- Sys.time()
            coords <- c(x, y)
            origin <- matrix(coords, ncol=2)
            message(paste(coords, collapse=" "))
            
            # constrain analysis to region around focal point (for computatioal speed)
            coords_ll <- SpatialPoints(origin, crs(windrose)) %>%
                  spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
                  coordinates()
            circle <- geo_buffer(coords_ll, width = radius * 1000) %>% # km to m
                  spTransform(crs(windrose))
            trans <- windrose %>% crop(circle) %>% mask(circle) %>% wind_trans()
            
            # calculate wind and climate analog surfaces
            windshed(trans, coords)
      }
      
      
      ############
      
      to_equal_area <- function(r, ...){
            projectRaster(r, crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), ...)}
      
      land <- raster("f:/cfsr/land.tif") %>% 
            rotate()# %>% to_equal_area(method="ngb")
      
      # load windrose data
      rose <- stack("data/roses_force/cfsr_climatology/roses_cfsr_2000s.tif") %>%
            rotate()# %>% to_equal_area()
      
      rose %>% mask(land) %>% values() %>% na.omit() %>% mean()
      rose[] <- 25
      
}






diffuse2 <- function(x){
      p <- x[c(1,3,5,7,9,11,13,15)] # SW ...cw... S
      g <- x[17:20]
      if(g[1]<g[2] & g[4]<g[3]) return(p[1]/sqrt(2)) #SW
      if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
      if(g[2]<g[1] & g[4]<g[3]) return(p[3]/sqrt(2)) #NW
      if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
      if(g[2]<g[1] & g[3]<g[4]) return(p[5]/sqrt(2)) #NE
      if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
      if(g[1]<g[2] & g[3]<g[4]) return(p[7]/sqrt(2)) #SE
      if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
      # dividing by 2 root 2 prevents the geocorrection from distorting the probability field
}


diffuse3 <- function(x){
      p <- x[c(1,3,5,7,9,11,13,15)] # SW ...cw... S
      g <- x[17:20]
      if(g[1]<g[2] & g[4]<g[3]) return(p[1]) #SW
      if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
      if(g[2]<g[1] & g[4]<g[3]) return(p[3]) #NW
      if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
      if(g[2]<g[1] & g[3]<g[4]) return(p[5]) #NE
      if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
      if(g[1]<g[2] & g[3]<g[4]) return(p[7]) #SE
      if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
      # dividing by 2 root 2 prevents the geocorrection from distorting the probability field
}

# width / height ratio of lat-lon grid cell at a given latitude
distortion <- function(lat, inc=.01){
      #distRhumb(c(0, lat), c(inc, lat)) / distRhumb(c(0, lat), c(0, lat+inc))
      distGeo(c(0, lat), c(inc, lat)) / distRhumb(c(0, lat), c(0, lat+inc))
}


# reallocate weight among windrose directions, to correct for assumption
# that diagonal neighbors are actually at 45 degrees
reallocate <- function(x, latitude, ...){
      #x <- rep(1, 8)
      #latitude <- 70
      #browser()
      theta <- atan(distortion(latitude, ...)) / pi * 180
      
      got <- 45 / 2 * .75
      shd <- theta / 2 * .75
      #ew_shd <- (90 - theta) / 2 * .75
      
      transfer <- (got - shd) / got
      trans_n <- x[4] * transfer
      trans_s <- x[8] * transfer
      x[c(4, 8)] <- x[c(4, 8)] * (1 - transfer)
      x[c(2, 6)] <- x[c(2, 6)] + (trans_n + trans_s)/2
      
      
      x
}

reallocate2 <- function(x, latitude, 
                        ...){
      #x <- rep(1, 8)
      #latitude <- 70
      #if(aspect > 1.5) browser()
      theta <- atan(distortion(latitude, ...) ) / pi * 180
      delta <- 45 - theta
      transfer <- delta / (delta + 45)
      
      x[c(1, 3, 5, 7)] <- x[c(1, 3, 5, 7)] * (1 - transfer)
      x[2] <- x[2] + sum(x[c(1, 3)]) * transfer
      x[6] <- x[6] + sum(x[c(5, 7)]) * transfer
      x
}


reallocate3 <- function(x, latitude, aspect, 
                        ...){
      #x <- rep(1, 8)
      #latitude <- 70
      #browser()
      theta <- atan(distortion(latitude, ...) * aspect) / pi * 180
      #delta <- 45 - theta
      transfer <- (theta - 45)/(90 + theta)
      #if(sample(100, 1) == 1) message(paste("     ", theta))
      x[c(1, 3, 5, 7)] <- x[c(1, 3, 5, 7)] * (1 - transfer)
      x[2] <- x[2] + sum(x[c(1, 3)]) * transfer
      x[6] <- x[6] + sum(x[c(5, 7)]) * transfer
      x
}

reallocate4 <- function(x, aspect, 
                        ...){
      #x <- rep(1, 8)
      #latitude <- 70
      #browser()
      #aspect <- aspect ^ .5
      theta <- atan(1 / aspect) / pi * 180
      
      trans <- theta / 90 - .5
      
      x[c(1, 3, 5, 7)] <- x[c(1, 3, 5, 7)] * (1 - abs(trans))
      if(trans > 0) {
            x[2] <- x[2] + sum(x[c(1, 3)]) * trans
            x[6] <- x[6] + sum(x[c(5, 7)]) * trans
      }     
      if(trans < 0) {
            x[4] <- x[4] + sum(x[c(3, 5)]) * trans
            x[8] <- x[8] + sum(x[c(1, 7)]) * trans
      } 
      
      
      x[c(2, 6)] <- x[c(2, 6)] / aspect
      #x[c(4, 8)] <- x[c(4, 8)] * aspect
      x[c(1, 3, 5, 7)] <- x[c(1, 3, 5, 7)] / (aspect)
      x
}



geo_correct <- function(x){
      #browser()
      
      # downweight diagonal travel in proportion to grid distance
      for(i in c(1,3,5,7)) x[[i]] <- x[[i]] / sqrt(2)
      
      
      # calculate aspect ratio (km wide/high) for each cell,
      # which is a function of latitude and raster resolution
      
      lat <- x[[1]]
      lat[] <- coordinates(lat)[,2]
      lat <- sapply(lat[,1], distortion)
      
      rsn <- res(x)
      rsn <- rsn[1] / rsn[2]
      
      aspect <- x[[1]]
      aspect[] <- rep(lat * rsn, each=ncol(aspect))
      
      
      
      # re-balance neighbor weights to reflect aspect ratio
      x <- calc(stack(x, aspect), function(x) reallocate4(x[1:8], x[9]) )
      
      
      
      
      # extent to which correction applies to each layer of windrose
      #easting <- abs(sin(windrose_bearings() / 180 * pi))
      ###for(i in c(1,3,5,7)) x[[i]] <- x[[i]] * (1 + (rat - 1) * sin(pi/4))
      ###for(i in c(1:3, 5:7)) x[[i]] <- x[[i]] * rat
      
      #
      
      # note: work region. trying to get reallocate2 to accomodate the 
      # jumps in cell resolution with latitude.
      # adjust that funciton. make sure it doesn't assume aspect <1
      
      
      
      # correct for weights allocation
      # x <- calc(stack(x, lat), function(x) reallocate2(x[1:8], x[9],
      #                                                  inc=mean(res(rose))) )
      
      
      #for(i in c(2, 6)) x[[i]] <- x[[i]] * rat
      #for(i in c(1,3,5,7)) x[[i]] <- x[[i]] * (1 + (rat - 1) * sin(pi/4))
      
      

      return(x)
}

wind_trans2 <- function(windrose, correction="c"){
      windrose <- geo_correct(windrose)
      windrose <- add_coords(windrose)
      trans <- transition_stack(windrose, diffuse3, directions=8, symm=F)
      #geoCorrection(trans, type=correction)
      trans
}

ws <- function(x, y, windrose, 
               radius = 1000){
      
      start <- Sys.time()
      coords <- c(x, y)
      origin <- matrix(coords, ncol=2)
      message(paste(coords, collapse=" "))
      
      # constrain analysis to region around focal point (for computatioal speed)
      coords_ll <- SpatialPoints(origin, crs(windrose)) %>%
            spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
            coordinates()
      circle <- geo_circle(coords_ll, width = radius * 1000)
      
      agg <- round(1/distortion(y))
      trans <- windrose %>% crop(circle) %>% aggregate(c(agg, 1)) %>% mask(circle) %>% wind_trans2()
      
      # calculate wind and climate analog surfaces
      wshd <- windshed(trans, coords)
      
      dst <- wshd
      dst[] <- pointDistance(coords, coordinates(dst), lonlat=T) / 1000
      stack(wshd, dst)
      
}

rose[] <- 5
rose2 <- rose
rose2[] <- 1
rose2[[2]][] <- 2
rose2[[8]][] <- 2
rose2[[3]][] <- 3
rose2[[7]][] <- 3
rose2[[4]][] <- 5
rose2[[6]][] <- 5
rose2[[5]][] <- 8

lats <- seq(0, 85, length.out=16)

d <- lats %>%
      lapply(function(x) ws(0, x, rose, 500)) %>%
      lapply(rasterToPoints) %>%
      lapply(as.data.frame) %>%
      lapply(function(x) mutate(x, id=y[1])) %>%
      bind_rows() %>%
      rename(wind = layer.1, dist = layer.2)


p <- ggplot(d %>% filter(is.finite(wind)), aes(x, y, fill=wind)) +
      facet_wrap(~id, scales="free") +
      geom_raster() +
      scale_fill_gradientn(colors=c("black", "white", "darkgreen", "yellow", 
                                    "darkred", "cyan", "darkblue", "black")) +
      theme_void()
p

stop("okay okay okayyyyy")
###

p <- ggplot(d, aes(dist, wind, color=factor(id))) +
      geom_smooth(se=F)

betas <- d %>% filter(is.finite(wind)) %>%
      split(.$id) %>%
      map(~ lm(wind ~ dist, data = .)) %>%
      map(coef) %>% map_dbl("dist") %>%
      data.frame(beta = .,
                 lat = lats)
ggplot(betas, aes(lat, beta)) + geom_point()




s <- d %>% group_by(id) %>% filter(is.finite(wind)) %>% summarize(wind=mean(wind)) %>% mutate(rose=1)
s2 <- d %>% group_by(id) %>% filter(is.finite(wind)) %>% summarize(wind=mean(wind)) %>% mutate(rose=2)

s <- bind_rows(s, s2)
ggplot(s %>% mutate(rose=paste0("r", rose)) %>% spread(rose, wind), 
       aes(r1, r2, color=id)) +
      geom_point(size=5) +
      scale_color_viridis_c()
      
      
# geo_correct the whole global windrose, not each landscape



# reallocate angular winds to E-W cardinals.
# assume half came from NS cardinals, reallocate the other half



