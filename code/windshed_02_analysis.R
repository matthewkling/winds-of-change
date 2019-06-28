

# to do: 
#     get gdistance wrapping on a sphere



library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)

select <- dplyr::select

analogs <- function(climate, # raster stack of first, second timesteps
                    coords, # lon-lat vector
                    reverse = FALSE, # should reverse analogs be calculated, instead of forward
                    sigma = 1 # standard deviation of similarity kernel, in degrees Celsius
){
      if(reverse) climate <- climate[[2:1]]
      coords <- matrix(coords, ncol=2)
      target <- raster::extract(climate[[1]], coords)
      kernel <- function(x) exp(-.5*(x/sigma)^2)
      calc(climate[[2]] - target, kernel)
}


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


################


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

# width / height ratio of lat-lon grid cell at a given latitude
distortion <- function(lat, inc=.01){
      distRhumb(c(0, lat), c(inc, lat)) / distRhumb(c(0, lat), c(0, lat+inc))
      #distGeo(c(0, lat), c(inc, lat)) / distRhumb(c(0, lat), c(0, lat+inc))
}


# reallocate weight among windrose directions, to correct for assumption
# that diagonal neighbors are actually at 45 degrees
reallocate2 <- function(x, latitude, ...){
      #x <- rep(1, 8)
      #latitude <- 70
      #browser()
      
      theta <- atan(distortion(latitude, ...)) / pi * 180
      delta <- 45 - theta
      transfer <- delta / (delta + 45)
      
      x[c(1, 3, 5, 7)] <- x[c(1, 3, 5, 7)] * (1 - transfer)
      x[2] <- x[2] + sum(x[c(1, 3)]) * transfer
      x[6] <- x[6] + sum(x[c(5, 7)]) * transfer
      x
}


geo_correct <- function(x){
      
      ## correction based on latitude
      lat <- x[[1]]
      lat[] <- coordinates(lat)[,2]
      rat <- lat
      ratio <- sapply(lat[,1], distortion)
      rat[] <- 1 / rep(ratio, each = ncol(lat))
      
      # extent to which correction applies to each layer of windrose
      for(i in c(2, 6)) x[[i]] <- x[[i]] * rat
      
      # correct for weights allocation
      y <- calc(stack(x, lat), 
                function(x) reallocate2(x[1:8], x[9], inc=.31))
      
      return(x)
}

wind_trans2 <- function(windrose, correction="c"){
      windrose <- geo_correct(windrose) # moved to global
      windrose <- add_coords(windrose)
      trans <- transition_stack(windrose, diffuse2, directions=8, symm=F)
      #geoCorrection(trans, type=correction)
      trans
}


#################


ws_summarize2 <- function(x, # raster layer of wind flow (where positive values are more accessible)
         origin # coordinates of center point (2-column matrix)
){
      browser()
      # transform everything to latlong
      p <- rasterToPoints(x)
      xy <- as.data.frame(rbind(origin, p[,1:2]))
      coordinates(xy) <- c("x", "y")
      crs(xy) <- crs(x)
      latlong <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      xy <- spTransform(xy, latlong)
      xy <- coordinates(xy)
      origin <- xy[1,]
      xy <- xy[2:nrow(xy),]
      
      # size of windshed
      w <- p[,3]
      ws_size <- mean(w)
      
      # distance and bearing to centroid of windshed
      ctd <- apply(xy, 2, function(z) weighted.mean(z, w, na.rm=T)) # not geodesically ideal
      ctd_dist <- distGeo(origin, ctd) / 1000
      ctd_brng <- bearing(origin, ctd)
      
      # distance and bearing to every cell
      dist <- distGeo(origin, xy) / 1000
      brng <- bearing(origin, xy)
      
      # mean distance, mean bearing
      ws_dist <- weighted.mean(dist, w, na.rm=T)
      iso <- circ_sd(brng, w, na.rm=T)
      ws_brng <- iso["bearing"]
      ws_iso <- 1 - anisotropy(brng, w)
      
      # convert centroid back to origin proj
      ctd <- as.data.frame(matrix(ctd, ncol=2))
      coordinates(ctd) <- c("V1", "V2")
      crs(ctd) <- latlong
      ctd <- spTransform(ctd, crs(x))
      ctd <- coordinates(ctd)
      
      names(ctd) <- names(ws_iso) <- names(ws_brng) <- NULL
      c(centroid_x = ctd[1],
        centroid_y = ctd[2],
        centroid_distance = ctd_dist,
        centroid_bearing = ctd_brng,
        windshed_distance = ws_dist,
        windshed_bearing = ws_brng,
        windshed_isotropy = ws_iso,
        windshed_size=ws_size)
}


woc <- function(x, y, windrose, climate, 
                radius = 1000, cost_to_flow, 
                output = "summary"){
      
      start <- Sys.time()
      coords <- c(x, y)
      origin <- matrix(coords, ncol=2)
      message(paste(coords, collapse=" "))
      
      # constrain analysis to region around focal point (for computatioal speed)
      coords_ll <- SpatialPoints(origin, crs(climate)) %>%
            spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
            coordinates()
      circle <- geo_buffer(coords_ll, width = radius * 1000) %>% # km to m
            spTransform(crs(climate))
      trans <- windrose %>% crop(circle) %>% mask(circle) %>% wind_trans2()
      clim <- climate %>% crop(circle) %>% mask(circle)
      
      # calculate wind and climate analog surfaces
      s <- list(wind_fwd = windshed(trans, coords),
                wind_rev = windshed(trans, coords, upwind = T),
                clim_fwd = analogs(clim, coords),
                clim_rev = analogs(clim, coords, reverse = T))
      s <- stack(s)
      
      # modify wind values
      #cost_to_flow <- function(cost) .25 ^ (cost * 1e7)
      s$wind_fwd <- calc(s$wind_fwd, cost_to_flow) %>% mask(clim[[1]])
      s$wind_rev <- calc(s$wind_rev, cost_to_flow) %>% mask(clim[[1]])
      
      # overlap between wind and climate
      s$overlap_fwd <- s$wind_fwd * s$clim_fwd
      s$overlap_rev <- s$wind_rev * s$clim_rev
      if(output == "rasters") return(s)
      
      # summary statistics of wind and climate surfaces
      ss <- s %>% as.list() %>% lapply(ws_summarize2, origin=origin)
      for(i in 1:length(ss)) names(ss[[i]]) <- paste0(names(s)[i], "_", names(ss[[i]]))
      ss <- unlist(ss)
      
      # total overlap across landscape 
      # (mean rather than sum, bc landscapes differ in number of cells due to geodesy)
      if(output == "summary") return(c(x=x, y=y, ss,
                                       runtime = difftime(Sys.time(), start, units="secs")))
      stop("invalid output argument")
}


############

# to_equal_area <- function(r, ...){
#       projectRaster(r, crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), ...)}

land <- raster("f:/cfsr/land.tif") %>% 
      rotate()# %>% to_equal_area(method="ngb")

# load windrose data
rose <- stack("data/roses_force/cfsr_climatology/roses_cfsr_1980s.tif") %>%
      rotate()# %>% to_equal_area()

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)

# correct geodesy
#rose <- geo_correct(rose)



# mean temperature (kelvin)
# climate <- list.files("data/cfsr_monthly/tmp2m", full.names=T)[1:24] %>%
#       stack() %>%
#       mean() %>%
#       rotate() %>%
#       #to_equal_area() %>%
#       mask(land)
#climate <- stack(climate, climate + 3)

climate <- stack("data/geographic/processed/temperature.tif")


# function to convert wind cost values to flow values
#cost_to_flow <- function(cost) .5 ^ (cost * 1e8)
cost_to_flow <- function(cost) (1/cost) ^ (1/3)

# test
coords <- c(-110, 45)
x <- woc(coords[1], coords[2], rose, climate, cost_to_flow=cost_to_flow)



ncores <- 6
pixels <- land %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(x < 1.5e7,
             y > -7.75e6) %>%
      #sample_n(1000) %>%
      mutate(batch = sample(1:ncores, nrow(.), replace=T)) %>%
      split(.$batch)


require(doParallel)
cl <- makeCluster(ncores)
registerDoParallel(cl)
start <- Sys.time()
d <- foreach(pts = pixels,
             .combine="rbind",
             .packages=(.packages())) %dopar% {
                   map2(pts$x, pts$y, possibly(woc, NULL), 
                        windrose=rose, climate=climate, cost_to_flow=cost_to_flow, radius=500) %>%
                         do.call("rbind", .) %>%
                         as.data.frame()
             }
Sys.time() - start
stopCluster(cl)
write_csv(d, "data/windshed/force_500km_v3.csv")


stop("wootwoot")




### some exploratory plots ###

pts <- pixels[[1]][1:49,]
r <- map2(pts$x, pts$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, cost_to_flow=cost_to_flow, radius=500,
          output="rasters") %>%
      lapply(function(x) x$wind_fwd)
po <- map2(pts$x, pts$y, function(x, y) matrix(c(x, y), ncol=2))
ri <- map2(r, po, ws_summarize) %>%
      map_dbl("isotropy") %>%
      round(3)
rd <- map2(r, po, ws_summarize) %>%
      map_dbl("mean_distance") %>%
      round(3)

r <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(r)) r[[i]]$id <- i
for(i in 1:length(r)) r[[i]]$isotropy <- ri[i]
for(i in 1:length(r)) r[[i]]$distance <- rd[i]
r <- bind_rows(r)

ggplot(r %>% mutate(dir = isotropy ^ distance), 
       aes(x, y, fill=wind_fwd)) +
      facet_wrap(~dir, scales="free") +
      geom_raster() +
      scale_fill_viridis_c() +
      theme_void()




