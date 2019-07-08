

library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)

select <- dplyr::select




# calculate the spatial distribution of climate analogs for a given location
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


unwrap <- function(x, width=20){
      x <- x %>% crop(extent(180-width, 180, -90, 90)) %>% shift(-360) %>% merge(x)
      x <- x %>% crop(extent(-180, -180+width, -90, 90)) %>% shift(360) %>% merge(x)
      x
}


geo_circle <- function(pts, width, thresh=100) {
      #pts = coords_ll, width = 500000
      # https://stackoverflow.com/questions/25411251/buffer-geospatial-points-in-r-with-gbuffer
      angles <- seq(from = 0, to = 360, by = 5)
      vertices <- geosphere::destPoint(p = pts, b = angles, d = width)
      
      if(pts[1,1] > thresh) vertices[,1] <- ifelse(vertices[,1] > 0, vertices[,1], vertices[,1] + 360)
      if(pts[1,1] < -thresh) vertices[,1] <- ifelse(vertices[,1] < 0, vertices[,1], vertices[,1] - 360)
      
      poly <- Polygon(vertices, hole=F)
      poly <- Polygons(list(poly), ID=NA)
      poly <- SpatialPolygons(Srl = list(poly), 
                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
      poly
}




ws_summarize2 <- function(x, # raster layer of wind flow (where positive values are more accessible)
                          origin, # coordinates of center point (2-column matrix)
                          latlon = T # is data already in lat-lon coordinates
){
      
      # transform everything to latlong
      p <- rasterToPoints(x)
      p <- p[is.finite(p[,3]), ]
      xy <- p[,1:2]
      
      if(!latlon){
            xy <- as.data.frame(rbind(origin, xy))
            coordinates(xy) <- c("x", "y")
            crs(xy) <- crs(x)
            latlong <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
            xy <- spTransform(xy, latlong)
            xy <- coordinates(xy)
            origin <- xy[1,]
            xy <- xy[2:nrow(xy),]
      }
      
      # size of windshed
      w <- round(p[,3], 3)
      ws_size <- mean(w)
      
      # distance and bearing to centroid of windshed
      ctd <- apply(xy, 2, function(z) weighted.mean(z, w, na.rm=T)) # not geodesically ideal
      ctd_dist <- distGeo(origin, ctd) / 1000
      ctd_brng <- bearing(origin, ctd)
      
      # distance and bearing to every cell
      dist <- distGeo(origin, xy) / 1000
      brng <- round(bearing(origin, xy))
      
      # mean distance, mean bearing
      ws_dist <- weighted.mean(dist, w, na.rm=T)
      iso <- circ_sd(brng, w, na.rm=T)
      ws_brng <- iso["bearing"]
      #ws_iso <- 1 - anisotropy(brng, w) # too slow
      ws_iso <- iso["iso"]
      
      # convert centroid back to origin proj
      if(!latlon){
            ctd <- as.data.frame(matrix(ctd, ncol=2))
            coordinates(ctd) <- c("V1", "V2")
            crs(ctd) <- latlong
            ctd <- spTransform(ctd, crs(x))
            ctd <- coordinates(ctd)
      }
      
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
      
      # constrain analysis to region around focal point
      coords_ll <- SpatialPoints(origin, crs(climate)) %>%
            spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
            coordinates()
      circle <- geo_circle(coords_ll, width = radius * 1000) # km to m
      
      # prepare wind and climate datasets       
      wind <- windrose %>% crop(circle) %>% mask(circle) %>% add_coords() 
      downtrans <- wind %>% transition_stack(downwind, directions=8, symm=F)
      uptrans <- wind %>% transition_stack(upwind, directions=8, symm=F)
      clim <- climate %>% crop(circle) %>% mask(circle)
      
      # calculate wind catchment and climate analog surfaces
      s <- list(wind_fwd = accCost(downtrans, origin),
                wind_rev = accCost(uptrans, origin),
                clim_fwd = analogs(clim, coords),
                clim_rev = analogs(clim, coords, reverse = T)) %>%
            stack()
      
      # modify wind values
      s$wind_fwd <- calc(s$wind_fwd, cost_to_flow) %>% mask(clim[[1]])
      s$wind_rev <- calc(s$wind_rev, cost_to_flow) %>% mask(clim[[1]])
      
      # overlap between wind and climate
      s$overlap_fwd <- s$wind_fwd * s$clim_fwd
      s$overlap_rev <- s$wind_rev * s$clim_rev
      if(output == "rasters") return(s)
      
      # summary statistics of wind and climate surfaces
      ex <- extent(s)
      if(ex@xmin < -180){
            shft <- -180 - ex@xmin
            s <- shift(s, shft) # ws_summarize needs real geography
            origin[1,1] <- origin[1,1] + shft
      }
      if(ex@xmax > 180){
            shft <- 180 - ex@xmax
            s <- shift(s, shft)
            origin[1,1] <- origin[1,1] + shft
      }
      
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


land <- raster("f:/cfsr/land.tif") %>% 
      rotate() %>% unwrap(180)

# load windrose data
rose <- stack("data/windrose/windrose_p2_2000s.tif") %>%
      rotate() %>% unwrap(180)

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)

# load current/future climate data
climate <- stack("data/geographic/processed/temperature.tif") %>% unwrap(180)


# function to convert wind cost values to flow values
#cost_to_flow <- function(cost) .5 ^ (cost * 1e8)
cost_to_flow <- function(cost) (1/cost) ^ (1/2)

# test
coords <- c(179, 67)
x <- woc(coords[1], coords[2], rose, climate, cost_to_flow=cost_to_flow)
#profvis({ x <- woc(coords[1], coords[2], rose, climate, cost_to_flow=cost_to_flow) })

ncores <- 6
pixels <- land %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(x < 180, x > -180) %>%
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
write_csv(d, "data/windshed/p2_500km.csv")


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




