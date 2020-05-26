

library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)

select <- dplyr::select




# calculate the spatial distribution of climate analogs for a given location
analogs <- function(climate, # raster stack of first, second timesteps, for 1+ climate variables
                    coords, # lon-lat vector
                    reverse = FALSE, # should reverse analogs be calculated, instead of forward
                    method = "gaussian", # c("gaussian", "triangular", "threshold)
                    sigma = 1 # vector of length 1+: width of similarity kernel; details vary by method
                    # gaussian = standard deviation, triangular = delta of half suitability,
                    # threshold = cutoff
){
      
      if(nlayers(climate) != length(sigma)*2) stop("must provide 2 raster layers and one sigma per climate variable")
      coords <- matrix(coords, ncol=2)
      #a <- climate[[1:length(sigma)]]
      a <- list()
      for(i in 1:length(sigma)){
            sig <- sigma[i]
            clim <- climate[[(i*2-1):(i*2)]]
            if(reverse) clim <- clim[[2:1]]
            target <- raster::extract(clim[[1]], coords)
            if(method == "gaussian") kernel <- function(x) exp(-.5*(x/sig)^2)
            if(method == "triangular") kernel <- function(x) pmax(0, 1 - (abs(x)/sig/2))
            if(method == "threshold") kernel <- function(x) as.integer(abs(x) <= sig)
            
            ai <- calc(clim[[2]] - target, kernel)
            a <- c(a, ai)
      }
      a <- stack(a)
      if(nlayers(a) > 1) a <- prod(a)
      if(class(a) == "RasterStack") a <- a[[1]]
      a
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

add_coords <- function(windrose){
      rows <- cols <- windrose[[1]]
      rows[] <- rep(1:nrow(rows), each=ncol(rows))
      cols[] <- rep(1:ncol(rows), nrow(rows))
      windrose <- stack(windrose, rows, cols)
      names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
      return(windrose)
}


woc <- function(x, y, windrose, climate,
                radius = 1000, time_conv=identity,
                method = "gaussian", sigma = 2,
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
      downtrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="downwind")
      uptrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="upwind")
      clim <- climate %>% crop(circle) %>% mask(circle)
      
      # calculate wind catchment and climate analog surfaces
      s <- list(wind_fwd = accCost(downtrans, origin) %>% "/"(3600) %>% calc(time_conv),
                wind_rev = accCost(uptrans, origin) %>% "/"(3600) %>% calc(time_conv),
                clim_fwd = analogs(clim, coords, method=method, sigma=sigma),
                clim_rev = analogs(clim, coords, method=method, sigma=sigma, reverse = T)) %>%
            stack()
      
      # climate similarity over water is 0
      n_cells <- length(na.omit(values(s$clim_fwd)))
      s[is.na(values(s))] <- 0
      s <- mask(s, circle)
      
      # overlap between wind and climate
      s$overlap_fwd <- s$wind_fwd * s$clim_fwd
      s$overlap_rev <- s$wind_rev * s$clim_rev
      
      if(output == "rasters") return(s)
      
      
      ### summary statistics of wind and climate surfaces
      
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
      
      ss <- s %>% as.list() %>% lapply(ws_summarize, origin=origin)
      for(i in 1:length(ss)) names(ss[[i]]) <- paste0(names(s)[i], "_", names(ss[[i]]))
      ss <- unlist(ss)
      
      if(output == "summary") return(c(x=x, y=y, ss,
                                       runtime = difftime(Sys.time(), start, units="secs")))
      stop("invalid output argument")
}

