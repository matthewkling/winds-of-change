


# to do: 
#     get gdistance wrapping on a sphere
#     downweight oceanic conductance







library(windshed)
library(tidyverse)
library(raster)
library(gdistance)




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

woc <- function(x, y, windrose, climate, radius = 1000, output = "summary"){
      
      start <- Sys.time()
      coords <- c(x, y)
      message(paste(coords, collapse=" "))
      
      # constrain analysis to region around focal point (for computatioal speed)
      coords_ll <- SpatialPoints(matrix(coords, ncol=2), crs(climate)) %>%
            spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
            coordinates()
      circle <- geo_buffer(coords_ll, width = radius * 1000) %>% # km to m
            spTransform(crs(climate))
      trans <- windrose %>% crop(circle) %>% mask(circle) %>% wind_trans()
      clim <- climate %>% crop(circle) %>% mask(circle)
      
      # calculate wind and climate analog surfaces
      s <- list(wind_fwd = windshed(trans, coords),
                wind_rev = windshed(trans, coords, upwind = T),
                clim_fwd = analogs(clim, coords),
                clim_rev = analogs(clim, coords, reverse = T))
      s <- stack(s)
      
      # modify wind values
      transform <- function(x) .05 ^ (x * 1e9)
      s$wind_fwd <- calc(s$wind_fwd, transform) %>% mask(clim[[1]])
      s$wind_rev <- calc(s$wind_rev, transform) %>% mask(clim[[1]])
      
      # overlap between wind and climate
      s$overlap_fwd <- s$wind_fwd * s$clim_fwd
      s$overlap_rev <- s$wind_rev * s$clim_rev
      if(output == "rasters") return(s)
      
      # total overlap across lanscape 
      # (mean rather than sum, bc landscapes differ in size due to geodesy)
      if(output == "summary") return(c(x=x, y=y, 
                                       cellStats(s, "mean"),
                                       runtime = difftime(Sys.time(), start, units="secs")))
      stop("invalid output argument")
}


############

to_equal_area <- function(r, ...){
      projectRaster(r, crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), ...)}

land <- raster("f:/cfsr/land.tif") %>% 
      rotate() %>% 
      to_equal_area(method="ngb")

# load windrose data
windrose <- stack("data/roses/cfsr_climatology/windrose_cfsr.tif") %>%
      rotate() %>%
      to_equal_area()
#trans <- wind_trans(windrose)

# mean temperature (kelvin)
climate <- list.files("data/cfsr_monthly/tmp2m", full.names=T)[1:24] %>%
      stack() %>%
      mean() %>%
      rotate() %>%
      to_equal_area() %>%
      mask(land)
climate <- stack(climate, climate + 3)

# test
#coords <- c(-118, 38)
coords <- c(-5e6, -2e6)
woc(coords[1], coords[2], windrose, climate)


# sample of land pixels
pts <- land %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(x < 1.5e7,
             y > -7.75e6) %>%
      sample_n(1000)

d <- map2(pts$x, pts$y, possibly(woc, NULL), 
          windrose=windrose, climate=climate, radius=500) %>%
      do.call("rbind", .) %>%
      as.data.frame()


# estimated hours to run all, on one core, at 500km radius
mean(d$runtime) * length(na.omit(values(land))) / 60 / 60


ggplot(d, aes(overlap_fwd, overlap_rev)) +
      geom_point()

ggplot(d, aes(overlap_fwd, overlap_rev)) +
      geom_point() +
      scale_x_log10() + scale_y_log10()

ggplot(d, aes(wind_fwd, overlap_fwd, color=clim_fwd)) +
      geom_point() +
      scale_color_viridis_c(trans="log10") +
      scale_x_log10() + scale_y_log10()

ggplot(d, aes(x, y, color=overlap_fwd)) +
      geom_point() +
      scale_color_viridis_c(trans="log10")


ggplot(d %>% filter(x < 1.5e7), aes(x, y, color=sqrt(runtime))) +
      geom_point() +
      scale_color_viridis_c(trans="log10")

d %>% filter(x < 1.5e7,
             y > -7.75e6) %>% 
      ggplot(aes(y, runtime)) + geom_point()




# get continentatlity and elevation covariates
# get biome covariates


