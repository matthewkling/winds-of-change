
library(windscape)
library(tidyverse)
library(raster)
library(gdistance)


# build windrose based on uniform distribution of wind vectors
if(F){
      f <- list.files("f:/CFSR/wnd10m", full.names=T)
      f <- f[!grepl("inv", f)]
      if(F) r <- raster(f[1])
      
      a <- r # %>% aggregate(2)
      a <- crop(a, extent(80, 90, 0, 10))
      s <- a
      s[] <- 0
      # a uniform distribution of wind angles
      for(i in seq(0, 359, 1)){
            a[] <- i
            s <- stack(s, a)
      }
      s <- s[[2:nlayers(s)]]
      
      # u,v components at 5 m/s
      u <- calc(s, function(x) sin(x / 180 * pi)) * 5
      v <- calc(s, function(x) cos(x / 180 * pi)) * 5
      w <- stack(u, v)
      
      p <- 1
      wr <- w %>%
            add_res() %>%
            add_lat() %>%
            raster::calc(fun=function(x) windrose_geo(x, p=p), 
                         forceapply=TRUE)
      
      # units are currently summed conductances in 1/sec^p
      # divide by number of timesteps, to give mean conductance,
      # and then convert back to 1/s
      wr <- wr / nlayers(u)
      wr <- wr ^ (1/p)
      names(wr) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S")
      
      #### checks ###
      # near equator, northeatward inter-cell distance in meters is
      icd <- distGeo(c(0,0), c(res(wr)[1],0))
      # and transit speed is 5 m/s divided by 8 neighbors
      vel <- 5^p / 8
      # check that they're similar
      plot(1/wr$E^p, main=icd/vel)
      # oh yeah
      
      
      saveRDS(wr, "data/windrose/isotropy_test_data.rds")
}
wr <- readRDS("data/windrose/isotropy_test_data.rds")


#####



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


ws <- function(x, y, windrose, 
               radius = 1000){
      
      coords <- c(x, y)
      origin <- matrix(coords, ncol=2)
      message(paste(coords, collapse=" "))
      
      # constrain analysis to region around focal point (for computatioal speed)
      coords_ll <- SpatialPoints(origin, crs(windrose)) %>%
            spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
            coordinates()
      circle <- geo_circle(coords_ll, width = radius * 1000)
      
      wshd <- windrose %>% crop(circle) %>% mask(circle) %>% 
            add_coords() %>% 
            transition_stack(windflow, directions=8, symm=F, direction="downwind") %>%
            accCost(coords)
      
      dst <- wshd
      dst[] <- pointDistance(coords, coordinates(dst), lonlat=T)
      stack(wshd, dst)
}


lats <- seq(1, 85, 1)
lats <- 5

d <- lats %>%
      lapply(function(x) ws(85, x, wr, 500)) %>%
      lapply(rasterToPoints) %>%
      lapply(as.data.frame) %>%
      lapply(function(x) mutate(x, id=y[1])) %>%
      bind_rows() %>%
      rename(wind = layer.1, dist = layer.2) %>%
      filter(is.finite(wind))


p <- ggplot(d, aes(x, y, fill=wind)) +
      facet_wrap(~id, scales="free") +
      geom_raster() +
      scale_fill_gradientn(colors=c("yellow", "red", "blue", "black")) +
      theme_void() +
      theme(strip.text=element_blank(),
            legend.position="none")
ggsave("figures/windsheds/test_windsheds.png", p, width=8, height=7, units="in")

p <- ggplot(d, aes(dist, wind, color=factor(id))) +
      geom_smooth(se=F)

betas <- d %>% filter(is.finite(wind)) %>%
      split(.$id) %>%
      map(~ lm(wind ~ dist, data = .)) %>%
      map(coef) %>% map_dbl("dist") %>%
      data.frame(beta = .,
                 lat = lats)
p <- ggplot(betas, aes(lat, beta)) + geom_point() + ylim(0, NA) +
      theme_minimal() +
      labs(x = "latitude",
           y = "beta: wind ~ distance")
ggsave("figures/windsheds/test_betas.png", p, width=6, height=6, units="in")




