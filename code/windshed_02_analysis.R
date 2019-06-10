

# to do: 
#     get gdistance wrapping on a sphere
#     downweight oceanic conductance


library(windshed)
library(tidyverse)
library(raster)
library(gdistance)
library(colormap)
library(grid)
library(gridExtra)

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
      transform <- function(x) .05 ^ (x * 1e9) # .05
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
rose <- stack("data/roses_force/cfsr_climatology/windrose_cfsr.tif") %>%
      rotate() %>%
      to_equal_area()

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)

# mean temperature (kelvin)
climate <- list.files("data/cfsr_monthly/tmp2m", full.names=T)[1:24] %>%
      stack() %>%
      mean() %>%
      rotate() %>%
      to_equal_area() %>%
      mask(land)
climate <- stack(climate, climate + 3)

# test
coords <- c(-5e6, -2e6)
woc(coords[1], coords[2], rose, climate)



pts <- land %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(x < 1.5e7,
             y > -7.75e6)
d <- map2(pts$x, pts$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, radius=500) %>%
      do.call("rbind", .) %>%
      as.data.frame()
write_csv(d, "data/windshed/force_500km.csv")




#####################


d <- read_csv("data/windshed/force_500km.csv")

dl <- d %>% gather(stat, value, wind_fwd:overlap_rev) %>%
      separate(stat, c("stat", "direction")) %>%
      mutate(stat = factor(stat, levels=c("wind", "clim", "overlap")),
             direction = factor(direction, levels=c("fwd", "rev"),
                                labels=c("forward", "reverse"))) %>%
      select(-runtime) %>%
      spread(stat, value) %>%
      mutate(windfill = overlap / clim)


#### plots of basic stats ####

for(var in c("wind", "clim", "overlap")){
      map <- ggplot(dl, aes_string("x", "y", fill=var)) +
            facet_grid(direction ~ .) +
            geom_raster() +
            scale_fill_viridis_c(trans="log10") +
            theme_void() +
            theme(text=element_text(size=20),
                  strip.text=element_text(size=20, angle=-90),
                  legend.position=c(.1, .7)) +
            guides(fill=guide_colorbar(barheight=15))
      ggsave(paste0("figures/windsheds/", var, ".png"), 
             width=12, height=12, units="in")
}



#### suitable area vs windfilling of suitable area ####

rnk <- function(x) ecdf(x)(x)
dd <- dl %>%
      filter(is.finite(windfill)) %>%
      mutate(color = colors2d(cbind(rnk(.$clim), rnk(.$windfill)),
                              c("forestgreen", "yellow", "red", "black")))

map <- ggplot(dd, aes(x, y)) +
      facet_grid(direction ~ .) +
      geom_raster(fill=dd$color) +
      theme_void() +
      theme(strip.text=element_text(size=50, angle=-90))
legend <- ggplot(dd, aes(clim, windfill)) +
      facet_grid(direction ~ .) +
      geom_point(color=dd$color, size=.2) +
      scale_y_log10(breaks=c(.001, .003, .01, .03, .1, .3, 1)) +
      xlim(0, .75) +
      theme_minimal() +
      theme(text=element_text(size = 45),
            strip.text=element_blank()) +
      labs(x = "suitable area",
           y = "proportion windfilling")
p <- arrangeGrob(legend, map, ncol=2, widths=c(1, 2))
png("figures/windsheds/windfill.png", width=3000, height=2000)
grid.draw(p)
dev.off()







#### fwd versus rev ####

for(var in c("wind", "clim", "overlap", "windfill")){
      dd <- dl %>%
            filter(is.finite(windfill)) %>%
            gather(stat, value, wind:windfill) %>%
            filter(stat == var) %>%
            spread(direction, value) %>%
            mutate(color = colors2d(cbind(rnk(.$forward), rnk(.$reverse)),
                                    c("black", "cyan", "gray85", "magenta")))
      map <- ggplot(dd, aes_string("x", "y")) +
            geom_raster(fill = dd$color) +
            theme_void() +
            theme(text=element_text(size=20),
                  strip.text=element_text(size=20, angle=-90),
                  legend.position=c(.1, .7)) +
            guides(fill=guide_colorbar(barheight=15))
      
      trans <- switch(var,
                      overlap="log10",
                      wind="log10",
                      windfill="log10",
                      clim="identity")
      
      legend <- ggplot(dd, aes(forward, reverse)) +
            geom_point(color=dd$color, size=.1) +
            scale_y_continuous(trans=trans) +
            scale_x_continuous(trans=trans) +
            theme_minimal() +
            theme(text=element_text(size = 25)) +
            labs(x = paste("forward", var),
                 y = paste("reverse", var))
      png(paste0("figures/windsheds/", var, "_FR.png"), 
          width=3000, height=2000)
      plot(map)
      plot(legend, vp=viewport(x=.15, y=.33, width=.2, height=.4))
      dev.off()
}


#####################

ggplot(d, aes(overlap_fwd, overlap_rev)) +
      geom_point()

ggplot(d, aes(overlap_fwd, overlap_rev)) +
      geom_point() +
      scale_x_log10() + scale_y_log10()

ggplot(d, aes(wind_fwd, overlap_fwd, color=clim_fwd)) +
      geom_point() +
      scale_color_viridis_c(trans="log10") +
      scale_x_log10() + scale_y_log10()


ggplot(d %>% filter(x < 1.5e7), aes(x, y, color=sqrt(runtime))) +
      geom_point() +
      scale_color_viridis_c(trans="log10")

d %>% filter(x < 1.5e7,
             y > -7.75e6) %>% 
      ggplot(aes(y, runtime)) + geom_point()




# get continentatlity and elevation covariates
# get biome covariates


