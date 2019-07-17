


library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(colormap)
library(grid)
library(gridExtra)
library(ecoclim)

select <- dplyr::select


unwrap <- function(x, width=20){
      x <- x %>% crop(extent(180-width, 180, -90, 90)) %>% shift(-360) %>% merge(x)
      x <- x %>% crop(extent(-180, -180+width, -90, 90)) %>% shift(360) %>% merge(x)
      x
}


# calculate the spatial distribution of climate analogs for a given location
analogs <- function(climate, # raster stack of first, second timesteps
                    coords, # lon-lat vector
                    reverse = FALSE, # should reverse analogs be calculated, instead of forward
                    sigma = 1 # standard deviation of similarity kernel, in degrees Celsius
){
      
      if(reverse) climate <- climate[[2:1]]
      target <- raster::extract(climate[[1]], coords)
      kernel <- function(x) exp(-.5*(x/sigma)^2)
      calc(climate[[2]] - target, kernel)
}




land <- raster("f:/cfsr/land.tif") %>% 
      rotate() %>% unwrap(180)

# load windrose data
rose <- stack("data/windrose/windrose_p1_2000s.tif") %>%
      rotate() %>% unwrap(180)

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)

# load current/future climate data
climate <- stack("data/geographic/processed/temperature.tif") %>% unwrap(180)

sigma <- 2

maps <- function(x){
      
      location <- x$location
      coords <- x$coords
      ext <- x$ext
      
      clim <- crop(climate, ext)
      wind <- rose %>% crop(ext) %>% add_coords() 
      downtrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="downwind")
      uptrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="upwind")
      
      #px <- land %>% crop(ext) %>% rasterToPoints() %>% as.data.frame()
      #pts <- sample_n(px, 4)
      
      pts <- data.frame(x=coords[1], y=coords[2])
      
      df <- map2_df(pts$x, pts$y, function(xx, yy){
            origin <- c(xx, yy) %>% matrix(ncol=2)
            stack(accCost(downtrans, origin) / 3600,
                  accCost(uptrans, origin) / 3600,
                  analogs(clim, origin, sigma=sigma),
                  analogs(clim, origin, sigma=sigma, reverse = T)) %>%
                  mask(crop(land, ext)) %>%
                  rasterToPoints() %>%
                  as.data.frame() %>%
                  rename(downwind = layer.1,
                         upwind = layer.2,
                         forward = layer.3,
                         reverse = layer.4) %>%
                  mutate(px = xx, py = yy)
      })
      
      d <- df  %>%
            select(-forward, -reverse) %>%
            gather(direction, hours, downwind, upwind)
      
      dd <- d %>% group_by(px, py, direction) %>% sample_n(1)
      
      m <- 2
      
      p <- ggplot(d,  aes(x, y, z=hours, fill=hours)) +
            facet_grid(px ~ direction) +
            geom_raster() +
            geom_contour(binwidth=100, color="black", alpha=.5, size=1) +
            geom_point(data=dd, aes(px, py), 
                       size=3, color="white") +
            annotate(geom="segment", 
                     x=dd$px+c(0, m, 0, -m), xend=dd$px+c(0, m*3, 0, -m*3),
                     y=dd$py+c(m, 0, -m, 0), yend=dd$py+c( m*3, 0, -m*3, 0), 
                     color="white", size=1.5) +
            geom_text(data=dd, aes(mean(range(d$x)), max(d$y), 
                                   label=paste0("\n", direction)), 
                      size=15, color="white", hjust=.5) +
            scale_fill_gradientn(colors=c("cyan", "purple", "darkred", "#290000")) +
            scale_y_continuous(expand=c(0,0), limits=c(min(d$y - 5), NA)) +
            theme_void() +
            theme(strip.text=element_blank(),
                  text = element_text(size=35, color="white"),
                  plot.background=element_rect(fill="black"),
                  legend.position="bottom") +
            guides(fill=guide_colorbar(barwidth=35)) +
            labs(fill = "wind-hours from origin\n\n")
      ggsave(paste0("figures/windsheds/examples/examples_", location, "_windtime.png"), p,
             width=16, height=9, units="in")
      
      
      ypg <- 1 # years / generation
      spg <- 3600 # seconds aloft / generation
      spy <- spg/ypg # wind-seconds traveled / species-year
      
      p <- ggplot(d,  aes(x, y, z=hours, fill=hours*3600/spy)) +
            facet_grid(px ~ direction) +
            geom_raster() +
            geom_contour(binwidth=100, color="black", alpha=.5, size=1) +
            geom_point(data=dd, aes(px, py), 
                       size=3, color="white") +
            annotate(geom="segment", 
                     x=dd$px+c(0, m, 0, -m), xend=dd$px+c(0, m*3, 0, -m*3),
                     y=dd$py+c(m, 0, -m, 0), yend=dd$py+c( m*3, 0, -m*3, 0), 
                     color="white", size=1.5) +
            geom_text(data=dd, aes(mean(range(d$x)), max(d$y), 
                                   label=paste0("\n", direction)), 
                      size=15, color="white", hjust=.5) +
            scale_fill_gradientn(colors=c("cyan", "purple", "darkred", "#290000")) +
            scale_y_continuous(expand=c(0,0), limits=c(min(d$y - 5), NA)) +
            theme_void() +
            theme(strip.text=element_blank(),
                  text = element_text(size=35, color="white"),
                  plot.background=element_rect(fill="black"),
                  legend.position="bottom") +
            guides(fill=guide_colorbar(barwidth=35)) +
            labs(fill = "years to colonization\n(@ 1 hr aloft / 1-yr generation)\n")
      ggsave(paste0("figures/windsheds/examples/examples_", location, "_sptime.png"), p,
             width=16, height=9, units="in")
      
      
      d <- df  %>%
            select(-downwind, -upwind) %>%
            gather(direction, analog, forward, reverse) %>%
            filter(is.finite(analog))
      
      dd <- d %>% group_by(px, py, direction) %>% sample_n(1)
      
      p <- ggplot(d,  aes(x, y, z=analog, fill=analog)) +
            facet_grid(px ~ direction) +
            geom_raster() +
            geom_contour(color="black", alpha=.5, size=1) +
            geom_point(data=dd, aes(px, py), 
                       size=3, color="white") +
            annotate(geom="segment", 
                     x=dd$px+c(0, m, 0, -m), xend=dd$px+c(0, m*3, 0, -m*3),
                     y=dd$py+c(m, 0, -m, 0), yend=dd$py+c( m*3, 0, -m*3, 0), 
                     color="white", size=1.5) +
            geom_text(data=dd, aes(mean(range(d$x)), max(d$y), 
                                   label=paste0("\n", direction)), 
                      size=15, color="white", hjust=.5) +
            scale_fill_gradientn(colors=c("cyan", "purple", "darkred") %>% rev()) +
            scale_y_continuous(expand=c(0,0), limits=c(min(d$y - 5), NA)) +
            theme_void() +
            theme(strip.text=element_blank(),
                  text = element_text(size=35, color="white"),
                  plot.background=element_rect(fill="black"),
                  legend.position="bottom") +
            guides(fill=guide_colorbar(barwidth=35)) +
            labs(fill = "temperature similarity\n\n")
      ggsave(paste0("figures/windsheds/examples/examples_", location, "_analogs.png"), p,
             width=16, height=9, units="in")
      
      
      
      mx <- max(c(df$downwind, df$upwind))
      d <- df  %>%
            mutate(emigration = forward * (1-downwind/mx),
                   immigration = reverse * (1-upwind/mx)) %>%
            select(-forward, -reverse, -upwind, -downwind) %>%
            gather(direction, hours, immigration, emigration)
      
      dd <- d %>% group_by(px, py, direction) %>% sample_n(1)
      
      p <- ggplot(d,  aes(x, y, z=hours, fill=hours)) +
            facet_grid(px ~ direction) +
            geom_raster() +
            geom_contour(color="black", alpha=.5, size=1) +
            geom_point(data=dd, aes(px, py), 
                       size=3, color="white") +
            annotate(geom="segment", 
                     x=dd$px+c(0, m, 0, -m), xend=dd$px+c(0, m*3, 0, -m*3),
                     y=dd$py+c(m, 0, -m, 0), yend=dd$py+c( m*3, 0, -m*3, 0), 
                     color="white", size=1.5) +
            geom_text(data=dd, aes(mean(range(d$x)), max(d$y), 
                                   label=paste0("\n", direction)), 
                      size=15, color="white", hjust=.5) +
            scale_fill_gradientn(colors=c("cyan", "purple", "darkred") %>% rev()) +
            scale_y_continuous(expand=c(0,0), limits=c(min(d$y - 5), NA)) +
            theme_void() +
            theme(strip.text=element_blank(),
                  text = element_text(size=35, color="white"),
                  plot.background=element_rect(fill="black"),
                  legend.position="bottom") +
            guides(fill=guide_colorbar(barwidth=35)) +
            labs(fill = "wind-temperature overlap\n\n")
      
      ggsave(paste0("figures/windsheds/examples/examples_", location, "_overlap.png"), p,
             width=16, height=9, units="in")
}


maps(list(location = "skukuza", 
          coords = c(31.6, -25),
          ext = extent(0, 41, -42, -0)))

maps(list(location = "kalahari", 
          coords = c(23.8, -22.4),
          ext = extent(0, 41, -42, -0)))

maps(list(location = "yellowstone", 
          coords = c(-110.6, 44.6),
          ext = extent(-130, -90, 25, 65)))

maps(list(location = "rockymtn", 
          coords = c(-105.7, 40.3),
          ext = extent(-130, -90, 25, 65)))

maps(list(location = "amazon", 
          coords = c(-63, -3), # c(-52.5, -2),
          ext = extent(-81, -34, -26, 13)))

maps(list(location = "blackforest", 
          coords = c(8.2, 48.2),
          ext = extent(-10, 30, 36, 66)))






# load windrose data
rose <- stack("data/windrose/windrose_p1_2000s.tif") %>%
      rotate() %>% unwrap(180)

ext = extent(climate)
wind <- rose %>% crop(ext) %>% add_coords() 
downtrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="downwind")



location = "krugerGlobal" 
coords = c(31.6, -25)

location = "yellowstoneGlobal"
coords = c(-110.6, 44.6)

location = "amazonGlobal"
coords = c(-68, -7)

location = "lhasaGlobal"
coords = c(91, 30)



downwind <- accCost(downtrans, matrix(coords, ncol=2))
d1 <- crop(downwind, extent(-360, 0, -90, 90)) %>% shift(360)
d2 <- crop(downwind, extent(0, 360, -90, 90))
downwind <- min(stack(d1, d2)) %>% rotate() %>%
      rasterToPoints() %>%
      as.data.frame()

dd <- data.frame(px=coords[1], py=coords[2])
m <- 5

p <- ggplot() +
      geom_raster(data=downwind, aes(x, y, fill=layer)) +
      geom_contour(data=downwind, aes(x, y, z=layer), binwidth=max(downwind$layer)/200,
                   color="black", alpha=.5, size=1) +
      geom_point(data=dd, aes(px, py), 
                 size=3, color="white") +
      annotate(geom="segment", 
               x=dd$px+c(0, m, 0, -m), xend=dd$px+c(0, m*3, 0, -m*3),
               y=dd$py+c(m, 0, -m, 0), yend=dd$py+c( m*3, 0, -m*3, 0), 
               color="white", size=1.5) +
      coord_fixed() +
      scale_fill_gradientn(colors=c("cyan", "purple", "darkred", "#290000")) +
      theme_void() +
      theme(legend.position = "none",
            plot.background = element_rect(fill='black', color="black"))


ggsave(paste0("figures/windsheds/examples/examples_", location, "_overlap.png"), p,
       width=32, height=18, units="in")
