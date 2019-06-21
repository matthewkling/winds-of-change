

library(raster)
library(tidyverse)
library(rgdal)
library(colormap)
library(grid)
library(gridExtra)
library(fasterize)
library(sf)



px <- read_csv("data/windshed/force_500km_v2.csv")

land <- raster("f:/cfsr/land.tif") %>% rotate()
latlon <- crs(land)
hd <- CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
land <- land %>% projectRaster(crs=hd, method="ngb")

world <- map_data("world")
coordinates(world) <- c("long", "lat")
crs(world) <- latlon
world <- spTransform(world, hd)
world <- as.data.frame(world)  

maps <- function(sp){
      #sp = "Sequoia sempervirens"
      range <- readOGR(paste0("F:/little_trees/raw_data/", sp))
      crs(range) <- latlon
      range <- spTransform(range, hd)
      template <- crop(land, range)
      rangel <- range %>% as("SpatialLines") %>% rasterize(template)
      range <- range %>% st_as_sf() %>% fasterize(template)
      range <- stack(range, rangel) %>% 
            sum(na.rm=T) %>% mask(template) %>%
            reclassify(c(-1, 0.01, NA)) %>%
            rasterToPoints() %>% as.data.frame()
      
      f <- left_join(range, px) %>%
            select(-layer, -runtime) %>%
            gather(var, value, -x, -y) %>%
            separate(var, c("property", "direction", "moment", "stat"), sep="_")
            
      d <- f %>%
            filter(moment == "windshed",
                   stat == "size") %>%
            spread(property, value) %>%
            mutate(windfill = overlap / clim) %>%
            gather(property, value, clim, wind, overlap, windfill) %>%
            mutate(direction = ifelse(direction == "fwd", "forward", "reverse"))
      
      dd <- d %>%
            select(-moment, -stat) %>%
            spread(property, value) %>%
            mutate(color = colors2d(cbind(rank(.$clim), rank(.$windfill)),
                                    c("forestgreen", "yellow", "red", "black")))
      
      map <- ggplot(dd, aes(x, y)) +
            facet_grid(. ~ direction) +
            geom_polygon(data=world, aes(long, lat, group=group), fill="gray80", color="white") +
            geom_raster(fill=dd$color) +
            theme_void() +
            theme(strip.text=element_text(size=15)) +
            coord_cartesian(xlim = range(dd$x), ylim = range(dd$y))
      legend <- ggplot(dd, aes(clim, windfill)) +
            #facet_grid(direction ~ .) +
            geom_point(color=dd$color, size=2) +
            scale_y_log10(breaks=c(.001, .003, .01, .03, .1, .3, 1)) +
            #xlim(0, .75) +
            theme_minimal() +
            theme(strip.text=element_blank()) +
            labs(x = "analogous area",
                 y = "proportion windfilling")
      p <- arrangeGrob(legend, map, nrow=1, widths=c(1, 2))
      ggsave(paste0("figures/species/windfill_clim/", sp, ".png"), p, width=9, height=3, units="in")
}

species <- c("Quercus lobata", "Quercus douglasii", "Sequoia sempervirens", 
             "Pinus ponderosa", "Abies balsamea", "Acer rubrum", 
             "Populus balsamifera", "Populus tremuloides",
             "Pseudotsuga menziesii")

lapply(species, maps)

