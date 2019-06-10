

library(raster)
library(tidyverse)
library(ineq)
library(grid)

# load data

to_equal_area <- function(r, ...){
      projectRaster(r, crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), ...)}

land <- raster("f:/cfsr/land.tif") %>% 
      rotate() %>% 
      to_equal_area(method="ngb")

windrose <- stack("data/roses_force/cfsr_climatology/windrose_cfsr.tif") %>%
      rotate() %>%
      to_equal_area() %>%
      mask(land)


# summarize

total <- sum(windrose) %>% log10()
directionality <- calc(windrose, Gini) # Gini: 1 = directionality, 0 = equality


# plot

s <- stack(total, directionality) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      rename(total = layer.1,
             directionality = layer.2)
s$color <- colors2d(dplyr::select(s, total, directionality),
                    c("yellow", 
                      "forestgreen", 
                      "darkblue",
                      "red"),
                    xtrans="rank", ytrans="rank")

map <- ggplot(s, aes(x, y)) +
      geom_raster(fill = s$color) +
      theme_void()

legend <- ggplot(s, aes(total, directionality)) +
      geom_point(color=s$color, size=.1) +
      theme_minimal() +
      theme(text=element_text(size = 25),
            axis.text.x=element_blank()) +
      labs(x = "log wind energy")

png("figures/tailwinds/speed_isotropy.png", width=2000, height=1000)
plot(map)
plot(legend, vp=viewport(x=.15, y=.33, width=.2, height=.4))
dev.off()
