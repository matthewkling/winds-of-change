

library(raster)
library(tidyverse)
library(ineq)
library(grid)
library(rNOMADS)
library(colormap)

# load data

to_equal_area <- function(r, ...){
      projectRaster(r, crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), ...)}
to_equal_area <- function(r, ...){
      projectRaster(r, crs=CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), ...)}


land <- raster("f:/cfsr/land.tif") %>% 
      rotate()# %>% 
      #to_equal_area(method="ngb")

windrose <- stack("data/roses_force/cfsr_climatology/roses_cfsr_1980s.tif") %>%
      rotate() %>%
      #to_equal_area() %>%
      mask(land)

# summarize

total <- sum(windrose)
directionality <- calc(windrose, Gini) # Gini: 1 = directionality, 0 = equality


s <- stack(total, directionality) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      rename(total = layer.1,
             directionality = layer.2) %>%
      filter(is.finite(total), is.finite(directionality))



### prevailing direction data
agg <- 10
modir <- "data/cfsr_monthly"
mf <- list.files(modir, full.names=T, recursive=T)
u <- mf[grepl("_u", mf)] %>% stack() %>% mean() %>% rotate() %>% #mask(land) %>% 
      aggregate(agg)
v <- mf[grepl("_v", mf)] %>% stack() %>% mean() %>% rotate() %>% #mask(land) %>% 
      aggregate(agg)
uv <- stack(u, v)
uvv <- values(uv)
nna <- !is.na(uvv[,1])
w <- MagnitudeAzimuth(uvv[nna,1], uvv[nna,2])
wm <- uv
wm[[1]][nna] <- w$magnitude
wm[[2]][nna] <- w$azimuth

a <- wm %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      rename(a = layer.2, m = layer.1) %>%
      mutate(radians = (90 - a) / 190 * pi)
aland <- wm %>% mask(aggregate(land, agg)) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      rename(a = layer.2, m = layer.1) %>%
      mutate(radians = (90 - a) / 190 * pi)


## plot

# s$color <- colors2d(dplyr::select(s, total, directionality),
#                     c("gold", "forestgreen", "darkblue", "red"),
#                     xtrans="rank", ytrans="rank")

s$color <- colors2d(dplyr::select(s, total, directionality),
                    c("yellow", "green", "dodgerblue", "red"),
                    xtrans="rank", ytrans="rank")

ocean <- land %>%
      reclassify(c(NA, NA, 1,
                   .5, 1.5, NA)) %>%
      rasterToPoints() %>%
      as.data.frame()

map <- ggplot() +
      geom_raster(data=s, aes(x, y), fill = s$color) +
      geom_spoke(data=aland, aes(x, y, angle=radians, radius=.1),
                 arrow = arrow(type="open", angle=10, length=unit(.2, "in")),
                 color = "black", size = .5) +
      #geom_raster(data=ocean, aes(x, y), fill = "white") +
      theme_void() +
      coord_cartesian(xlim=c(-180, 180),
                      ylim=c(-90, 90),
                      expand = 0)

legend <- ggplot(s, aes(total, directionality)) +
      geom_point(color=s$color, size=.1) +
      theme_minimal() +
      scale_x_log10() +
      theme(text=element_text(size = 25)) +
      labs(x = "mean wind force")

png("figures/tailwinds/speed_isotropy_direction.png", width=2000, height=1000)
plot(map)
plot(legend, vp=viewport(x=.12, y=.35, width=.22, height=.44))
dev.off()





climate <- list.files("data/cfsr_monthly/tmp2m", full.names=T)[1:24] %>%
      stack() %>%
      mean() %>%
      rotate() %>%
      mask(land) %>%
      rasterToPoints() %>%
      as.data.frame()

map <- ggplot() +
      geom_raster(data=climate, aes(x, y, fill = layer)) +
      geom_spoke(data=a, aes(x, y, angle=radians, radius=.1),
                 arrow = arrow(type="open", angle=10, length=unit(.2, "in")),
                 color = "black", size = .5) +
      scale_fill_gradientn(colors=c("gray97", "cyan", "dodgerblue", "darkblue",
                                    "darkorchid", "darkred", "red", "yellow")) +
      theme_void() +
      theme(legend.position="none") +
      coord_cartesian(xlim=c(-180, 180),
                      ylim=c(-90, 90),
                      expand=0)

png("figures/tailwinds/temperature_direction.png", width=2000, height=1000)
plot(map)
dev.off()
