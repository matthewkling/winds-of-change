

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


# higher definition version


### prevailing direction data
agg <- 3
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

map <- ggplot() +
      geom_raster(data=climate, aes(x, y, fill = layer)) +
      geom_spoke(data=a, aes(x, y, angle=radians, radius=.01),
                 arrow = arrow(type="open", angle=10, length=unit(.1, "in")),
                 color = "black", size = .25) +
      scale_fill_gradientn(colors=c("gray97", "cyan", "dodgerblue", "darkblue",
                                    "darkorchid", "darkred", "red", "yellow")) +
      theme_void() +
      theme(legend.position="none") +
      coord_cartesian(xlim=c(-180, 180),
                      ylim=c(-90, 90),
                      expand=0)

png("figures/tailwinds/temperature_direction_hd.png", width=3000, height=1500)
plot(map)
dev.off()



### latitudinal patterns in meridional winds ###

v <- mf[grepl("_v", mf)] %>% stack() %>% mean() %>% rotate() %>% stack(land) %>%
      rasterToPoints() %>% as.data.frame()

v <- v %>%
      mutate(latitude = abs(y),
             polarity = ifelse(y>0, layer, -layer),
             land = ifelse(is.na(land), "water", "land"))

p <- ggplot(v, aes(latitude, polarity, color=land)) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=c(0, 30, 60, 90)) +
      geom_smooth(se=F) +
      annotate(geom="text", x=c(15, 45, 75), y=1,
              label=c("Hadley cell", "Ferrel cell", "Polar cell")) +
      scale_x_continuous(limits=c(0,90), breaks=c(0, 30, 60, 90), expand=c(0,0)) +
      theme_minimal() +
      theme(legend.position=c(.5, .1)) +
      labs(x = "degrees poleward",
           y = "poleward wind speed (m/s)",
           color = NULL)
ggsave("figures/tailwinds/meridional.png", p, width=8, height=6, units="in")


vs <- v %>%
      bind_rows(mutate(v, land="all")) %>%
      mutate(lat = plyr::round_any(latitude, 1)) %>%
      group_by(land, lat) %>%
      summarize(v = mean(polarity),
                v3 = mean(polarity^3))

p <- ggplot(vs, aes(lat, ymin=0, ymax=v, y=v, color=land, fill=land)) +
      geom_ribbon(alpha=.1, position="identity") +
      geom_line() +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=c(0, 30, 60, 90)) +
      annotate(geom="text", x=c(15, 45, 75), y=1.05,
               label=c("Hadley cell", "Ferrel cell", "Polar cell")) +
      scale_color_manual(values=c("red", "forestgreen", "blue")) +
      scale_fill_manual(values=c("red", "forestgreen", "blue")) +
      scale_x_continuous(limits=c(0,90), breaks=c(0, 30, 60, 90), expand=c(0,0)) +
      theme_minimal() +
      theme(legend.position=c(.5, .1)) +
      labs(x = "degrees poleward",
           y = "poleward wind speed (m/s)",
           color=NULL, fill=NULL)
ggsave("figures/tailwinds/meridional_bins.png", p, width=8, height=6, units="in")





vs <- v %>%
      bind_rows(mutate(v, land="all")) %>%
      mutate(lat = plyr::round_any(y, 5)) %>%
      group_by(land, lat) %>%
      summarize(v = mean(polarity),
                v3 = mean(polarity^3))

ymn <- -5

p <- ggplot(vs, aes(lat, ymin=0, ymax=v, y=v, color=land, fill=land)) +
      geom_ribbon(alpha=.1, position="identity") +
      geom_line() +
      geom_hline(yintercept=0) +
      annotate(geom="segment", x=xb, xend=xb, y=0, yend=1) +
      annotate(geom="text", x=c(-75, -45, -15, 15, 45, 75), y=1,
               angle=c(-75, -45, -15, 15, 45, 75) + c(90, 90, 90, -90, -90, -90), 
               lineheight=.7, hjust=0.5,
               label=c("Polar\ncell", "Ferrel\ncell", "Hadley\ncell", 
                       "Hadley\ncell", "Ferrel\ncell", "Polar\ncell")) +
      annotate(geom="text", x=xb, y=-.25, angle=xb, hjust=1, size=3,
               label=c("Polar high", "Subpolar low", "Subtropical high", "Equatorial low", 
                       "Subtropical high", "Subpolar low", "Polar high")) +
      scale_x_continuous(limits=c(-90, 270), breaks=xb, expand=c(0,0)) +
      scale_color_manual(values=c("red", "forestgreen", "blue")) +
      scale_fill_manual(values=c("red", "forestgreen", "blue")) +
      theme_minimal() +
      theme(legend.position=c(.55, .52),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank()) +
      labs(x = NULL,
           y = "poleward windspeed (m3/s)",
           color = NULL, fill = NULL) +
      coord_polar(start=pi, direction=-1) +
      ylim(ymn, NA)
ggsave("figures/tailwinds/meridional_bins_polar.png", p, width=8, height=8, units="in")

