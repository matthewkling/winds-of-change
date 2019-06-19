
# visualize stats across a (mountain range) transect

# elevation
# temperature
# wind speed/direction

# local, forward, reverse velocity
# local, forward, reverse windshed
# local alignment
# forward, reverse windfilling


library(windscape)
library(raster)
library(tidyverse)
library(geosphere)
library(grid)
library(gridExtra)


land <- raster("f:/cfsr/land.tif") %>% 
      rotate()

# load windrose data
rose <- stack("data/roses_velocity/cfsr_climatology/roses_cfsr_1980s.tif") %>%
      rotate()

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)
names(rose) <- windrose_names()

# mean temperature (kelvin)
climate <- list.files("data/cfsr_monthly/tmp2m", full.names=T)[1:24] %>%
      stack() %>%
      mean() %>%
      rotate() %>%
      mask(land)
climate <- stack(climate, climate + 3)
names(climate) <- c("clim0", "clim1")




transect_plot <- function(transect = NULL, n=25){
      
      if(is.null(transect)){
            message("Click the map twice to draw an EAST-WEST transect, then hit ESC.")
            plot(climate[[1]], 
                 col=colorRampPalette(c("gray98" "cyan", "dodgerblue", "darkorchid", 
                                        "red", "yellow"))(100))
            transect <- drawLine()
      }
      
      coords <- coordinates(transect)[[1]][[1]]
      tpoints <- apply(coords, 2, function(x) seq(x[1], x[2], length.out=n))
      dst <- distGeo(tpoints[1,], tpoints) / 1000
      
      
      ros <- rose %>% 
            subset(c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) %>%
            stack(., .) %>%
            raster::extract(tpoints, method="bilinear")
      clim <- raster::extract(climate, tpoints, method="bilinear")
      
      
      bearing <- bearingRhumb(coords[1,], coords[2,])
      octant <- floor(bearing / 45) + 1
      remainder <- bearing %% 45
      
      fx <- function(x) sum(x * c(remainder, 45-remainder))/45
      tw <- ros[, c(octant, octant+1)] %>% apply(1, fx)
      hw <- ros[, c(octant+4, octant+5)] %>% apply(1, fx)
      w <- tw - hw
      
      d <- cbind(dst, clim, w) %>%
            as.data.frame() %>%
            filter(!is.na(clim0)) %>%
            mutate(uclim = 1/(lag(clim0) - lead(clim0)),
                   uwind = w)
      
      
      mx <- diff(range(na.omit(d$dst))) / 5
      mn <- mx / 5
      norm <- function(x){
            x <- x / max(abs(x), na.rm=T) * mx
            inc <- is.finite(x) & abs(x) < mn
            x[inc] <- mn * sign(x[inc])
            x
      }
      
      d <- d %>% 
            arrange(dst) %>%
            mutate(uclim = sign(uclim) * scales::rescale(log(abs(uclim)), c(.01, 1)),
                  uwind = norm(uwind),
                  uclim = norm(uclim),
                  uwind = ifelse(is.na(uclim), NA, uwind))
      
      arw <- arrow(type="closed", angle=15, length = unit(0.1, "inches"))
      ynudge <- 0
      
      wind <- ggplot(d) +
            geom_ribbon(aes(dst, ymin=clim0, ymax=max(clim0 + 1, na.rm=T)), 
                        stat="identity", fill="gray90") +
            geom_segment(aes(x=dst, xend=dst+uwind,
                             y=clim0+ynudge, yend=clim0+ynudge),
                         arrow=arw, color="dodgerblue") +
            geom_segment(aes(x=dst, xend=dst+uclim,
                             y=clim0-ynudge, yend=clim0-ynudge),
                         arrow=arw, color="black") +
            geom_point(aes(dst, clim0, color=sign(uwind)==sign(uclim)), size=4) +
            scale_color_manual(values=c("darkred", "forestgreen")) +
            annotate(geom="text", 
                     x=range(d$dst, na.rm=T) + diff(range(d$dst, na.rm=T)) / 20 * c(1,-1), 
                     y=max(d$clim0 + .5), hjust=c(0, 1),
                     label=c("climate tracking", "prevailing wind"), 
                     color=c("black", "dodgerblue")) +
            theme_minimal() +
            theme(legend.position="none") +
            scale_y_reverse() +
            labs(y = "temperature (deg C)", 
                 x = "distance along transect (km)")
      print(wind)
      
      return(list(transect = transect,
           data = d,
           wind = wind))
      
      ####################
      # deprecated analog plot code below -- fix
      
      tol <- 2
      for(i in 1:nrow(d)){
            
            # forward
            #analog <- d$clim1 > (d$clim0[i]-tol) & d$clim1 < (d$clim0[i]+tol)
            analog <- d$clim1 < (d$clim0[i]+tol)
            distance <- abs(d$x - d$x[i])
            distance[!analog] <- NA
            a <- which.min(distance)
            xf <- d$x[a]
            cf <- d$clim1[a]
            if(length(a)==0) xf <- cf <- NA
            d$xf[i] <- xf
            d$cf[i] <- cf
            
            # reverse
            #analog <- d$clim0 > (d$clim1[i]-tol) & d$clim0 < (d$clim1[i]+tol)
            analog <- d$clim0 > (d$clim1[i]-tol)
            distance <- abs(d$x - d$x[i])
            distance[!analog] <- NA
            a <- which.min(distance)
            xr <- d$x[a]
            cr <- d$clim0[a]
            if(length(a)==0) xr <- cr <- NA
            d$xr[i] <- xr
            d$cr[i] <- cr
      }
      
      forward <- ggplot(d) +
            geom_line(aes(x, clim0)) +
            geom_line(aes(x, clim1)) +
            geom_point(aes(x, clim0), color="red") +
            geom_point(aes(x, clim1)) +
            geom_segment(aes(x, clim0, xend=xf, yend=cf), 
                         color="red", arrow=arw) +
            scale_y_reverse()
      
      reverse <- ggplot(d) +
            geom_line(aes(x, clim0)) +
            geom_line(aes(x, clim1)) +
            geom_point(aes(x, clim0)) +
            geom_point(aes(x, clim1), color="dodgerblue") +
            geom_segment(aes(xr, cr, xend=x, yend=clim1), 
                         color="dodgerblue", arrow=arw) +
            scale_y_reverse()
      
      analogs <- arrangeGrob(forward, reverse, ncol=1)
      
      list(transect = transect,
           data = d,
           wind = wind, 
           analogs = analogs)
}

transect_plot()



add_tag <- function(x, tag){
      x$wind <- x$wind +
            annotate(geom="text", label=tag,
                     x=min(x$data$x, na.rm=T),
                     y=min(x$data$clim0, na.rm=T),
                     hjust=0, size=6)
      x
}

# good ones
brz <- transect_plot(brz$transect) %>% add_tag("\nBrazilian\ncerrado")
aus <- transect_plot(aus$transect) %>% add_tag("\nAustralian\ncontinent")
app <- transect_plot(app$transect) %>% add_tag("\nSouthern\nAppalachia")
bc <- transect_plot(bc$transect) %>% add_tag("\nBritish\nColumbia")

p <- arrangeGrob(bc$wind, brz$wind, app$wind, aus$wind, ncol=2)
ggsave("figures/transects/transects.png", p, width=12, height=9, units="in")
