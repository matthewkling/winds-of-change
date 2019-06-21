
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
            message("Click the map twice to draw a bounding box around the area of interest.")
            plot(climate[[1]], 
                 col=colorRampPalette(c("gray98", "cyan", "dodgerblue", "darkorchid", 
                                        "red", "yellow"))(100))
            ext <- drawExtent()
            plot(climate[[1]] %>% crop(ext), 
                 col=colorRampPalette(c("gray98", "cyan", "dodgerblue", "darkorchid", 
                                        "red", "yellow"))(100))
            message("Click the map twice to draw an EAST-WEST transect, then hit ESC.")
            transect <- drawLine()
      }
      
      coords <- coordinates(transect)[[1]][[1]]
      tpoints <- apply(coords, 2, function(x) seq(x[1], x[2], length.out=n))
      colnames(tpoints) <- c("x", "y")
      dst <- distGeo(tpoints[1,], tpoints) / 1000
      
      
      ros <- rose %>% 
            subset(c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) %>%
            stack(., .) %>%
            raster::extract(tpoints, method="bilinear")
      clim <- raster::extract(climate, tpoints, method="bilinear")
      
      
      bearing <- bearingRhumb(coords[1,], coords[2,])
      octant <- floor(bearing / 45) + 1
      remainder <- bearing %% 45
      
      # net wind in transect direction
      fx <- function(x) sum(x * c(remainder, 45-remainder))/45
      tw <- ros[, c(octant, octant+1)] %>% apply(1, fx)
      hw <- ros[, c(octant+4, octant+5)] %>% apply(1, fx)
      w <- tw - hw
      
      total <- ros[,1:8] %>% apply(1, sum)
      prop <- (tw + hw) / total
      message(mean(prop, na.rm=T))
      
      
      d <- cbind(dst, clim, w, tpoints) %>%
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
                   uwind = ifelse(is.na(uclim), NA, uwind),
                   aligned = sign(uwind)==sign(uclim),
                   aligned = factor(aligned, levels=c(T, F)))
      
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
            geom_point(aes(dst, clim0, color=aligned), size=4) +
            scale_color_manual(values=c("forestgreen", "darkred"), drop=F) +
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
                  prop = prop,
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

#tp <- transect_plot(tp$transect)


##################################

australia <- transect_plot(australia$transect)
greenland <- transect_plot(greenland$transect)
rockies <- transect_plot(rockies$transect)
cali <- transect_plot(cali$transect)
russia <- transect_plot(russia$transect)
sahara <- transect_plot(sahara$transect)

d <- list(rockies, cali, greenland, 
          australia, russia, sahara)

saveRDS(d, "data/transect_data.rds")
d <- readRDS("data/transect_data.rds")

labels <- c("Transverse mountain wind", "(California)", "Katabatic wind",
            "Cross-desert wind", "Southwesterly", "Thermal low")
for(i in 1:length(d)) d[[i]]$data$id <- labels[i]

d <- lapply(d, function(x) x$data) %>%
      bind_rows() %>%
      mutate(id=factor(id, levels=labels))

d <- d %>%
      group_by(id) %>%
      mutate(ymax = max(clim0) + .1 * diff(range(clim0, na.rm=T))) %>%
      mutate(i = as.integer(id),
             idi = paste0(i, ": ", id))

ynugde <- 0
arw <- arrow(type="closed", angle=15, length = unit(0.05, "inches"))
p <- ggplot(d) +
      facet_wrap(~idi, scales="free") +
      geom_ribbon(aes(dst, ymin=clim0, ymax=ymax), 
                  stat="identity", fill="gray90") +
      geom_point(aes(dst, clim0, color=aligned), size=2) +
      geom_segment(aes(x=dst, xend=dst+uwind,
                       y=clim0+ynudge, yend=clim0+ynudge, color=aligned),
                   arrow=arw, #color="blue", 
                   size=.3) +
      # geom_segment(aes(x=dst, xend=dst+uclim,
      #                  y=clim0-ynudge, yend=clim0-ynudge),
      #              arrow=arw, color="black") +
      scale_color_manual(values=c("forestgreen", "red"), drop=F) +
      # annotate(geom="text", 
      #          x=range(d$dst, na.rm=T) + diff(range(d$dst, na.rm=T)) / 20 * c(1,-1), 
      #          y=max(d$clim0 + .5), hjust=c(0, 1),
      #          label=c("climate tracking", "prevailing wind"), 
      #          color=c("black", "dodgerblue")) +
      theme_minimal() +
      theme(legend.position="none",
            strip.text=element_text(size=12)) +
      scale_y_reverse() +
      labs(y = "temperature (°C; note flipped scale)", 
           x = "distance along transect (km)")
#ggsave("figures/transects/transects.png", p, width=8, height=6, units="in")


climate <- list.files("data/cfsr_monthly/tmp2m", full.names=T)[1:24] %>%
      stack() %>%
      mean() %>%
      rotate() %>%
      mask(land) %>%
      rasterToPoints() %>%
      as.data.frame()

dp <- d %>% group_by(id) %>% arrange(dst) %>% slice(1)
#world <- map_data("world")
map <- ggplot() +
      #geom_polygon(data=world, aes(long, lat, group=group), fill="gray80") +
      geom_raster(data=climate, aes(x, y, fill = layer)) +
      geom_path(data=d, aes(x, y, group=id), color="black", size=.35,
                arrow=arrow(type="closed", angle=20, length=unit(.05, "in"))) +
      geom_point(data=dp, aes(x, y), color="black") +
      geom_text(data=dp, aes(x, y-7, label=i), color="black", size=4) +
      # scale_fill_gradientn(colors=c("gray97", "azure", "cyan", "dodgerblue", "blue",
      #                               "purple", "red", "orange", "yellow")) +
      scale_fill_gradientn(colors=c("black", "black", "black", "darkmagenta", "#6b49ff", 
                                    "dodgerblue", "turquoise", "#48c13f", 
                                    "gold", "red", "#660000")) +
      #theme_void() +
      theme(legend.position="none",
            axis.title.y=element_blank(),
            axis.title.x=element_text(color="white"), # all this whiteness makes it align with latitude plot
            axis.text=element_text(color="white"),
            axis.ticks=element_line(color="white"),
            panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
            panel.grid.minor.y=element_blank(), panel.grid.major.y=element_line(size=1),
            panel.background = element_rect(fill="gray90")) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(-90, 90), breaks=xb)
#ggsave("figures/transects/map.png", p, width=6, height=3, units="in")




### latitudinal patterns in meridional winds ###

modir <- "data/cfsr_monthly"
mf <- list.files(modir, full.names=T, recursive=T)

v <- mf[grepl("_v", mf)] %>% stack() %>% mean() %>% rotate() %>% stack(land)
v <- mf[grepl("tmp2m", mf)][1:24] %>% stack() %>% mean() %>% rotate() %>% stack(v) %>%
      rasterToPoints() %>% as.data.frame() %>%
      rename(temp = layer.1,
             v = layer.2) %>%
      mutate(latitude = abs(y),
             polarity = ifelse(y>0, v, -v),
             land = ifelse(is.na(land), "water", "land"))

vs <- v %>%
      bind_rows(mutate(v, land="all")) %>%
      mutate(lat = plyr::round_any(y, 5)) %>%
      group_by(land, lat) %>%
      summarize(v = mean(v),
                temp =  mean(temp)) %>%
      mutate(aligned_lat = sign(lat) == sign(v),
            dtemp = lead(temp) - lag(temp),
            aligned_temp = sign(v) != sign(dtemp),
            aligned_temp = ifelse(is.na(aligned_temp), aligned_lat, aligned_temp),
            
            xend = lat+sign(v)*15,
            xend = ifelse(xend < -90, -90, xend),
            xend = ifelse(xend > 90, 90, xend)  )

xb <- seq(-90, 90, 30)

bands <- data.frame(x0 = xb[1:6], x1 = xb[2:7],
                    name = c("Polar\ncell", "Ferrel\ncell", "Hadley\ncell")[c(1:3,3:1)],
                    winds = c("e", "p", "e", "e", "p", "e"))


mn <- -55
mx <- 30
latitude <- ggplot(vs %>% filter(land=="land")) +
      geom_rect(data=bands, aes(xmin=x0, xmax=x1, ymin=mn, ymax=mx, fill=winds),
                alpha=.15) +
      geom_ribbon(aes(lat, ymin=temp, ymax=mx), 
                  stat="identity", fill="gray90") +
      geom_vline(xintercept=xb, color="white", size=2) +
      geom_line(aes(lat, temp)) +
      geom_segment(aes(x=lat, xend=xend, y=temp, yend=temp,
                       color = aligned_temp),
                   arrow = arw, size=.25) +
      geom_point(aes(x=lat, y=temp, color = aligned_temp), size=2) +
      geom_text(data=bands, aes(x=(x0+x1)/2, y=mn, label=name), 
                hjust=.5, nudge_y=-8, lineheight=.7, size=3) +
      scale_color_manual(values=c("red", "forestgreen"), drop=F)  +
      scale_fill_manual(values=c("red", "forestgreen"), drop=F)  +
      scale_x_continuous(breaks=xb, limits=c(-90, 90), expand=c(0,0),
                         position="top") +
      theme_minimal() +
      theme(legend.position="none",
            panel.grid=element_blank()) +
      scale_y_reverse(limits=c(NA, mn), expand=c(0,0)) +
      labs(y = "temperature (°C; note flipped scale)", 
           x = "latitude") +
      coord_flip()

pc <- arrangeGrob(map, latitude, nrow=1, widths=c(2, 1))
pc <- arrangeGrob(pc, p, ncol=1, heights=c(1, 4/3))
ggsave("figures/transects/transects.png", pc, width=9, height=7, units="in")

