

library(windscape)
library(raster)
library(tidyverse)
library(geosphere)
library(grid)
library(gridExtra)


land <- raster("f:/cfsr/land.tif") %>% 
      rotate()

# load windrose data
rose <- stack("data/windrose/windrose_p1_2000s.tif") %>%
      rotate()

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)
names(rose) <- windrose_names()

# mean temperature
climate <- stack("data/geographic/processed/temperature.tif")
names(climate) <- c("clim0", "clim1")

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

d <- readRDS("data/transect_data.rds")

if(F){
      australia <- transect_plot(d$australia$transect)
      greenland <- transect_plot(d$greenland$transect)
      rockies <- transect_plot(d$rockies$transect)
      cali <- transect_plot(d$cali$transect)
      russia <- transect_plot(d$russia$transect)
      sahara <- transect_plot(d$sahara$transect)
      
      d <- list(rockies, cali, greenland, 
                australia, russia, sahara)
      names(d) <- c("rockies", "cali", "greenland",
                    "australia", "russia", "sahara")
      #saveRDS(d, "data/transect_data.rds")
}


labels <- c("Transverse mountain wind", "Sea breeze", "Katabatic wind",
            "Cross-desert wind", "Southwesterly", "Thermal low")
for(i in 1:length(d)) d[[i]]$data$id <- labels[i]

d <- lapply(d, function(x) x$data) %>%
      bind_rows() %>%
      mutate(id=factor(id, levels=labels))

plett <- c("c", "e", "g", "d", "f", "h")

d <- d %>%
      group_by(id) %>%
      mutate(ymax = max(clim0) + .1 * diff(range(clim0, na.rm=T)),
             i = as.integer(id),
             idi = paste0(i, ": ", id),
             idi2 = paste0("(", plett[i], ") ", id),
             scalar = diff(range(clim0, na.rm=T)) / max(abs(uwind), na.rm=T) / 5 ,
             yend = ifelse(aligned==FALSE, clim0 + uwind * scalar, clim0 - uwind * scalar),
             ymax = ifelse(ymax > yend, ymax, yend),
             ymax = max(ymax, na.rm=T),
             xend = dst + sign(uwind) * max(dst) / 50)

dp <- d %>% group_by(idi) %>% arrange(dst) %>% slice(1)
de <- d %>% group_by(idi) %>% arrange(desc(dst)) %>% slice(1)

ynudge <- 0
arw <- arrow(type="closed", angle=15, length = unit(0.1, "in"))
p <- ggplot(d) +
      facet_wrap(~id, scales="free") +
      geom_ribbon(aes(dst, ymin=clim0, ymax=ymax), 
                  stat="identity", fill="gray90") +
      #geom_point(aes(dst, clim0, color=aligned), size=2) +
      geom_segment(aes(x=dst, xend=xend,
                       y=clim0, yend=clim0, 
                       color=aligned),
                   arrow=arw, #color="blue", 
                   size=.3) +
      geom_line(aes(dst, ymax)) +
      geom_point(data=dp, aes(dst, ymax)) +
      geom_point(data=de, aes(dst, ymax), color="black", fill="white", shape=21) +
      #geom_segment(aes(x=dst, xend=dst+uclim,
      #                 y=clim0-ynudge, yend=clim0-ynudge),
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


clim <- climate[[1]] %>%
      rasterToPoints() %>%
      as.data.frame()

dp <- d %>% group_by(id) %>% arrange(dst) %>% slice(1)
de <- d %>% group_by(id) %>% arrange(desc(dst)) %>% slice(1)

x0 <- -113
y0 <- 39
r0 <- 15
x1 <- -130
y1 <- -25
r1 <- 45
circle_sm <- swfscMisc::circle.polygon(x0, y0, r0, poly.type="cartesian") %>%
      mapview::coords2Polygons(ID = "A")
circle_lg <- swfscMisc::circle.polygon(x1, y1, r1, poly.type="cartesian") %>%
      mapview::coords2Polygons(ID = "A")
clim_circ <- climate[[1]] %>%
      crop(circle_sm) %>% mask(circle_sm) %>%
      rasterToPoints() %>% as.data.frame() %>%
      mutate(x=((x-x0)/r0*r1)+x1, 
             y=((y-y0)/r0*r1)+y1)
d_circ <- d %>%
      filter(i %in% 1:2) %>%
      mutate(x=((x-x0)/r0*r1)+x1, 
             y=((y-y0)/r0*r1)+y1)
dp_circ <- d_circ %>% group_by(id) %>% arrange(dst) %>% slice(1)
de_circ <- d_circ %>% group_by(id) %>% arrange(desc(dst)) %>% slice(1)

# parameters for segment joining the circles
dx <- x1 - x0
dy <- y1 - y0
dh <- sqrt((dx)^2 + (dy)^2)
sz <- .25


xb <- seq(-90, 90, 30)


#world <- map_data("world")
map <- ggplot() +
      geom_raster(data=clim, aes(x, y, fill = clim0)) +
      geom_path(data=d %>% filter(! i %in% 1:2), aes(x, y, group=id), color="black", size=.35) +
      geom_point(data=dp %>% filter(! i %in% 1:2), aes(x, y), color="black") +
      geom_point(data=de %>% filter(! i %in% 1:2), aes(x, y), color="black", fill="white", shape=21) +
      geom_text(data=dp %>% filter(! i %in% 1:2), 
                aes(x, y-7, label=plett[i]), color="black", size=4) +
      
      annotate(geom="segment", size=sz, x=x1, y=y1, 
               xend=x0 + dx * (r0 / dh), 
               yend=y0 + dy * (r0 / dh)) +
      geom_polygon(data=fortify(circle_lg), aes(long, lat),
                   color="black", fill="white", size=sz) +
      geom_raster(data=clim_circ, aes(x, y, fill = clim0)) +
      geom_polygon(data=fortify(circle_sm), aes(long, lat),
                   color="black", fill=NA, size=sz) +
      geom_polygon(data=fortify(circle_lg), aes(long, lat),
                   color="black", fill=NA, size=sz) +
      geom_path(data=d_circ, aes(x, y, group=id), 
                color="black", size=sz) +
      geom_point(data=dp_circ, aes(x, y), color="black") +
      geom_point(data=de_circ, aes(x, y), color="black", fill="white", shape=21) +
      geom_text(data=dp_circ, aes(x, y-7, label=plett[i]), color="black", size=4) +
      
      
      
      scale_fill_gradientn(colors=c("black", "black", "black", "darkmagenta", "#6b49ff", 
                                    "dodgerblue", "turquoise", "#48c13f", 
                                    "gold", "red", "#660000")) +
      #theme_void() +
      labs(fill="°C") +
      guides(fill=guide_colorbar(barwidth=7, barheight=.3)) +
      theme(legend.position=c(.6, .22),
            legend.direction="horizontal",
            legend.background=element_blank(),
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
             
             xend = lat+sign(v)*5,
             xend = ifelse(xend < -90, -90, xend),
             xend = ifelse(xend > 90, 90, xend),
             
             yend = temp + (v)*sign(dtemp)*15)


bands <- data.frame(x0 = xb[1:6], x1 = xb[2:7],
                    name = c("Polar\ncell", "Ferrel\ncell", "Hadley\ncell")[c(1:3,3:1)],
                    winds = c("e", "p", "e", "e", "p", "e"))


mn <- -55
mx <- 35
arw <- arrow(type="closed", angle=15, length=unit(.1, "in"))
latitude <- ggplot(vs %>% filter(land=="land")) +
      geom_rect(data=bands, aes(xmin=x0, xmax=x1, ymin=mn, ymax=mx, fill=winds),
                alpha=.15) +
      geom_ribbon(aes(lat, ymin=temp, ymax=mx), 
                  stat="identity", fill="gray90") +
      geom_vline(xintercept=xb, color="white", size=1) +
      geom_line(aes(lat, temp), color="white", size=1) +
      # geom_segment(aes(x=lat, xend=lat, y=temp, yend=yend,
      #                  color = aligned_temp),
      #              arrow = arw, size=.25) +
      geom_segment(aes(x=lat, xend=xend, y=temp, yend=temp,
                       color = aligned_temp),
                   arrow = arw, size=.25) +
      #geom_point(aes(x=lat, y=temp, color = aligned_temp), size=1.5) +
      geom_text(data=bands, aes(x=(x0+x1)/2, y=mn, label=name), 
                hjust=.5, nudge_y=-8, lineheight=.7, size=3) +
      #annotate(geom="text", x=75, y=28, label="b", size=5) +
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

pc <- arrangeGrob(map, latitude, nrow=1, widths=c(1.5, 1))
pc <- arrangeGrob(pc, p, ncol=1, heights=c(1, 4/3))


source("E:/edges/range-edges/code/utilities.r")
ggs("figures/manuscript/fig_2.png", pc, width=9, height=7, units="in", 
    add=grid.text(letters[1:8], 
                  x=c(.05, .63,
                      .05, .05, 
                      .37, .37, 
                      .69, .69), 
                  y=c(.965, .965,
                      .54, .27, 
                      .54, .27, 
                      .54, .27),
                  gp=gpar(fontsize=20, fontface="bold", col="black")))

ggsave("figures/transects/transects.png", pc, width=9, height=7, units="in")





darkness <- theme(plot.background = element_rect(fill="black", color="black"),
                  text = element_text(color="white"),
                  strip.text = element_text(color="white"),
                  axis.text = element_text(color="white"))

pc <- arrangeGrob(map + darkness + 
                        theme(axis.text = element_text(color="black"),
                              axis.title = element_text(color="black"),
                              axis.ticks = element_line(color="black")), 
                  latitude + darkness + 
                        theme(panel.background = element_rect(fill="white")), 
                  nrow=1, widths=c(1.5, 1))
pc <- arrangeGrob(pc, 
                  p + darkness, 
                  ncol=1, heights=c(1, 4/3))
ggsave("figures/transects/transects_black.png", pc, width=9, height=7, units="in", bg="black")
