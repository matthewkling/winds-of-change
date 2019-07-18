


library(raster)
library(tidyverse)
library(rNOMADS)
library(grid)



# monthly summaries of CFSR
modir <- "data/cfsr_monthly"
mf <- list.files(modir, full.names=T, recursive=T)


# long-term means
t <- mf[grepl("tmp2m", mf)][1:24] %>% stack() %>% mean()
u <- mf[grepl("_u", mf)][1:72] %>% stack() %>% mean()
v <- mf[grepl("_v", mf)][1:72] %>% stack() %>% mean()


# convert wind components to speed and direction
w <- MagnitudeAzimuth(values(u), values(v))
wm <- wa <- u
wm[] <- w$magnitude
wa[] <- w$azimuth


# climate change
gradient <- function(x, component="x"){
      y_cellsize <- 0.8 # guessing maybe PRISM quotes latitude as its raster size
      x_cellsize <- y_cellsize * 78846.81 / 111131.75 # this is the ratio at 45N
      browser()
      dzdx <- ((x[3] + 2*x[6] + x[9]) - (x[1] + 2*x[4] + x[7])) / (8 * x_cellsize)
      dzdy <- ((x[7] + 2*x[8] + x[9]) - (x[1] + 2*x[2] + x[3])) / (8 * y_cellsize)
      if(component=="x"){return(dzdx)}
      if(component=="y"){return(dzdy)}
}

grangle <- function(x){
      n <- sqrt(length(x))
      if(is.na(x[(length(x)+1)/2])) return(NA)
      x <- matrix(x, n)
      m <- lm(y ~ x, data = data.frame(y = apply(x, 1, mean, na.rm=T), x = (1:n)))
      dx <- coef(m)["x"]
      m <- lm(y ~ x, data = data.frame(y = apply(x, 2, mean, na.rm=T), x = (1:n)))
      dy <- coef(m)["x"]
      MagnitudeAzimuth(dx, dy)["azimuth"]
}


neighborhoods <- c(3, 9, 27)
for(neigh in neighborhoods){
      
      land <- raster("f:/cfsr/land.tif")
      
      wmat <- matrix(1, neigh, neigh)
      
      ta <- focal(mask(t, land), w=wmat, fun=grangle) %>% rotate()
      taa <- focal(mask(t, land) %>% rotate(), w=wmat, fun=grangle)
      ta <- merge(ta, taa)
      
      
      # vx1 <- rotate(-1 * focal(t, w=wmat, fun=function(x) gradient(x, component="x")))
      # vx2 <- -1 * focal(rotate(t), w=wmat, fun=function(x) gradient(x, component="x"))
      # vx <- merge(vx1, vx2)
      # vy1 <- rotate(focal(t, w=wmat, fun=function(x) gradient(x, component="y")))
      # vy2 <- focal(rotate(t), w=wmat, fun=function(x) gradient(x, component="y"))
      # vy <- merge(vy1, vy2)
      # vx <- trim(vx)
      # vy <- crop(vy, vx)
      # tv <- MagnitudeAzimuth(values(vx), values(vy))
      # tm <- ta <- vx
      # tm[] <- log10(1 / tv$magnitude)
      # ta[] <- tv$azimuth
      
      cr <- function(x) crop(rotate(x), ta)
      s <- stack(cr(u), cr(v), cr(wa), cr(wm), 
                 ta) #vx, vy, ta, tm)
      names(s) <- c("wx", "wy", "wa", "wm", 
                    "ta")# "tx", "ty", "ta", "tm")
      
      land <- raster("f:/cfsr/land.tif")
      s <- stack(s, cr(land))
      
      
      
      ###########
      
      # project to hobo-dyer equal area
      #s <- projectRaster(s, crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
      
      #continents <- map_data("world")
      #states <- map_data("state")
      
      
      ########################
      
      
      d <- s %>% 
            rasterToPoints %>% 
            as.data.frame()
      d <- d %>%
            mutate(theta = abs(wa - ta),
                   theta = ifelse(theta > 180, 360-theta, theta))
      
      d$tm <- 1
      
      dl <- filter(d, is.finite(land))
      
      
      
      compass <- ggplot() +
            geom_point(data=data.frame(x=0:360), 
                       aes(x, .9, color=x),
                       size=5) +
            geom_segment(data=data.frame(x=seq(0, 360, 45)), 
                         aes(x=x, xend=x, y=.5, yend=1.25, color=x),
                         arrow=arrow(type="closed", angle=12,
                                     length=unit(0.5, "inches")),
                         size=2) +
            coord_polar() +
            scale_x_continuous(breaks=seq(0, 360, 45),
                               limits=c(0, 360)) +
            scale_color_gradientn(colors=rainbow(20)) +
            theme_minimal() +
            theme(legend.position="none",
                  axis.text.y=element_blank(),
                  axis.text.x=element_blank(),
                  panel.grid=element_blank(),
                  axis.title=element_blank())
      
      
      map <- ggplot() +
            geom_raster(data=dl, 
                        aes(x, y), fill="white") +
            geom_raster(data=filter(d, is.finite(land)), 
                        aes(x, y, fill=ta#, alpha=tm
                            )) +
            scale_alpha_continuous(trans="log", range=0:1) +
            scale_fill_gradientn(colors=rainbow(20)) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            theme_void() +
            theme(legend.position="none",
                  panel.background=element_rect(fill="black"))
      png(paste0("figures/tailwinds/temperature_angle_s", neigh, ".png"), width=2000, height=1000)
      plot(map)
      plot(compass, vp=viewport(x=.13, y=.33, width=.25, height=.35))
      dev.off()
      
      
      
      da <- dl %>%
            mutate(tm=ifelse(tm > 1, 1, tm),
                   tm=ifelse(tm < -.8, -.8, tm))
      ld <- expand.grid(theta=seq(0, 180, length.out=20)-360/80, 
                        velocity=seq(min(da$tm), max(da$tm), length.out=20))
      pal <- c("dodgerblue", "blue", "purple", "red", "orange")
      col_v0 <- "white"
      tach <- ggplot(ld, aes(theta, velocity, 
                             fill=theta, alpha=velocity)) +
            geom_tile(fill=col_v0, alpha=1) +
            geom_tile() +
            annotate(geom="segment", 
                     x=60, xend=60, y=min(da$tm), yend=max(da$tm)*.7,
                     arrow=arrow(type="closed", angle=15), size=1) +
            annotate(geom="segment", 
                     x=59, xend=60, y=mean(da$tm), yend=mean(da$tm),
                     arrow=arrow(type="closed", angle=15)) +
            annotate(geom="segment", 
                     x=1, xend=0, y=mean(da$tm), yend=mean(da$tm),
                     arrow=arrow(type="closed", angle=15)) +
            annotate(geom="line", size=1,
                     x=0:60, y=mean(da$tm)) +
            annotate(geom="text",
                     x=c(30, 75), y=c(mean(da$tm)*2.2, 0),
                     label=c(expression("Theta"), "v"), 
                     size=14, fontface="bold", parse=T) +
            scale_x_continuous(limits=c(0, 360),
                               breaks=seq(0, 180, 45),
                               labels = function(x) paste0(x, "°")) +
            scale_fill_gradientn(colours=pal,
                                 breaks=seq(0, 180, 45), limits=c(0, 180)) +
            coord_polar() +
            theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text.x=element_text(color="white", size=25),
                  #axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  plot.background=element_blank(),
                  panel.background=element_blank(),
                  panel.grid=element_blank())
      
      brk <- data.frame(x=c(0, 45, 90, 135, 180),
                        y=1.5)
      legend <- ggplot() +
            geom_segment(data=brk, aes(x=x, xend=x, 
                     y=0, yend=1, color=x), size=2) +
            geom_text(data=brk, aes(x, y, label=paste0(x, "°"), 
                                    color=x),
                       size=15) +
            geom_point(data=data.frame(x=seq(0, 180, 1)),
                       aes(x, 1, color=x), size=10) +
            coord_polar() +
            xlim(0, 360) +
            ylim(0, NA) +
            scale_color_gradientn(colors=(pal)) +
            theme_void() +
            theme(legend.position="none")
      
      map <- ggplot() +
            geom_raster(data=dl, 
                        aes(x, y), fill=col_v0) +
            geom_raster(data=da, 
                        aes(x, y, fill=theta#, alpha=tm
                            )) +
            scale_fill_gradientn(colours=pal,
                                 breaks=seq(0, 180, 45), limits=c(0, 180)) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            theme_void() +
            theme(legend.position="none",
                  panel.background=element_rect(fill="black"))
      png(paste0("figures/tailwinds/alignment_s", neigh, ".png"), width=2000, height=1000)
      nudgey <- .08
      plot(map)
      plot(legend, vp=viewport(x=.07, y=.33+nudgey, width=.3, height=.5))
      grid.text("tailwind", hjust=0.5, x=.07, y=.58+nudgey, gp=gpar(col=pal[1], cex=3.5, fontface="bold"))
      #grid.text("crosswind", hjust=0, x=.11, y=.335+nudgey, gp=gpar(col=pal[3], cex=3.5, fontface="bold"))
      grid.text("headwind", hjust=0.5, x=.07, y=.09+nudgey, gp=gpar(col=pal[5], cex=3.5, fontface="bold"))
      dev.off()
}


###


map <- ggplot() +
      geom_raster(data=dl, 
                  aes(x, y), fill="white") +
      geom_raster(data=filter(d, is.finite(land)), 
                  aes(x, y, fill=wa, alpha=wm)) +
      scale_alpha_continuous(trans="log", range=0:1) +
      scale_fill_gradientn(colors=rainbow(20)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_void() +
      theme(legend.position="none",
            panel.background=element_rect(fill="black"))
png("figures/tailwinds/wind_direction.png", width=2000, height=1000)
plot(map)
plot(compass, vp=viewport(x=.13, y=.33, width=.25, height=.35))
dev.off()


p <- ggplot(dl %>% sample_n(10000), 
            aes(wa, wm, color=wa)) +
      geom_point(shape=16) +
      coord_polar() +
      scale_x_continuous(breaks=seq(0, 360, 45),
                         limits=c(0, 360)) +
      ylim(-500, NA) +
      scale_color_gradientn(colors=rainbow(20)) +
      theme_minimal() +
      theme(legend.position="none") +
      labs("mean wind azimuth and magnitude (terrestrial only)")
ggsave("figures/tailwinds/scatter_azimuth_magnitude.png", p, width=8, height=8, units="in")

p <- ggplot(d %>% sample_n(10000), 
            aes(wa, wm, color=wa)) +
      geom_point(shape=16) +
      coord_polar() +
      scale_x_continuous(breaks=seq(0, 360, 45),
                         limits=c(0, 360)) +
      ylim(-500, NA) +
      scale_color_gradientn(colors=rainbow(20)) +
      theme_minimal() +
      theme(legend.position="none") +
      labs("mean wind azimuth and magnitude (terrestrial and ocean)")
ggsave("figures/tailwinds/scatter_azimuth_magnitude_all.png", p, width=8, height=8, units="in")


map <- ggplot() +
      geom_raster(data=d, 
                  aes(x, y, fill=wa, alpha=wm)) +
      #geom_path(data=continents, 
      #          aes(long, lat, group=group),
      #          color="white", size=.1) +
      scale_alpha_continuous(trans="sqrt", range=0:1) +
      scale_fill_gradientn(colors=rainbow(20)) +
      scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
      scale_y_continuous(expand=c(0,0), limits=c(-80, NA)) +
      theme_void() +
      theme(legend.position="none")
png("figures/tailwinds/wind_direction_full.png", width=2000, height=1000)
plot(map)
plot(compass, vp=viewport(x=.13, y=.33, width=.25, height=.35))
dev.off()





ggplot() +
      geom_raster(data=d, aes(x, y, fill=theta)) +
      geom_path(data=continents, 
                aes(long, lat, group=group),
                color="white", size=.1) +
      geom_path(data=states, 
                aes(long, lat, group=group),
                color="white", size=.1) +
      coord_cartesian(xlim=c(-125, -100),
                      ylim=c(25, 50)) +
      scale_fill_viridis_c() +
      theme_void()



ggplot() +
      geom_raster(data=d %>%
                        filter(is.finite(land)) %>%
                        select(x, y, wa, ta) %>%
                        gather(stat, value, wa, ta), 
                  aes(x, y, fill=value)) +
      geom_path(data=continents, 
                aes(long, lat, group=group),
                color="black", size=.1) +
      geom_path(data=states, 
                aes(long, lat, group=group),
                color="black", size=.1) +
      facet_wrap(~stat, nrow=1) +
      coord_cartesian(xlim=c(-125, -100),
                      ylim=c(25, 50)) +
      scale_fill_gradientn(colors=rainbow(20)) +
      theme_void()


ggplot(dl %>% sample_n(5000), 
       aes(wm, tm, color=theta)) +
      geom_point() +
      geom_smooth(method=lm, se=F, color="black") +
      scale_x_log10() +
      scale_color_gradientn(colours=c("forestgreen", "darkred"))


# percent of grid cells >90 degrees
mean(d$theta[is.finite(d$land)]> 90)

mean(dl$theta > 135) / mean(dl$theta < 45)


dlecdf <- dl %>% mutate(thetaq = ecdf(theta)(theta),
                        theta = round(theta, 2)) %>%
      filter(theta %in% c(45, 135)) %>%
      group_by(theta) %>%
      sample_n(1) %>%
      arrange(theta)

p <- ggplot(dl, aes(theta)) +
      geom_rect(data=dlecdf, 
                aes(xmin=c(0, 180), xmax=theta, 
                    ymin=c(0,1), ymax=thetaq),
                alpha=.5) +
      geom_abline(slope=1/180, intercept=0, linetype=2) +
      stat_ecdf() +
      scale_x_continuous(breaks=seq(0, 180, 45)) +
      scale_y_continuous(breaks=seq(0, 1, .1)) +
      theme_minimal()
ggsave("figures/tailwinds/ecdf_theta.png", p, width=8, height=8, units="in")


p <- ggplot(dl, 
            aes(y, theta)) +
      geom_point(size=.05) + 
      geom_smooth(se=F, color="red") +
      scale_y_continuous(breaks=seq(0, 180, 45), expand=c(0,0)) +
      scale_x_continuous(breaks=c(-90, -67, -23.5, 0, 23.5, 67, 90), expand=c(0,0)) +
      theme_minimal() +
      labs(x="latitude")
ggsave("figures/tailwinds/scater_latitude_theta.png", p, width=8, height=8, units="in")
