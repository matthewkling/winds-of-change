


library(raster)
library(tidyverse)
library(rNOMADS)
library(grid)


# process CFSR data:

# 1. calculate monnhly means of U and V wind components
# 2. calculate long-term means
# 2. calculate wind direction


# binary land-ocean layer
land <- raster("f:/cfsr/soilm1.gdas.197901.grb2") %>% 
      reclassify(c(-Inf, Inf, 1)) %>%
      writeRaster("f:/cfsr/land.tif")


modir <- "data/cfsr_monthly"



#### calculate monthly means from raw hourly data

# wind
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
means <- function(file){
      #file <- f[1]
      message(basename(file))
      out <- paste0(modir, "/wnd10m/", sub("grb2", "tif", basename(file)))
      if(file.exists(sub("\\.tif", "_v.tif", out))) return("pass")
      w <- stack(file)
      even <- function(x) x %% 2 == 0
      u <- w[[which(!even(1:nlayers(w)))]] %>% calc(mean, filename=sub("\\.tif", "_u.tif", out))
      v <- w[[which(even(1:nlayers(w)))]] %>% calc(mean, filename=sub("\\.tif", "_v.tif", out))
}
if(T) lapply(f[241:360], means)

# temperature
f <- list.files("f:/CFSR/tmp2m", full.names=T)
f <- f[!grepl("inv", f)]
means <- function(file){
      #file <- f[1]
      message(basename(file))
      out <- paste0(modir, "/tmp2m/", sub("grb2", "tif", basename(file)))
      if(file.exists(out)) return("pass")
      w <- stack(file) %>% calc(mean, filename=out)
}
if(F) lapply(f, means)






# load monthly summaries of CFSR generated above
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

stack(u, v, wm, wa) %>% writeRaster("data/geographic/processed/wind_uvma.tif")


# climate change
gradient <- function(x, component="x"){
      y_cellsize <- 0.8 # PRISM quotes latitude as its raster size
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
      
      cr <- function(x) crop(rotate(x), ta)
      s <- stack(cr(u), cr(v), cr(wa), cr(wm), 
                 ta) #vx, vy, ta, tm)
      names(s) <- c("wx", "wy", "wa", "wm", 
                    "ta")# "tx", "ty", "ta", "tm")
      
      land <- raster("f:/cfsr/land.tif")
      s <- stack(s, cr(land))
      
      writeRaster(s, paste0("data/alignment/alignment_s", neigh, ".tif"), overwrite=T)
      
      
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
      grid.text("headwind", hjust=0.5, x=.07, y=.09+nudgey, gp=gpar(col=pal[5], cex=3.5, fontface="bold"))
      dev.off()
}




### wind direction figure ###

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


### composite figure for manuscript ####

library(png)
p1 <- readPNG("figures/tailwinds/wind_direction.png") %>% rasterGrob(interpolate=TRUE)
p2 <- readPNG("figures/tailwinds/temperature_angle_s3.png") %>% rasterGrob(interpolate=TRUE)
p3 <- readPNG("figures/tailwinds/alignment_s3.png") %>% rasterGrob(interpolate=TRUE)

pad <- rectGrob(gp=gpar(fill="white", col="white"))
p <- arrangeGrob(p1, pad, p2, nrow=1, widths=c(60,1,60))
p <- arrangeGrob(p, pad, p3, ncol=1, heights=c(40, 1, 80))

source("E:/edges/range-edges/code/utilities.r")
ggs("figures/manuscript/SI_fig_alignment.png",
    p, width=12, height=9.1, units="in",
    add = grid.text(letters[1:3], 
                    x=c(.03, .53, .03), 
                    y=c(.92, .92, .52),
                    gp=gpar(fontsize=30, fontface="bold", col="white")))
