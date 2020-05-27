


library(windscape)
library(tidyverse)
library(raster)


# empirical examples ##################################################

samples <- function(infile, outdir, points, ncores=1){
      
      # file admin
      message(infile)
      outfile <- paste0(outdir, "/", basename(infile))
      outfile <- sub("grb2", "rds", outfile)
      if(any(file.exists(c(outfile,
                           sub("\\.rds", paste0("_", ncores, ".rds"), outfile))))){
            return("skipping")}
      
      # load and organize hourly wind data
      inv <- paste0(infile, ".inv") %>%
            readLines()
      process <- which(!grepl("6 hour fcst", inv))
      w <- infile %>%
            brick() %>%
            subset(process)
      even <- function(x) x %% 2 == 0
      w <- w[[c(which(!even(1:nlayers(w))),
                which(even(1:nlayers(w))))]]
      
      p <- extract(w, points)
      
      uv <- p %>%
            split(1:nrow(.)) %>%
            lapply(function(x) matrix(x, ncol=2, byrow=F))
      for(i in 1:length(uv)){
            x <- as.data.frame(uv[[i]])
            names(x) <- c("u", "v")
            x$x <- points$x[i]
            x$y <- points$y[i]
            uv[[i]] <- x
      }
      
      saveRDS(uv, outfile)
      
}


# hourly input data
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]

set.seed(123)
land <- raster("f:/cfsr/land.tif") %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(y > -60) %>%
      sample_n(1000) %>%
      dplyr::select(x, y)

map(f, samples, outdir="data/wind_regime/samples", points=land, ncores=6)

# select examples
ex <- c(601, 668, 818, 562)
ed <- d %>%
      filter(id %in% ex) %>%
      mutate(id = factor(id, levels = ex))
write_csv(ed, "figures/regimes/example_data.csv")


ed <- read_csv("figures/regimes/example_data.csv") %>%
      mutate(id = factor(id, levels = ex,
                         labels = c("\n[i]\nlow speeds &\nconsistent direction", "\n[ii]\nhigh speeds &\nconsistent direction", 
                                    "\n[iii]\nhigh speeds &\nvariable direction", "\n[iv]\nlow speeds &\nvariable direction")))


regime <- function(u, v, return="speed"){
      uv <- cbind(u, v)
      colnames(uv) <- NULL
      u <- mean(uv[,1])
      v <- mean(uv[,2])
      
      speeds <- sqrt(uv[,1]^2 + uv[,2]^2)
      speed <- mean(speeds)
      if(return=="speed") return(speed)
      
      dir <- do.call(atan2, as.list(colMeans(uv))) / pi * 180
      
      dirs <- atan2(uv[,1], uv[,2])
      x <- sum(sin(dirs) * speeds) / sum(speeds)
      y <- sum(cos(dirs) * speeds) / sum(speeds)
      
      rbar <- sqrt(x^2 + y^2) # mean resultant vector
      iso <- sqrt(1-rbar) # circular standard deviation
      
      if(return=="anisotropy") return(1 - iso)
}

centroids <- ed %>%
      group_by(id) %>%
      summarize(speed = regime(u, v, "speed"),
                aniso = regime(u, v, "anisotropy"),
                x = x[1], y = y[1],
                u = mean(u), v = mean(v)) %>%
      mutate(letter = c("i", "ii", "iii", "iv")) %>%
      mutate_at(vars(speed, aniso), signif, digits=2)

centroids <- sa %>% dplyr::select(speed, anisotropy, color) %>%
      rename(aniso = anisotropy) %>%
      mutate_at(vars(speed, aniso), signif, digits=2) %>%
      group_by(speed,aniso) %>% sample_n(1) %>%
      right_join(centroids) %>%
      mutate(dir = atan2(v, u))



rad <- 5
circle <- data.frame(angle=seq(0, pi*2, .01)) %>%
      mutate(x=cos(angle) * rad, 
             y=sin(angle) * rad)
acol <- "black"
lim <- 7.5


examples <- ed %>% filter(abs(u) < lim, abs(v) < lim) %>%
      ggplot(aes(u, v)) +
      facet_grid(.~id, space = "free", scales = "free") +
      geom_point(aes(color=id), alpha=.01, shape=16, size=1) +
      geom_path(data=circle, aes(x, y), color=acol) +
      geom_path(data=circle, aes(x*.6, y*.6), color=acol) +
      annotate(geom="segment", x=0, xend=0, y=-rad*.1, yend=rad*.1, color=acol) +
      annotate(geom="segment", x=-rad*.1, xend=rad*.1, y=0, yend=0, color=acol) +
      annotate(geom="segment",
               x=rad*c(0,1,-0,-1), xend=rad*c(0,1.15,-0,-1.15),
               y=rad*c(1,0,-1,0), yend=rad*c(1.15,0,-1.15,0),
               arrow=grid::arrow(type="closed", angle=15, length=unit(.1, "in")), color=acol) +
      geom_segment(data=centroids, aes(x=0, y=0, xend=cos(dir)*rad, yend=sin(dir)*rad)) +
      annotate(geom="text", x=0, y=rad * 1.23, label="N", color=acol, size=4) +
      annotate(geom="text", x=0, y=-rad * .85, label=paste0(rad, " m/s"), color=acol, size=4) +
      annotate(geom="text", x=0, y=-rad * .42, label=paste0(rad*.6, " m/s"), color=acol, size=4) +
      scale_color_manual(values=centroids$color) +
      theme_void(base_size = 18) + 
      theme(legend.position="none")





# legend #############################################################


land <- raster("data/land.tif")
w <- stack("data/regimes.tif") %>%
      mask(land) %>%
      rotate()
names(w) <- c("speed", "direction", "anisotropy", "u", "v")


# direction dataset
agg <- 10
wd <- w[[c("u", "v")]] %>%
      aggregate(agg)
drn <- wd[[1]]
drn[] <- atan2(values(wd)[,1], values(wd)[,2])
rsn <- mean(res(wd))
drn <- drn %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      rename(direction = u) %>% 
      mutate(xa = sin(direction) * rsn / 2,
             ya = cos(direction) * rsn / 2)


# speed and anisotropy dataset
sa <- w %>%
      rasterToPoints() %>%
      as.data.frame()

library(geosphere)
dg <- function(lat) distGeo(c(0, lat), c(1, lat))/1000
areas <- sa %>% dplyr::select(y) %>% distinct()
areas$area <- sapply(areas$y, dg)

sa <- sa %>%
      left_join(areas) %>%
      arrange(speed) %>%
      mutate(speedrnk = cumsum(area)) %>%
      arrange(anisotropy) %>%
      mutate(anisornk = cumsum(area))

sa$color <- colormap::colors2d(dplyr::select(sa, speedrnk, anisornk),
                               c("red", "yellow", "cyan", "blue"))

lt <- summarize_at(sa, vars(speed, anisotropy), funs(min, max)) %>%
      gather(stat, value) %>%
      separate(stat, c("var", "stat")) %>%
      spread(var, value)
lt <- expand.grid(anisotropy = lt$anisotropy,
                  speed = lt$speed)
lt <- lt[c(1, 2, 4, 3),]

library(ggtext)
legend <- ggplot(sa, aes(speed, anisotropy)) +
      geom_point(color=sa$color, size=.1) +
      geom_label(data=centroids, aes(x=speed, y=aniso, label=letter), size = 5) +
      theme_minimal(base_size = 16) + 
      scale_x_log10(breaks=c(1,2,5,10)) +
      scale_y_sqrt(breaks=c(.02, .2, .5), labels = c(".02", "0.2", "0.5")) +
      labs(x = "speed (m/s)",
           y = "<span style='font-size:16pt'><br>directional anisotropy</span>
                  <span style='color:white'>-------<br></span>
                  <span style='font-size:12pt'>variable   <--------->   consistent</span>
                  <span style='color:white'>--------</span>") +
      theme(axis.title.y = element_markdown())






####### map ###############################################################

library(tidyverse)
library(sf)        # for manipulation of simple features objects
library(lwgeom)    # for st_transform_proj()
library(rworldmap) # for getMap()

##### simulate particle trails ######

library(raster)

# load raster dataset with 2 layers: zonal and meiridional windspeeds
w <- stack("data/regimes.tif") %>% stack(raster("data/land.tif")) %>% rotate()
names(w) <- c("speed", "direction", "anisotropy", "u", "v", "land")


# latitude-based adjustment factor to correct for geodesic distortion
distortion <- function(lat){
      xd <- geosphere::distGeo(c(0, lat), c(.0001, lat))
      yd <- geosphere::distGeo(c(0, lat), c(0, lat - sign(lat) * .0001))
      yd / xd
}

# function to simulate particle movements across vector field
particle_trails <- function(d, # vector field rasters
                            np, # number of particles
                            nt, # number of timesteps
                            r, # smoothing neighborhood radius
                            scale){ # multiplier converting velocity to lat-lon
      #browser()
      x <- sort(unique(coordinates(d)[,1]))
      y <- desc(sort(unique(coordinates(d)[,2])))
      dx <- as.matrix(d[[1]])
      dy <- as.matrix(d[[2]])
      distort <- sapply(y, distortion)
      
      
      # sample starting locations, in proportion to area and 1/speed
      speed <- sqrt(sum(d^2))
      lat <- speed
      lat[] <- coordinates(lat)[,2]
      weight <- cos(lat/180*pi) / speed
      start <- coordinates(d)[sample(1:ncell(speed), np, replace=T, 
                                     prob = values(weight)),]
      
      # add some noise
      rs <- mean(res(d)) / 2
      start <- start +  runif(nrow(start)*2, -rs, rs)
      
      ptx <- pty <- matrix(NA, nt, np)
      for(p in 1:np){
            px <- start[p,1]
            py <- start[p,2]
            
            for(t in 1:nt){
                  xi <- which(x > px - r & x < px + r)
                  yi <- which(y > py - r & y < py + r)
                  if(length(xi)==0 | length(yi)==0) browser()
                  px <- px + mean(dx[yi, xi]) * scale * mean(distort[yi])
                  py <- py + mean(dy[yi, xi]) * scale
                  ptx[t, p] <- px
                  pty[t, p] <- py
                  
                  if(py > 90 | py < -90) break()
                  if(px > 180) px <- px - 360
                  if(px < -180) px <- px + 360
                  if(px > 180 | px < -180) break()
            }
      }
      
      frame <- function(pd){
            colnames(pd) <- paste0("p", 1:ncol(pd))
            as.data.frame(pd) %>%
                  mutate(time = 1:nrow(.)) %>%
                  gather(id, z, -time)
      }
      
      left_join(frame(ptx) %>% rename(x = z),
                frame(pty) %>% rename(y = z))
}

if(F){
      f <- particle_trails(d = w[[c("u", "v")]], np = 100000, nt = 50, r = 1, scale = .1)
      saveRDS(f, "data/data_speedsampling.rds")
}


##### plotting ####

world <- getMap(resolution = "low") %>% rgeos::gBuffer(width = 0) %>% rgeos::gUnaryUnion()
world_sf <- st_as_sf(world)

crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
world_wintri <- st_transform_proj(world_sf, crs = crs_wintri)

grat_wintri <- 
      st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
      st_transform_proj(crs = crs_wintri)

# vectors of latitudes and longitudes that go once around the 
# globe in 1-degree steps
lats <- c(90:-90, -90:90, 90)
longs <- c(rep(c(180, -180), each = 181), 180)

# turn into correctly projected sf collection
wintri_outline <- 
      list(cbind(longs, lats)) %>%
      st_polygon() %>%
      st_sfc( # create sf geometry list column
            crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      ) %>% 
      st_sf() %>%
      st_transform_proj(crs = crs_wintri) # transform to Winkel tripel



f <- readRDS("data/data_speedsampling.rds")

ff <- f %>%
      filter(abs(x)<180, abs(y)<90) %>%
      group_by(id) %>%
      mutate(span = any(between(x, -180, -150)) &
                   any(between(x, 150, 180)),
             sid = case_when(span & x > 0 ~ id,
                             span & x < 0 ~ paste0(id, "b"),
                             TRUE ~ id))

fwt <- ff
coordinates(fwt) <- c("x", "y")
we <- raster::extract(w, fwt)
lnd <- raster::extract(land %>% rotate(), fwt)
fwt <- st_as_sf(fwt)
st_crs(fwt) <- st_crs(world_sf)
fwt <- st_transform_proj(fwt, crs = crs_wintri)



fc <- st_coordinates(fwt) %>% as.data.frame()
fwt <- bind_cols(fwt, fc)
fwt <- bind_cols(fwt, as.data.frame(we))
fwt <- bind_cols(fwt, data.frame(land = lnd))



fwt <- fwt %>% mutate(anisoland = anisotropy * land,
                      speedland = speed * land)

wf <- which(is.finite(fwt$speedland))
fwt$color <- "white"
fwt$color[wf] <- colormap::colors2d(cbind(log10(fwt$speedland[wf]), 
                                          sqrt(fwt$anisoland[wf])),
                                    colors=c("red", "yellow", "cyan", "blue"),
                                    xtrans="rank", ytrans="rank")




### arrows

arw <- expand_grid(x = seq(-165, 165, 30),
                   y = seq(-75, 75, 15))
ffa <- ff %>%
      mutate(xr = round(x),
             yr = round(y))

arw <- ffa %>%
      filter(paste(xr, yr) %in% paste(arw$x, arw$y)) %>%
      group_by(id) %>%
      mutate(n = n()) %>%
      #filter(n > 1) %>%
      filter(time != 50) %>%
      group_by(xr, yr) %>%
      arrange(id) %>%
      slice(1)

arw <- left_join(ffa, arw) %>%
      group_by(id) %>%
      filter(any(is.finite(n))) %>%
      filter(is.finite(n) |
                   time == time[is.finite(n)]+1) %>%
      mutate(n = n()) %>%
      filter(n == 2) %>%
      mutate(x = c(x[1], x[1] + (x[2] - x[1])/1000),
             y = c(y[1], y[1] + (y[2] - y[1])/1000))

coordinates(arw) <- c("x", "y")
land <- getMap(resolution = "low") %>% rgeos::gBuffer(width = 0) %>% rgeos::gUnaryUnion()
crs(arw) <- crs(land)

arw <- st_as_sf(arw)
st_crs(arw) <- st_crs(world_sf)
arw <- st_transform_proj(arw, crs = crs_wintri)
arw$land <- st_distance(world_wintri, arw)[1,]
arw <- bind_cols(arw, as.data.frame(st_coordinates(arw)))


### plot ###

map_particles <- ggplot() + 
      geom_sf(data = wintri_outline, fill = "gray20", color = NA) + # ocean
      geom_sf(data = world_wintri, color = "black", fill="black", size = 0.5/.pt) + # land
      geom_path(data = fwt, aes(X, Y, group=sid, alpha=time), color=fwt$color, size=.1) + # wind
      geom_sf(data = world_wintri, color = "black", fill=NA, size = .25) + # coastlines
      geom_path(data = arw, aes(X, Y, group = id), 
                arrow = arrow(type = "closed", length = unit(3, "mm"), angle = 10), 
                color = "white", size = .25) +
      coord_sf(datum = NULL, expand=F) +
      theme_void() +
      theme(legend.position="none")





# final composite figure ##################################################

devtools::source_url("https://raw.githubusercontent.com/matthewkling/range-edges/master/code/utilities.R")
library(gridExtra)
library(grid)

p <- arrangeGrob(examples, legend, nrow=1, widths=c(3, 1.4))
p <- arrangeGrob(p, map_particles, ncol=1, heights=c(.825, 2))

ggs("figures/manuscript/fig_1.png", p, width=12, height=10.5, units="in", dpi=1000,
    add=grid.text(letters[1:3], 
                  x=c(.02, .72, .02), 
                  y=c(.985, .985, .69),
                  gp=gpar(fontsize=25, fontface="bold", col="black")))

