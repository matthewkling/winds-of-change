

library(windscape)
library(tidyverse)
library(raster)


# derive wind regime stats ##########################################################

regime <- function(x, ...){
      uv <- matrix(x, ncol=2, byrow=F)
      u <- mean(uv[,1])
      v <- mean(uv[,2])
      
      speeds <- sqrt(uv[,1]^2 + uv[,2]^2)
      speed <- mean(speeds)
      dir <- do.call(atan2, as.list(colMeans(uv))) / pi * 180
      
      dirs <- atan2(uv[,1], uv[,2])
      x <- sum(sin(dirs) * speeds) / sum(speeds)
      y <- sum(cos(dirs) * speeds) / sum(speeds)
      
      rbar <- sqrt(x^2 + y^2) # mean resultant vector
      iso <- sqrt(1-rbar) # circular standard deviation
      
      c(speed = speed, dir = dir, 
        rbar = rbar, aniso = 1 - iso,
        u = u, v = v, x = x, y = y)
}

windregime <- function(infile, outdir, ncores=1){
      
      # file admin
      message(infile)
      outfile <- paste0(outdir, "/", basename(infile))
      outfile <- sub("grb2", "tif", outfile)
      if(any(file.exists(c(outfile,
                           sub("\\.tif", paste0("_", ncores, ".tif"), outfile))))){
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
      
      if(ncores == 1){
            wr <- raster::calc(w, fun=regime, forceapply=TRUE, filename = outfile)
            return(wr)
      } else {
            
            # split dataset into batches
            nlayer <- nlayers(w)/2
            core <- rep(1:ncores, each=floor(nlayer/ncores))
            core <- c(core, rep(ncores, nlayer %% ncores))
            w <- lapply(1:ncores, function(x) list(i=x,
                                                   data=w[[c(which(core == x),
                                                             which(core == x) + nlayer)]]))
            
            # process batches in parallel
            require(doParallel)
            cl <- makeCluster(ncores)
            registerDoParallel(cl)
            wr <- foreach(x = w,
                          .export = "regime",
                          .packages=c("raster", "windscape")) %dopar% {
                                raster::calc(x$data, fun=regime, forceapply=TRUE, 
                                             filename=paste0(substr(outfile, 1, nchar(outfile)-4),
                                                             "_", x$i,
                                                             substr(outfile, nchar(outfile)-3, nchar(outfile))))
                          }
            stopCluster(cl)
      }
}


# hourly input data
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
f <- rev(f)
map(f, windregime, outdir="data/wind_regime/monthly", ncores=6)



# compile long-term regimes from monthlies

mean_regime <- function(x){
      
      # average speed
      speed <- lapply(x, function(y) stack(y)[[1]]) %>%
            Reduce("+", .) %>%
            "/"(length(x))
      
      # prevailing direction
      u <- lapply(x, function(y) stack(y)[[5]]) %>%
            Reduce("+", .) %>%
            "/"(length(x))
      v <- lapply(x, function(y) stack(y)[[6]]) %>%
            Reduce("+", .) %>%
            "/"(length(x))
      dir <- u
      dir[] <- atan2(values(u), values(v)) / pi * 180
      
      # anisotropy (mean of partial mean resultant vectors, weighted by speed of each)
      ss <- lapply(x, function(y) stack(y)[[1]]) %>%
            Reduce("+", .)
      xs <- lapply(x, function(y) stack(y)[[7]] * stack(y)[[1]]) %>%
            Reduce("+", .) %>%
            "/"(ss)
      ys <- lapply(x, function(y) stack(y)[[8]] * stack(y)[[1]]) %>%
            Reduce("+", .) %>%
            "/"(ss)
      rbar <- sqrt(xs^2 + ys^2)
      aniso <- 1 - sqrt(1 - rbar)
      
      stack(speed, dir, aniso, u, v)
}

wr <- list.files("data/wind_regime/monthly", full.names=T) %>%
      mean_regime() %>%
      writeRaster("data/wind_regime/regimes.tif", overwrite=T)





# global wind regime plots ###################################################

land <- raster("f:/cfsr/land.tif")
w <- stack("data/wind_regime/regimes.tif") %>%
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
sa$color <- colormap::colors2d(dplyr::select(sa, speed, anisotropy),
                                 c("magenta", "cyan", "darkblue", "darkred"),
                                 xtrans="rank", ytrans="rank")

map <- ggplot() +
      geom_raster(data=sa, aes(x, y), fill = sa$color) +
      geom_segment(data=drn, aes(x=x+xa*.999, xend=x+xa, y=y+ya*.999, yend=y+ya),
                 arrow = arrow(type="open", angle=10, length=unit(.1, "in")),
                 color = "white", size = .5) +
      #geom_raster(data=ocean, aes(x, y), fill = "white") +
      theme_void() +
      coord_cartesian(xlim=c(-180, 180),
                      ylim=c(-90, 90),
                      expand = 0)

legend <- ggplot(sa, aes(speed, anisotropy)) +
      geom_point(color=sa$color, size=.1) +
      theme_minimal() +
      scale_x_log10(breaks=c(1,2,5,10)) +
      ylim(0:1) +
      #theme(text=element_text(size = 25)) +
      labs(x = "speed (m/s)")



# wind regime cartoons ####################################################

library(tidyverse)
library(grid)
library(circular)
library(windscape)
library(geosphere)

n <- 1000
mid <- pi/4
spread <- 1.5
speed_spread <- .5

circ_diff <- function(x, y){
      d <- abs(x - y)
      d[d > pi] <- 2*pi - d[d > pi]
      d
}

ref <- data.frame(angle = rvonmises(n, mid, spread) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid, speed_spread) * runif(n),
             speed = speed/mean(speed),
             panel = "ref")

stren <- data.frame(angle = rvonmises(n, mid, spread) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid, speed_spread) * runif(n),
             speed = speed/mean(speed) * .5,
             panel = "stren")

iso <- data.frame(angle = rvonmises(n, mid, spread * 4) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid, speed_spread * 4) * runif(n) / 2,
             speed = speed/mean(speed),
             panel = "iso")

dir <- data.frame(angle = rvonmises(n, mid + 1.5, spread) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid + 1.5, speed_spread) * runif(n),
             speed = speed/mean(speed),
             panel = "dir")


d <- bind_rows(ref, dir, stren, iso) %>%
      mutate(panel = factor(panel, levels=c("ref", "dir", "iso", "stren"))) %>%
      mutate(y = cos(angle) * speed,
             x = sin(angle) * speed)


mx <- max(d$speed) * .5

circle <- data.frame(angle=seq(0, pi*2, .01)) %>%
      mutate(x=cos(angle) * mx, 
             y=sin(angle) * mx)

txt <- data.frame(panel = factor(c("ref", "dir", "iso", "stren"), 
                                 levels=c("ref", "dir", "iso", "stren")),
                  x = 0, y = -1.25 * mx,
                  label = c("a local\nwind regime", "direction\ndifference", 
                            "directionality\ndifference", "strength\ndifference"))

brk <- seq(-1, 1, .1)

eb <- element_blank()

aniso <- as.character(round(1 - circ_sd(ref$angle / pi * 180, ref$speed)["iso"], 2))
aniso2 <- as.character(round(1 - circ_sd(iso$angle / pi * 180, iso$speed)["iso"], 2))

stats <- data.frame(panel = c("ref", "dir", "iso", "stren"),
                    x = 0, y = -2,
                    direction = c("45°", "135°", "45°", "45°"),
                    anisotropy = c(aniso, aniso, aniso2, aniso),
                    speed = c("5 m/s", "5 m/s", "5 m/s", "2.5 m/s")) %>%
      gather(stat, value, direction:speed) %>%
      mutate(label = paste0(stat, " = ", value),
             y = case_when(stat == "direction" ~ y,
                           stat == "anisotropy" ~ y - .25,
                           TRUE ~ y - .5),
             color = case_when(panel == "ref" ~ "red",
                               panel == "dir" & stat == "direction" ~ "red",
                               panel == "iso" & stat == "anisotropy" ~ "red",
                               panel == "stren" & stat == "speed" ~ "red",
                               TRUE ~ "black"))

speeds <- data.frame(panel = "ref", x = 0, y = c(-.5, -1) * mx,
                     label = c("2.5 m/s", "5 m/s"))

segment <- data.frame(panel = factor("ref", levels=levels(d$panel)))

cartoons <- ggplot(d) +
      facet_grid(.~panel) +
      geom_path(data=circle, aes(x, y), color="red") +
      geom_path(data=circle, 
                aes(x/2, y/2), color="red") +
      geom_segment(aes(x=x*.999, xend=x, y=y*.999, yend=y),
                   arrow = arrow(type="closed", angle=10, length=unit(.1, "in")),
                   alpha = .25) +
      geom_segment(data=segment, x=0, xend=0, y=0, yend=.7 * mx, color="red") +
      annotate(geom="point", x=0, y=0, color="red", size=3) +
      shadowtext::geom_shadowtext(data=stats, aes(x, y, label=label, alpha=color), 
                                  color="black", fontface="bold", lineheight=.75, bg.colour='white') +
      shadowtext::geom_shadowtext(data=segment, x=0, y=.85 * mx, label="N", 
                                  color="red", bg.colour='white', fontface="bold") +
      shadowtext::geom_shadowtext(data=speeds, aes(x=x, y=y, label=label), 
                                  color="red", bg.colour='white', size=3, fontface="bold") +
      scale_alpha_discrete(range=c(.3, 1)) +
      scale_x_continuous(expand=c(0,0), breaks=brk) +
      scale_y_continuous(expand=c(0,0), breaks=brk, limits=c(-3, NA)) +
      coord_fixed() +
      theme(axis.title=eb, axis.text=eb, axis.ticks=eb,
            panel.background = eb, strip.text=eb, strip.background = eb,
            legend.position="none")



# final composite figure ##################################################

library(gridExtra)
p <- arrangeGrob(cartoons, legend, nrow=1, widths=c(3, 1.2))
p <- arrangeGrob(p, map, ncol=1, heights=c(1, 2))
ggsave("figures/tailwinds/speed_isotropy_direction.png", p, width=12, height=9, units="in")

source("E:/edges/range-edges/code/utilities.r")
ggs("figures/manuscript/fig_1.png", p, width=12, height=9, units="in",
    add=grid.text(letters[1:3], 
                  x=c(.02, .72, .02), 
                  y=c(.965, .965, .62),
                  gp=gpar(fontsize=25, fontface="bold", col="black")))
