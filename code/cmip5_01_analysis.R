
library(ecoclim)
library(raster)
library(tidyverse)
library(rNOMADS)

map <- purrr::map

####### calculate long-term means #######

# derive three sets of ensemble means, each using only the models that are available
# for all the eras being compared: historic-future, lgm-historic, lgm-historic-future.


### derive climatology per model-run

average <- function(x, # a data frame with "path", and "set" variables
                        outdir){
      #x <- f[[18]] 
      message(x$set[1])
      outfile <- paste0(outdir, "/", x$set[1], ".tif")
      if(file.exists(outfile)) return("skipping")
      
      # pull out the right years, and average them
      s <- map(x$path, stack)
      if(x$era[1] != "lgm"){
            s <- map(s, function(x){
                  y <- as.integer(substr(names(x), 2, 5)) %in% c(1901:2000, 2071:2100) %>% which()
                  if(length(y) > 0) x <- x[[y]]
            })
      }
      s <- s[!sapply(s, is.null)] %>%
            stack() %>% 
            mean() %>%
            writeRaster(outfile)
      return("averaged")
}

ff <- parseMetadata("E:/cmip5_wind",
                   keys=list(era = c("historical", "rcp85", "lgm"),
                             var = c("uas", "vas"),
                             mod = basename(list.dirs("E:/cmip5_wind/historical/Amon/uas", recursive=F)))) %>%
      mutate(run = basename(dirname(path)))

# future
f <- ff %>%
      filter(era != "lgm") %>%
      group_by(var, mod, run) %>%
      mutate(n = length(unique(era))) %>%
      filter(n == 2) %>% # models with historic and future data
      mutate(set = paste(era, var, mod, run, sep="_")) %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/historic_future/model-runs/")

# lgm
f <- ff %>%
      filter(era != "rcp85") %>%
      group_by(var, mod, run) %>%
      mutate(n = length(unique(era))) %>%
      filter(n == 2) %>% # models with lgm and historic data
      mutate(set = paste(era, var, mod, run, sep="_")) %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/lgm_historic/model-runs/")

# all
f <- ff %>%
      group_by(var, mod, run) %>%
      mutate(n = length(unique(era))) %>%
      filter(n == 3) %>% # models with data for all three eras
      mutate(set = paste(era, var, mod, run, sep="_")) %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/lgm_historic_future/model-runs/")



### average per model, across runs

average <- function(x, # a data frame with "path", and "set" variables
                    outdir){
      message(x$set[1])
      outfile <- paste0(outdir, "/", x$set[1], ".tif")
      if(file.exists(outfile)) return("skipping")
      s <- map(x$path, stack) %>%
            stack() %>% 
            mean() %>%
            writeRaster(outfile)
      return("averaged")
}

# future
f <- data.frame(path = list.files("E:/wind/cmip5_climatologies/historic_future/model-runs/", full.names=T),
                stringsAsFactors = F) %>%
      mutate(file = sub("\\.tif", "", basename(path))) %>%
      separate(file, c("era", "var", "mod", "run"), sep="_") %>%
      mutate(set = paste(era, var, mod, sep="_")) %>%
      
      group_by(var, mod, run) %>%
      mutate(n = length(unique(era))) %>%
      filter(n == 2) %>% ungroup() %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/historic_future/models/")

# lgm
f <- data.frame(path = list.files("E:/wind/cmip5_climatologies/lgm_historic/model-runs/", full.names=T),
                stringsAsFactors = F) %>%
      mutate(file = sub("\\.tif", "", basename(path))) %>%
      separate(file, c("era", "var", "mod", "run"), sep="_") %>%
      mutate(set = paste(era, var, mod, sep="_")) %>%
      
      group_by(var, mod, run) %>%
      mutate(n = length(unique(era))) %>%
      filter(n == 2) %>% ungroup() %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/lgm_historic/models/")

# all
f <- data.frame(path = list.files("E:/wind/cmip5_climatologies/lgm_historic_future/model-runs/", full.names=T),
                stringsAsFactors = F) %>%
      mutate(file = sub("\\.tif", "", basename(path))) %>%
      separate(file, c("era", "var", "mod", "run"), sep="_") %>%
      mutate(set = paste(era, var, mod, sep="_")) %>%
      
      group_by(var, mod, run) %>%
      mutate(n = length(unique(era))) %>%
      filter(n == 3) %>% ungroup() %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/lgm_historic_future/models/")





### average across all models

average <- function(x, # a data frame with "path", and "set" variables
                    outdir){
      message(x$set[1])
      outfile <- paste0(outdir, "/", x$set[1], ".tif")
      if(file.exists(outfile)) return("skipping")
      
      pd <- pts <- expand.grid(x=0:360, y=-90:90)
      coordinates(pts) <- c("x", "y")
      
      pd$w <- map(x$path, stack) %>% 
            sapply(raster::extract, y=pts) %>%
            apply(1, mean)
      
      w <- rasterFromXYZ(pd) %>%
            writeRaster(outfile)
      return("averaged")
}

# future
f <- data.frame(path = list.files("E:/wind/cmip5_climatologies/historic_future/models/", full.names=T),
                stringsAsFactors = F) %>%
      mutate(file = sub("\\.tif", "", basename(path))) %>%
      separate(file, c("era", "var", "mod"), sep="_") %>%
      mutate(set = paste(era, var, sep="_")) %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/historic_future/ensemble/")

# lgm
f <- data.frame(path = list.files("E:/wind/cmip5_climatologies/lgm_historic/models/", full.names=T),
                stringsAsFactors = F) %>%
      mutate(file = sub("\\.tif", "", basename(path))) %>%
      separate(file, c("era", "var", "mod"), sep="_") %>%
      mutate(set = paste(era, var, sep="_")) %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/lgm_historic/ensemble/")

# all
f <- data.frame(path = list.files("E:/wind/cmip5_climatologies/lgm_historic_future/models/", full.names=T),
                stringsAsFactors = F) %>%
      mutate(file = sub("\\.tif", "", basename(path))) %>%
      separate(file, c("era", "var", "mod"), sep="_") %>%
      mutate(set = paste(era, var, sep="_")) %>%
      split(.$set)
m <- purrr::map(f, possibly(average, NA), 
                outdir="E:/wind/cmip5_climatologies/lgm_historic_future/ensemble/")






##########################

# land
pd <- pts <- expand.grid(x=0:360, y=-90:90)
coordinates(pts) <- c("x", "y")
pd$land <- raster("f:/cfsr/land.tif") %>% 
      raster::extract(pts)
land <- rasterFromXYZ(pd) %>%
      rotate()

continents <- map_data("world")

cell_area <- function(lon, lat){
      geosphere::areaPolygon(data.frame(lon = lon + c(.5, .5, -.5, -.5),
                             lat = lat+ c(.5, -.5, -.5, -.5)))
}



w <- wr <- list.files("E:/wind/winds_of_change/cmip5_climatologies/lgm_historic_future/ensemble/", full.names=T) %>%
      stack() %>%
      rotate()
for(i in 1:nlayers(w)) w[[i]] <- raster::focal(w[[i]], 
                                               matrix(c(0,0,0,3,3,3,0,0,0), nrow=3, byrow=T), 
                                               mean, NAonly=T, na.rm=T)
names(w) <- names(wr)
crs(w) <- list.files("E:/wind/winds_of_change/cmip5_climatologies/lgm_historic_future/models/", full.names=T)[1] %>%
      raster() %>% crs()
w <- w %>%
      stack(land) %>%
      projectRaster(crs=CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")) %>%
      trim() %>%
      rasterToPoints() %>%
      as.data.frame()

d <- w %>%
      na.omit() %>%
      dplyr::select(-land) %>%
      gather(stat, value, historical_uas:rcp85_vas) %>%
      separate(stat, c("era", "var")) %>%
      spread(var, value) %>%
      mutate(magnitude = MagnitudeAzimuth(uas, vas)$magnitude,
             azimuth = MagnitudeAzimuth(uas, vas)$azimuth) %>%
      gather(var, value, uas:azimuth) %>%
      spread(era, value)

pd <- pairsData(d, c("lgm", "historical", "rcp85"), "var")

ggplot(d, aes(historical, lgm)) +
      facet_wrap(~var, scales="free") +
      geom_point()

pd %>% 
      filter(!(x_var=="rcp85" & y_var=="lgm")) %>%
      ggplot(aes(x_value, y_value)) +
      facet_wrap(~ var + x_var + y_var, scales="free", ncol=2) +
      geom_point()



complots <- function(variable, d){
      
      require(tidyverse)
      require(grid)
      require(gridExtra)
      select <- dplyr::select
      dv <- filter(d, var==variable)
      
      unit <- switch(variable,
                     uas="(m/s)",
                     var="(m/s)",
                     magnitude="(m/s)",
                     azimuth="(degrees)")
      
      
      #style <- theme(text=element_text(size=20, color="white"),
      #               plot.background = element_rect(fill="black", color="black"))
      style <- theme(text=element_text(size=20))
      
      bg <- annotate(geom="rect", fill="gray90", 
                     xmin=min(dv$x), xmax=max(dv$x),
                     ymin=min(dv$y), ymax=max(dv$y))
      
      # maps of individual time periods
      td <- dv %>% 
            gather(era, value, historical:rcp85) %>%
            mutate(era = factor(era, 
                                levels=c("lgm", "historical", "rcp85"),
                                labels=c("LGM", "20th century", "late 21st century")))
      col_vals <- quantile(abs(td$value), c(0, .5, 1))
      col_vals <- col_vals[c(3, 2, 1, 2, 3)] * c(-1, -1, 1, 1, 1)
      col_lims <- max(abs(td$value))*c(-1, 1)
      pal <- c("yellow", "darkred", "black", "darkblue", "cyan")
      if(variable=="azimuth"){
            col_lims <- c(0, 360)
            col_vals <- 1:20
            pal <- rainbow(20)
      }
      if(variable=="magnitude"){
            col_vals <- quantile(abs(td$value), c(0, .5, 1))
            col_lims <- c(0, max(td$value))
            pal <- c("yellow", "darkred", "black") %>% rev()
      }
      tmaps <- td %>%
            ggplot(aes(x, y, fill=value)) +
            facet_grid(.~era) +
            bg +
            geom_raster() +
            scale_fill_gradientn(colors=pal,
                                 values=scales::rescale(col_vals),
                                 limits=col_lims) +
            labs(fill=paste(variable, unit)) +
            guides(fill=guide_colorbar(barwidth=20)) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            theme_void() +
            theme(legend.position="top") +
            style
      
      # maps of deltas between periods
      dd <- dv %>%
            mutate(lgm_historical = historical - lgm,
                   historical_rcp85 = rcp85 - historical) %>%
            select(x, y, lgm_historical, historical_rcp85) %>%
            gather(stat, value, lgm_historical, historical_rcp85) %>%
            mutate(stat = factor(stat, 
                                 levels=c("lgm_historical", "historical_rcp85"),
                                 labels=c("LGM to 20th century", "20th cenury to late 21st century")))
      col_vals <- quantile(abs(dd$value), c(0, .9, 1))
      col_vals <- col_vals[c(3, 2, 1, 2, 3)] * c(-1, -1, 1, 1, 1)
      col_lims <- max(abs(td$value))*c(-1, 1)
      pal <- c("magenta", "darkorchid", "black", "darkgreen", "green")
      if(variable=="azimuth"){
            dd <- mutate(dd, 
                         value = abs(value),
                         value = pmin(value, 360-value))
            col_vals <- 1:2
            col_lims <- range(dd$value)
            pal <- c("black", "magenta")
      }
      
      dmaps <- dd %>%
            ggplot(aes(x, y, fill=value)) +
            facet_grid(.~stat) +
            bg +
            geom_raster() +
            scale_fill_gradientn(colors=pal,
                                 values=scales::rescale(col_vals),
                                 limits=col_lims) +
            labs(fill=paste("change in", variable, unit)) +
            guides(fill=guide_colorbar(barwidth=20)) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            theme_void() +
            theme(legend.position="top") +
            theme(plot.margin = unit(c(.5,0,0,0), "cm")) +
            style
      
      # scatterplots comparing time periods
      sd <- dv %>%
            mutate(lgm_historical = historical - lgm,
                   historical_rcp85 = rcp85 - historical)
      if(variable=="azimuth"){
            sd <- mutate(sd, 
                         lgm_historical = abs(lgm_historical),
                         lgm_historical = pmin(lgm_historical, 360-lgm_historical), 
                         historical_rcp85 = abs(historical_rcp85),
                         historical_rcp85 = pmin(historical_rcp85, 360-historical_rcp85))
      }
      scatter1 <- sd %>%
            ggplot(aes(lgm, historical, color=lgm_historical)) +
            geom_point() +
            scale_color_gradientn(colors=pal,
                                  values=scales::rescale(col_vals),
                                  limits=col_lims) +
            labs(x=paste0(variable, " ", unit, ", LGM"),
                 y=paste0(variable, " ", unit, ", 20th century")) +
            theme_minimal() +
            theme(legend.position = "none") +
            style
      scatter2 <- sd %>%
            ggplot(aes(historical, rcp85, color=historical_rcp85)) +
            geom_point() +
            scale_color_gradientn(colors=pal,
                                  values=scales::rescale(col_vals),
                                  limits=col_lims) +
            labs(y=paste0(variable, " ", unit, ", late 21st century"),
                 x=paste0(variable, " ", unit, ", 20th century")) +
            theme_minimal() +
            theme(legend.position = "none") +
            style
      
      p <- arrangeGrob(scatter1, scatter2, nrow=1)
      p <- arrangeGrob(tmaps, dmaps, p, ncol=1, heights=c(1.25, 1.75, 3))
      
      png(paste0("figures/cmip5/deltas_", variable, ".png"),
          width=1000, height=1000)
      grid.draw(p)
      dev.off()
      
      source("E:/edges/range-edges/code/utilities.r")
      ggs(paste0("figures/manuscript/SI_fig_cmip5_", variable, ".png"),
          p, width=15, height=15, units="in",
          add = grid.text(letters[1:7], 
                          x=c(.02, .3533, .6967,
                              .02, .52,
                              .02, .52), 
                          y=c(.91, .91, .91,
                              .68, .68,
                              .47, .47),
                          gp=gpar(fontsize=30, fontface="bold", col="black")))
}

lapply(unique(d$var), complots, d=d)



#####################################
#####################################

w <- wr <- list.files("E:/wind/cmip5_climatologies/historic_future/ensemble/", full.names=T) %>%
      stack() %>%
      rotate()
for(i in 1:nlayers(w)) w[[i]] <- raster::focal(w[[i]], 
                                               matrix(c(0,0,0,3,3,3,0,0,0), nrow=3, byrow=T), 
                                               mean, NAonly=T, na.rm=T)
names(w) <- names(wr)
w <- w %>%
      stack(land) %>%
      trim() %>%
      rasterToPoints() %>%
      as.data.frame()

d <- w %>%
      filter(is.finite(land),
             is.finite(historical_uas)) %>%
      rename(u0 = historical_uas, v0 = historical_vas,
             u1 = rcp85_uas, v1 = rcp85_vas) %>%
      mutate(m0 =  MagnitudeAzimuth(u0, v0)$magnitude,
             a0 =  MagnitudeAzimuth(u0, v0)$azimuth,
             m1 =  MagnitudeAzimuth(u1, v1)$magnitude,
             a1 =  MagnitudeAzimuth(u1, v1)$azimuth) %>%
      mutate(dm = m1 - m0,
             rm = m1/m0,
             pcm = dm / m0,
             da = abs(a1-a0),
             da = pmin(da, 360-da)) %>%
      rowwise() %>%
      mutate(area = cell_area(x, y))


p <- ggplot(d, aes(x, y, fill=rm)) +
      geom_raster() +
      scale_fill_gradient2(high="darkorchid", mid="white", low="darkgreen",
                           trans="log10", breaks=c(.1, .3, 1, 3, 10)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(-80, 90)) +
      guides(fill=guide_colorbar(ticks.colour="black")) +
      theme_void() +
      theme(legend.position=c(.1, .3),
            text=element_text(color="white")) +
      labs(fill="change in\nmagnitude")
ggsave("figures/cmip5_magnitude_change_map.png", 
       p, width=8, height=4, units="in", bg="black")


p <- ggplot(d, aes(x, y, fill=da)) +
      geom_raster() +
      scale_fill_gradient(high="darkorchid", low="white") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(-80, 90)) +
      guides(fill=guide_colorbar(ticks.colour="black")) +
      theme_void() +
      theme(legend.position=c(.1, .3),
            text=element_text(color="white")) +
      labs(fill="change in\nmagnitude")
ggsave("figures/cmip5_azimuth_change_map.png", 
       p, width=8, height=4, units="in", bg="black")


ed <- d %>% ungroup() %>%
      arrange(da) %>%
      mutate(ecdf = cumsum(area)/sum(area))

p <- ggplot(ed) +
      geom_point(aes(da, m0), size=.5) +
      geom_line(aes(da, ecdf*max(m0)), color="red") +
      scale_x_sqrt(breaks=c(0, 1, 2, 5, 10, 20, 50, 100, 180)) +
      scale_y_continuous(sec.axis = sec_axis(~./max(ed$m0), 
                                             name="cumulative proportion of land area\n")) +
      theme_minimal() +
      theme(axis.text.y.right = element_text(color = "red"),
            axis.title.y.right = element_text(color = "red"),
            axis.ticks.y.right = element_line(color = "red")) +
      labs(y = "historic wind magnitude (m/s)\n",
           x = "projected change in mean wind azimuth (degrees)")
ggsave("figures/cmip5_azimuth_change_scatter.png", p, width=6, height=6, units="in")


p <- ggplot(d, 
            aes(da, abs(dm)/m0)) +
      geom_point(size=.5) +
      geom_smooth(se=F, color="red") +
      scale_y_log10() +
      theme_minimal() +
      labs(title="relative magnitude and azimuth are changing in the same places",
           y = "absolute percent change in wind magnitude",
           x = "projected change in mean wind azimuth (degrees)")
ggsave("figures/cmip5_azimuth_change_scatter_2.png", p, width=6, height=6, units="in")

p <- ggplot(d, 
       aes(abs(dm)/m0, da)) +
      geom_point(size=.5) +
      geom_smooth(se=F, color="red") +
      scale_x_log10() +
      theme_minimal() +
      labs(title="relative magnitude and azimuth are changing in the same places",
           x = "absolute percent change in wind magnitude",
           y = "projected change in mean wind azimuth (degrees)")
ggsave("figures/cmip5_azimuth_change_scatter_2b.png", p, width=6, height=6, units="in")








ggplot(d, 
       aes(rm, da)) +
      geom_point() +
      geom_smooth(se=F, color="red") +
      scale_x_log10() +
      theme_minimal() +
      labs(title="relative magnitude and azimuth are changing in the same places")



ggplot(d, aes(abs(dm), da)) +
      geom_point() +
      labs(title="absolute magnitude and azimuth are not changing in the same places")


ggplot(d, aes(m0, da)) +
      geom_point() +
      geom_smooth(se=F, color="red") +
      theme_minimal()

ggplot(d, aes(m0, dm)) +
      geom_point() +
      geom_smooth(se=F, color="red") +
      theme_minimal()


