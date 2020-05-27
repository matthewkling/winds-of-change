
library(tidyverse)
library(raster)
library(dismo)

source("code/05_windshed_functions.R")


# custom windroses
if(F){
      
      intdir <- "data/windrose/monthly_p1_wnd10m"
      f <- list.files(intdir, full.names=T)
      f <- f[grepl("05_|06_|07_", f)]
      years <- 30/4
      p <- 1
      outfile <- "data/case_study/windrose_p1_1980_2009_mayjunejuly.tif"
      f <- f %>% 
            lapply(stack) %>%
            Reduce("+", .) %>%
            "/"(24 * 365 * years) %>% # number of hours
            "^"(1/p) %>%
            writeRaster(outfile, overwrite=T)
      
      
      intdir <- "data/windrose/monthly_p1_wnd10m"
      f <- list.files(intdir, full.names=T)
      f <- f[grepl("09_|10_|11_", f)]
      years <- 30/4
      p <- 1
      outfile <- "data/case_study/windrose_p1_1980_2009_septoctnov.tif"
      f <- f %>% 
            lapply(stack) %>%
            Reduce("+", .) %>%
            "/"(24 * 365 * years) %>% # number of hours
            "^"(1/p) %>%
            writeRaster(outfile, overwrite=T)
      
      
      intdir <- "data/windrose/monthly_p2_wnd10m"
      f <- list.files(intdir, full.names=T)
      f <- f[grepl("09_|10_|11_", f)]
      years <- 30/4
      p <- 2
      outfile <- "data/case_study/windrose_p2_1980_2009_septoctnov.tif"
      f <- f %>% 
            lapply(stack) %>%
            Reduce("+", .) %>%
            "/"(24 * 365 * years) %>% # number of hours
            "^"(1/p) %>%
            writeRaster(outfile, overwrite=T)
      
      
      intdir <- "data/windrose/monthly_p3_wnd10m"
      f <- list.files(intdir, full.names=T)
      f <- f[grepl("09_|10_|11_", f)]
      years <- 30/4
      p <- 3
      outfile <- "data/case_study/windrose_p3_1980_2009_septoctnov.tif"
      f <- f %>% 
            lapply(stack) %>%
            Reduce("+", .) %>%
            "/"(24 * 365 * years) %>% # number of hours
            "^"(1/p) %>%
            writeRaster(outfile, overwrite=T)
}


# species range data
spdirs <- list.dirs("F:/little_trees/raw_data", full.names=T)

sp <- "Pinus contorta"

gene <- function(sp){
      
      range <- rgdal::readOGR(spdirs[grepl(sp, spdirs)])
      ext <- extent(range)
      
      # load baseline and future temperature data
      climate <- stack("data/geographic/processed/temperature.tif")
      climate0 <- climate[[1]]
      climate <- crop(climate, ext)
      
      # rasterize species range, ensuring holes are handled properly
      rng <- climate[[1]]
      rng[] <- NA
      rrange <- fasterize::fasterize(sf::st_as_sf(range), rng)
      if(any(range$CODE==0)){
            rng0 <- fasterize::fasterize(sf::st_as_sf(range[range$CODE==0,]), rng) %>%
                  reclassify(c(NA, NA, 1, 0, 1, NA))
            rrange <- mask(rrange, rng0)
      }
      
      
      # clip climate to species range
      climate <- mask(climate, rrange)
      names(climate) <- c("clim0", "clim1")
      
      # prep wind conductance data
      land <- raster("f:/cfsr/land.tif") %>% 
            rotate()
      rose <- stack("data/case_study/windrose_p1_1980_2009_mayjunejuly.tif") %>%
            rotate()
      rose <- land %>%
            reclassify(c(NA, NA, .05)) %>% # weight of .1 makes it possible to cross narrow waterways
            "*"(rose)
      rose <- crop(rose, extent(range)*1.5)
      rose <- climate %>%
            reclassify(c(NA, NA, .25,
                         -Inf, Inf, 1)) %>% # weight of .1 for unpopulated regions
            "*"(rose)
      
      wind <- rose %>% add_coords() 
      downtrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="downwind")
      uptrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="upwind")
      
      
      # function to calculate wind and climate surfaces for a pixel
      woc <- function(x, y, downtrans, uptrans, clim, 
                      cost_to_flow, 
                      output = "summary"){
            #browser()
            coords <- c(x, y)
            origin <- matrix(coords, ncol=2)
            message(paste(coords, collapse=" "))
            
            # calculate wind catchment and climate analog surfaces
            s <- list(wind_rev = accCost(uptrans, origin),
                      clim_rev = analogs(clim, coords, sigma=2, reverse=T)) %>%
                  stack()
            
            # modify wind values
            s$wind_rev <- calc(s$wind_rev, cost_to_flow) %>% mask(clim[[1]])
            
            # overlap between wind and climate
            s$overlap_rev <- s$wind_rev * s$clim_rev
            
            # set water values to zero
            s[is.na(s[])] <- 0
            
            ss <- s %>% as.list() %>% 
                  sapply(function(x) sum(values(x)[is.finite(values(x))]))
            names(ss) <- names(s)
            return(c(x=x, y=y, ss))
      }
      
      # grid cells across species range
      pts0 <- climate %>% crop(rrange) %>% mask(rrange) %>% 
            stack(crop(land, rrange)) %>%
            rasterToPoints() %>% 
            as.data.frame()
      pts <- filter(pts0, is.finite(temperature.1))
      
      # windscapes for every cell
      cost_to_flow <- function(x) 1/x
      d <- map2(pts$x, pts$y, possibly(woc, NULL), 
                downtrans=downtrans, uptrans=uptrans, clim=climate, 
                cost_to_flow=cost_to_flow) %>%
            do.call("rbind", .) %>%
            as.data.frame()
      
      d <- select(pts0, x, y) %>% full_join(d)
      saveRDS(d, paste0("data/case_study/data_", sp, ".rds"))
      
      ld <- land %>%
            crop(extent(rrange) * 1.1) %>%
            rasterToPoints() %>%
            as.data.frame()
      
      p <- ggplot(d, aes(x, y, fill=overlap_rev)) +
            geom_raster(data=ld, aes(x, y), fill="gray90") +
            geom_raster() +
            scale_fill_gradientn(colors=c("red", "yellow", "forestgreen")) +
            theme_void() +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            coord_fixed() +
            theme(legend.position=c(.2, .3)) +
            labs(fill = "inbound\nwind-analog\noverlap")
      ggsave(paste0("data/case_study/exploratory_plots_geneflow/", sp, ".png"), 
             p, width=5, height=5, units="in")
      
      return(d)
}





##### range expansion outside range ######

prep_climate <- function(){
      
      # US & canada
      ext <- extent(-170, -52, 25, 72)
      
      template <- raster("f:/cfsr/land.tif") %>% 
            rotate() %>%
            crop(ext)
      
      # prep baseline and future climate layers
      vars <- c("prec_1_", "prec_4_", "prec_7_", "prec_10_",
                "tmin10_1_", "tmin10_4_", "tmin10_7_", "tmin10_10_",
                "tmax10_1_", "tmax10_4_", "tmax10_7_", "tmax10_10_")
      clim <- list.files("F:/chelsa/monthly48", full.names=T)
      clim <- clim[grepl(paste(vars, collapse="|"), clim)]
      clim <- stack(clim) %>%
            crop(ext) %>%
            aggregate(round(res(template)[1]/res(.)[1])) %>%
            resample(template)
      
      fmos <- c("_1_", "_10_", "_4_", "_7_")
      fvars <- c("_pr_", "_tasmax_", "_tasmin_")
      fclim <- list.files("F:/chelsa/cmip5", full.names=T)
      fclim <- fclim[grepl("2080", fclim) & grepl("rcp85", fclim) & !grepl("_tas_", fclim)]
      fclim <- fclim[grepl(paste(fmos, collapse="|"), fclim)]
      clim2 <- list()
      for(v in fvars){
            for(m in fmos){
                  message(paste(v, m))
                  fc <- fclim[grepl(v, fclim) & grepl(m, fclim)] %>%
                        stack() %>%
                        crop(ext) %>%
                        mean() %>%
                        aggregate(round(res(template)[1]/res(.)[1])) %>%
                        resample(template)
                  clim2[paste0(v, m)] <- fc
            }
      }
      fclim <- stack(clim2)
      names(fclim) <- names(clim)
      fclim <- mask(fclim, clim[[1]])
      
      saveRDS(list(clim, fclim), "data/case_study/climate/north_america.rds")
}
if(F) prep_climate()






# fit a SDM and project it to late century
# map connectivity from current range to future suitability surface
sdm <- function(sp, rosefile){
      
      climate <- readRDS("data/case_study/climate/north_america.rds")
      clim <- climate[[1]]
      fclim <- climate[[2]]
      
      # load Little range polygon
      range <- rgdal::readOGR(spdirs[grepl(sp, spdirs)])
      ext <- extent(range) * 1.5
      ext@ymax <- ext@ymax + 20
      
      # rasterize species range, ensuring holes are handled properly
      rng <- crop(clim[[1]], ext)
      rng[] <- NA
      rrange <- fasterize::fasterize(sf::st_as_sf(range), rng)
      if(any(range$CODE==0)){
            rng0 <- fasterize::fasterize(sf::st_as_sf(range[range$CODE==0,]), rng) %>%
                  reclassify(c(NA, NA, 1, 0, 1, NA))
            rrange <- mask(rrange, rng0)
      }
      
      # fit maxent model
      clim <- crop(clim, rrange)
      fclim <- crop(fclim, rrange)
      fit <- maxent(crop(clim, extent(rrange)*1.25), 
                    p=coordinates(rrange)[is.finite(values(rrange)),])
      pred0 <- predict(fit, clim)
      pred1 <- predict(fit, fclim)
      
      # prep wind conductance data
      land <- raster("f:/cfsr/land.tif") %>% rotate()
      rose <- stack(rosefile) %>% rotate()
      
      # downweight conductance over water
      rose <- land %>%
            reclassify(c(NA, NA, .05)) %>% 
            "*"(rose) %>%
            crop(rrange)
      land <- crop(land, rrange)
      
      wind <- rose %>% add_coords() 
      downtrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="downwind")
      
      # wind accessibility of future suitable area from current range
      xy <- rasterToPoints(rrange)[,1:2]
      cost <- accCost(downtrans, xy) %>%
            mask(clim[[1]])
      
      # sum of accessibility from every location in range
      ccost <- cost
      ccost[] <- 0
      for(i in 1:nrow(xy)) ccost <- ccost + 1/accCost(downtrans, xy[i,])
      ccost <- mask(ccost, clim[[1]])
      
      # geographic distance
      dst <- distance(rrange)
      
      # plot
      pn <- pmin(as.vector(extent(trim(rrange))),
                 as.vector(extent(trim(clim[[1]]))))
      px <- pmax(as.vector(extent(trim(rrange))),
                 as.vector(extent(trim(clim[[1]]))))
      plot_ext <- extent(pn[1], px[2], pn[3], px[4]) * 1.5
      
      unocc <- reclassify(rrange, c(NA, NA, 1, 0, 2, NA))
      cost <- mask(cost, unocc)
      ccost <- mask(ccost, unocc)
      dst <- mask(dst, unocc)
      
      d <- stack(clim[[1]], cost*1000, ccost, rrange, pred1, dst) %>%
            mask(land) %>%
            crop(plot_ext) %>%
            rasterToPoints() %>%
            as.data.frame()
      names(d)[3:ncol(d)] <- c("land", "wind", "windsum", "range", "suit1", "distance")
      d <- mutate(d, range = ifelse(range==0, NA, range))
      
      p <- ggplot(d, aes(x, y)) +
            geom_raster(data=filter(d, is.finite(land)), fill="gray80") +
            geom_raster(data=filter(d, is.finite(wind)), aes(fill=(wind/distance))) +
            geom_raster(data=filter(d, is.finite(range)), fill="black") +
            scale_fill_gradientn(colors=c("forestgreen", "yellow", "red"),
                                 trans="log10") +
            theme_void() +
            coord_fixed() +
            labs(title=sp) +
            theme(plot.title=element_text(hjust=.5))
      ggsave(paste0("data/case_study/exploratory_plots/", sp, ".png"), p, width=8, height=6, units="in")
      
      p <- ggplot(d, aes(x, y)) +
            geom_raster(data=filter(d, is.finite(land)), fill="gray80") +
            geom_raster(data=filter(d, is.finite(wind)), aes(fill=windsum)) +
            geom_raster(data=filter(d, is.finite(range)), fill="black") +
            scale_fill_gradientn(colors=c("red", "yellow", "forestgreen")) +
            theme_void() +
            coord_fixed() +
            labs(title=sp) +
            theme(plot.title=element_text(hjust=.5))
      ggsave(paste0("data/case_study/exploratory_plots_pressure/", sp, ".png"), p, width=8, height=6, units="in")
      
      
      return(d)
}


###########

library(gridExtra)
library(grid)

sp <- "Pinus contorta"

if(F){
      saveRDS(sdm(sp, "data/case_study/windrose_p1_1980_2009_septoctnov.tif"), "data/case_study/sdm_data_p1.rds")
      saveRDS(sdm(sp, "data/case_study/windrose_p2_1980_2009_septoctnov.tif"), "data/case_study/sdm_data_p2.rds")
      saveRDS(sdm(sp, "data/case_study/windrose_p3_1980_2009_septoctnov.tif"), "data/case_study/sdm_data_p3.rds")
      saveRDS(gene(sp), "data/case_study/geneflow_data.rds")
}

dsdm <- readRDS("data/case_study/sdm_data_p2.rds")
dgene <- readRDS("data/case_study/geneflow_data.rds")

xlims <- c(range(dsdm$x) + c(0, -7))
ylims <- c(range(dgene$y) + c(0, 3))

palette <- c("#FFB000", "#FE6100", "#DC267F", "#785EF0", "#648FFF", "#64D3FF")

p1 <- ggplot(dgene, aes(x, y, fill=overlap_rev)) +
      geom_raster(data=filter(dsdm, is.finite(land)), fill="black") +
      geom_raster() +
      annotate(geom="text", x=xlims[1]+4, y=ylims[1]+1.5, 
               size=4, lineheight=.75, hjust=0, vjust=0,
               label="genetic\nrescue") +
      scale_fill_gradientn(colors = palette, na.value = "black") +
      theme_void() +
      scale_x_continuous(expand=c(0,0), limits=xlims) +
      scale_y_continuous(expand=c(0,0), limits=ylims) +
      coord_fixed() +
      theme(legend.position="none")

p2 <- ggplot(dsdm, aes(x, y)) +
      geom_raster(data=filter(dsdm, is.finite(land)), fill="black") +
      geom_raster(data=filter(dsdm, is.finite(wind)), aes(fill=windsum*suit1)) +
      annotate(geom="text", x=xlims[1]+4, y=ylims[1]+1.5, 
               size=4, lineheight=.75, hjust=0, vjust=0,
               label="range\nexpansion") +
      scale_fill_gradientn(colors = palette) +
      scale_x_continuous(expand=c(0,0), limits=xlims) +
      scale_y_continuous(expand=c(0,0), limits=ylims) +
      theme_void() +
      coord_fixed() +
      theme(plot.title=element_text(hjust=.5),
            legend.position="none")


dx <- dsdm %>%
      select(-wind, -range, -distance) %>%
      filter(!is.na(windsum),
             x >= min(xlims), x <= max(xlims),
             y >= min(ylims), y <= max(ylims)) %>%
      mutate(windsum = log10(windsum)) %>%
      gather(stat, value, windsum, suit1) %>%
      mutate(stat = factor(stat, levels=c("suit1", "windsum"))) %>%
      group_by(stat) %>% 
      mutate(value = scales::rescale(value))

dxt <- data.frame(x=min(dx$x), y=min(dx$y),
                  stat = c("suit1", "windsum"),
                  label=c("future\nsuitability", "wind\nshadow"))

p2c <- ggplot(dx, aes(x, y)) +
      geom_raster(data=filter(dsdm, is.finite(land)), fill="black") +
      geom_raster(data=dx, aes(fill=value)) +
      geom_text(data=dxt, aes(x, y, label=label),
                nudge_x = 8, nudge_y = 6, hjust=0, lineheight=.75) +
      facet_wrap(~stat) +
      scale_fill_gradientn(colors = palette) +
      scale_x_continuous(expand=c(0,0), limits=xlims) +
      scale_y_continuous(expand=c(0,0), limits=ylims) +
      theme_void() +
      coord_fixed() +
      theme(legend.position="none",
            strip.text=element_blank())

dx1 <- dgene %>%
      select(-overlap_rev) %>%
      gather(stat, value, wind_rev, clim_rev) %>%
      group_by(stat) %>% 
      mutate(value = scales::rescale(value))


dxt1 <- data.frame(x=min(dx$x), y=min(dx$y),
                   stat = c("clim_rev", "wind_rev"),
                   label=c("analog\navailability", "wind\nshadow"))

p1c <- ggplot(dx1, aes(x, y)) +
      geom_raster(data=filter(dsdm, is.finite(land)), fill="black") +
      geom_raster(data=dx1, aes(fill=value)) +
      geom_text(data=dxt1, aes(x, y, label=label),
                nudge_x = 8, nudge_y = 6, hjust=0, lineheight=.75) +
      facet_wrap(~stat) +
      scale_fill_gradientn(colors=palette, na.value="black") +
      scale_x_continuous(expand=c(0,0), limits=xlims) +
      scale_y_continuous(expand=c(0,0), limits=ylims) +
      theme_void() +
      coord_fixed() +
      theme(legend.position="none",
            strip.text=element_blank())


pa <- arrangeGrob(p1c, p2c, nrow=1)
pb <- arrangeGrob(p1, p2, nrow=1)
p <- arrangeGrob(pb, pa, ncol=1, heights = c(2, 1))

source("E:/edges/range-edges/code/utilities.r")
ggs("figures/manuscript/fig_5.png", p, width=9, height=5.2, units="in",
    add = grid.text(letters[1:6], 
                    x=c(0, 0, .25, .5, .5, .75) + .02, 
                    y=c(.33, 0, 0, .33, 0, 0) + .06,
                    gp=gpar(fontsize=20, fontface="bold", col="black")))


