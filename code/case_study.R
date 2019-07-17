

spdirs <- list.dirs("F:/little_trees/raw_data", full.names=T)


sp <- "Pinus contorta"
range <- rgdal::readOGR(spdirs[grepl(sp, spdirs)])



##### gene flow within species range ######

# run the first parts of windshed_02_analysis

plot(crop(climate[[1]], range), col="black")
plot(crop(climate[[1]], range) %>% mask(range), add=T)



cost_to_flow <- function(cost) (1/cost) ^ .5



climate <- stack("data/geographic/processed/temperature.tif")
climate <- crop(climate, extent(range)*1.5)


rng <- climate
rng[] <- NA
rng <- fasterize::fasterize(sf::st_as_sf(range), rng)



climate <- mask(climate, range)
names(climate) <- c("clim0", "clim1")





land <- raster("f:/cfsr/land.tif") %>% 
      rotate()
rose <- stack("data/windrose/windrose_p2_2000s.tif") %>%
      rotate()
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight of .1 makes it possible to cross narrow waterways
      "*"(rose)
rose <- crop(rose, extent(range)*1.5)
rose <- climate %>%
      reclassify(c(NA, NA, .1,
                   -Inf, Inf, 1)) %>% # weight of .1 for unpopulated regions
      "*"(rose)


# prepare wind and climate datasets       
wind <- rose %>% add_coords() 
downtrans <- wind %>% transition_stack(downwind, directions=8, symm=F)
uptrans <- wind %>% transition_stack(upwind, directions=8, symm=F)


woc <- function(x, y, downtrans, uptrans, clim, 
                radius = 1000, cost_to_flow, 
                output = "summary"){
      #browser()
      coords <- c(x, y)
      origin <- matrix(coords, ncol=2)
      message(paste(coords, collapse=" "))
      
      # calculate wind catchment and climate analog surfaces
      s <- list(wind_fwd = accCost(downtrans, origin),
                wind_rev = accCost(uptrans, origin),
                clim_fwd = analogs(clim, coords, sigma=2), # sigma=2 inispired by Rehfeldt 1999
                clim_rev = analogs(clim, coords, sigma=2, reverse=T)) %>%
            stack()
      
      # modify wind values
      s$wind_fwd <- calc(s$wind_fwd, cost_to_flow) %>% mask(clim[[1]])
      s$wind_rev <- calc(s$wind_rev, cost_to_flow) %>% mask(clim[[1]])
      
      # overlap between wind and climate
      s$overlap_fwd <- s$wind_fwd * s$clim_fwd
      s$overlap_rev <- s$wind_rev * s$clim_rev
      
      # set water values to zero
      s[is.na(s[])] <- 0
      
      # summary statistics of wind and climate surfaces
      ss <- s %>% as.list() %>% lapply(ws_summarize, origin=origin)
      for(i in 1:length(ss)) names(ss[[i]]) <- paste0(names(s)[i], "_", names(ss[[i]]))
      ss <- unlist(ss)
      return(c(x=x, y=y, ss))
}


pts <- climate %>% crop(range) %>% mask(range) %>% 
      rasterToPoints() %>% 
      as.data.frame()

d <- map2(pts$x[1:5], pts$y[1:5], possibly(woc, NULL), 
     downtrans=downtrans, uptrans=uptrans, clim=climate, 
     cost_to_flow=cost_to_flow, radius=500) %>%
      do.call("rbind", .) %>%
      as.data.frame()

saveRDS(d, "data/case_study/data.rds")



f <- read_csv("data/case_study/data.rds") %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_")



ggplot(d, aes(x, y, fill=overlap_fwd)) +
      geom_raster() +
      scale_fill_gradientn(colors=c("lavender", "#6b49ff", 
                                    "dodgerblue", "turquoise", "#48c13f", 
                                    "gold", "red", "#660000")) +
      theme_void()
      


ggplot(pts %>% gather(time, value, clim0, clim1), 
       aes(x, y, fill=value)) +
      facet_grid(.~time) +
      geom_raster() +
      scale_fill_gradientn(colors=c("lavender", "#6b49ff", 
                                    "dodgerblue", "turquoise", "#48c13f", 
                                    "gold", "red", "#660000")) +
      theme_void()


# figure out how to make the range holes actual holes


