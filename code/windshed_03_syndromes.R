

# examples of landscapes falling in 4 quadrants:
# hindrance vs facilitation, isotropic vs anisotropic


library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(colormap)
library(grid)
library(gridExtra)
library(ecoclim)

select <- dplyr::select




################################# windshed model setup #####################################

# adapted from windshed_03_examples.r

stop("get functions and input datasets from windshed_02_analysis script")


# f0 <- read_csv("data/windshed/p1_30y_250km_inv.csv") %>%
#       select(-runtime) %>%
#       gather(var, value, -x, -y) %>%
#       separate(var, c("property", "direction", "moment", "stat"), sep="_")



woc <- function(x, y, windrose, climate, 
                radius = 250, time_conv=identity,
                sigma = 2,
                output = "summary"){
      #browser()
      start <- Sys.time()
      coords <- c(x, y)
      origin <- matrix(coords, ncol=2)
      message(paste(coords, collapse=" "))
      
      # constrain analysis to region around focal point
      coords_ll <- SpatialPoints(origin, crs(climate)) %>%
            spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
            coordinates()
      circle <- geo_circle(coords_ll, width = radius * 1000) # km to m
      
      # prepare wind and climate datasets       
      wind <- windrose %>% crop(circle) %>% mask(circle) %>% add_coords() 
      downtrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="downwind")
      uptrans <- wind %>% transition_stack(windflow, directions=8, symm=F, direction="upwind")
      clim <- climate %>% crop(circle) %>% mask(circle)
      
      # calculate wind catchment and climate analog surfaces
      s <- list(wind_fwd = accCost(downtrans, origin) / 3600,
                wind_rev = accCost(uptrans, origin) / 3600,
                clim_fwd = analogs(clim, coords, sigma=sigma),
                clim_rev = analogs(clim, coords, sigma=sigma, reverse = T)) %>%
            stack()
      
      # modify wind values
      s$wind_fwd <- calc(s$wind_fwd, time_conv) %>% mask(clim[[1]])
      s$wind_rev <- calc(s$wind_rev, time_conv) %>% mask(clim[[1]])
      
      # overlap between wind and climate
      s$overlap_fwd <- s$wind_fwd * s$clim_fwd
      s$overlap_rev <- s$wind_rev * s$clim_rev
      
      # set water values to zero
      #s$wind_fwd[is.na(s$wind_fwd[])] <- max(na.omit(values(s$wind_fwd)))
      #s[is.na(s[])] <- 0
      s <- mask(s, circle)
      if(output == "rasters") return(s)
      
      # summary statistics of wind and climate surfaces
      ex <- extent(s)
      if(ex@xmin < -180){
            shft <- -180 - ex@xmin
            s <- shift(s, shft) # ws_summarize needs real geography
            origin[1,1] <- origin[1,1] + shft
      }
      if(ex@xmax > 180){
            shft <- 180 - ex@xmax
            s <- shift(s, shft)
            origin[1,1] <- origin[1,1] + shft
      }
      
      ss <- s %>% as.list() %>% lapply(ws_summarize, origin=origin)
      for(i in 1:length(ss)) names(ss[[i]]) <- paste0(names(s)[i], "_", names(ss[[i]]))
      ss <- unlist(ss)
      
      if(output == "summary") return(c(x=x, y=y, ss,
                                       runtime = difftime(Sys.time(), start, units="secs")))
      stop("invalid output argument")
}


# load current/future climate data
climate <- stack("data/geographic/processed/temperature.tif") %>% unwrap(180)

# load land data
land <- raster("f:/cfsr/land.tif") %>% 
      rotate() %>% unwrap(180)

# load windrose data
rose <- stack("data/windrose/windrose_p1_wnd10m.tif") %>%
      rotate() %>% unwrap(180)

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, water)) %>% # weight=.1 makes it possible to cross narrow waterways
      "*"(rose)





################################# classify syndromes #####################################

# adapted from windshed_04_figures.r


truncate <- function(x, q=.005, sides=c("high")){
      q <- quantile(x, c(q, 1-q), na.rm=T)
      if("low" %in% sides) x[x<q[1]] <- q[1]
      if("high" %in% sides) x[x>q[2]] <- q[2]
      x
}


infile <- "data/windshed/p1_30y_250km_inv.csv"

f <- read_csv(infile) %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_") %>%
      mutate(direction = ifelse(direction=="fwd", "outbound", "inbound"))


drn <- "outbound"

dg <- function(lat) distGeo(c(0, lat), c(1, lat))/1000
areas <- f %>% select(y) %>% distinct()
areas$area <- sapply(areas$y, dg)

bearings <- f %>%
      filter(direction == drn,
             moment == "centroid",
             stat == "bearing") %>%
      spread(property, value) %>%
      select(x, y, clim, wind) %>%
      rename(clim_bearing = clim, wind_bearing = wind) %>%
      mutate(divergence = abs((clim_bearing - wind_bearing + 540) %% 360 - 180))

d <- f %>%
      filter(direction == drn,
             moment == "windshed",
             stat == "size") %>%
      spread(property, value) %>%
      mutate(windfill = overlap / clim)
d <- f %>%
      filter(direction == drn,
             moment == "windshed",
             property == "wind",
             stat == "isotropy") %>%
      rename(isotropy = value) %>%
      select(x, y, isotropy) %>%
      left_join(d, .) %>%
      
      mutate(cell_id = 1:nrow(.)) %>%
      
      filter(is.finite(windfill),
             is.finite(clim),
             is.finite(isotropy)) %>%
      mutate(isotropy=truncate(isotropy),
             windfill = truncate(windfill, .001, c("high", "low"))) %>%
      
      left_join(bearings) %>%
      mutate(andiv = (1 - isotropy) * divergence) %>%
      
      left_join(areas) %>%
      arrange(clim) %>%
      mutate(climrnk = cumsum(area) / sum(area)) %>%
      arrange(windfill) %>%
      mutate(windfillrnk = cumsum(area) / sum(area)) %>%
      arrange(overlap) %>%
      mutate(overlaprnk = cumsum(area) / sum(area)) %>%
      arrange(isotropy) %>%
      mutate(isotropyrnk = cumsum(area) / sum(area)) %>%
      
      arrange(andiv) %>%
      mutate(andivrnk = cumsum(area) / sum(area)) %>%
      
      # mutate(syndrome = case_when(climrnk < mean(climrnk) ~ "climate-limited",
      #                             windfillrnk > mean(windfillrnk) ~ "wind-facilitated",
      #                             isotropyrnk > mean(isotropyrnk) ~ "speed-hindered",
      #                             TRUE ~ "direction-hindered"))
      
      # mutate(syndrome = case_when(climrnk < .25 ~ "climate-limited",
      #                             #windfillrnk > mean(windfillrnk) ~ "wind-facilitated",
      #                             overlaprnk > .66 ~ "wind-facilitated",
      #                             andivrnk > .5 ~ "direction-hindered",
      #                             TRUE ~ "speed-hindered"))

mutate(syndrome = ifelse(climrnk < .25,  "climate-limited", "other")) %>%
      group_by(syndrome) %>%
      arrange(syndrome, windfill) %>% ##### use overlap or windfill here #####
      mutate(overlaprnk = cumsum(area) / sum(area)) %>%
      ungroup() %>%
      mutate(syndrome = case_when(syndrome != "other" ~ syndrome,
                                  overlaprnk > .66 ~ "wind-facilitated",
                                  TRUE ~ "other")) %>%
      group_by(syndrome) %>%
      arrange(syndrome, andiv) %>%
      mutate(andivrnk = cumsum(area) / sum(area)) %>%
      ungroup() %>%
      mutate(syndrome = case_when(syndrome != "other" ~ syndrome,
                                  andivrnk > .5 ~ "direction-hindered",
                                  TRUE ~ "speed-hindered"))


d %>% sample_n(5000) %>% 
      ggplot(aes(clim, windfill, color=syndrome)) + 
      geom_point() +
      scale_color_manual(values=c("white", "red", "gold", "black"))

d %>% group_by(syndrome) %>% summarize(area = sum(area)) %>% mutate(area = area/sum(area))



################################# maps for syndrome exemplars ######################################

inverse <- function(x) 1/x

e <- d %>%
      na.omit() %>%
      filter(abs(y) < 60) %>%
      filter(
            (syndrome != "climate-limited" | climrnk < quantile(climrnk, .1)),
            #(syndrome != "wind-facilitated" | windfillrnk > quantile(windfillrnk, .8)),
            (syndrome != "wind-facilitated" | overlaprnk > quantile(overlaprnk, .9)),
            (syndrome != "direction-hindered" | andivrnk > quantile(andivrnk, .9)),
            (syndrome != "speed-hindered" | wind < .01)
      ) %>%
      filter(syndrome == "climate-limited" | 
                   between(clim, .35, .45)) %>% # control for climate or facilitation ratio hard to intuit
      group_by(syndrome) %>%
      sample_n(10) %>%
      ungroup()

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=inverse, radius=250,
          output="rasters")


dd <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(dd)) dd[[i]]$cell_id <- e$cell_id[i]
dd <- bind_rows(dd) %>%
      filter(is.finite(clim_fwd), is.finite(wind_fwd)) %>%
      left_join(select(e, cell_id, syndrome))

dd$color <- colors2d(cbind((dd$clim_fwd), rank(dd$wind_fwd)^4),
                     c("green", "yellow", "black", "cyan"))

maps <- ggplot(dd, aes(x, y)) +
      geom_point(data=e, color="red", size=10) +
      geom_raster(fill=dd$color) +
      facet_wrap(~ syndrome + cell_id, nrow=4, scales="free") +
      theme_void()
ggsave("figures/windsheds/archetype_candidates_v2.png", maps, width=10, height=5, units="in")



###########

e <- filter(d, cell_id %in% c(148858, 104183, 111686, 178733))

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=inverse, radius=250,
          output="rasters")

dd <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(dd)) dd[[i]]$cell_id <- e$cell_id[i]
dd <- bind_rows(dd) %>%
      filter(is.finite(clim_fwd), is.finite(wind_fwd)) %>%
      left_join(select(e, cell_id, syndrome))

dd$color <- colors2d(cbind((dd$clim_fwd), rank(dd$wind_fwd)^4),
                     c("green", "yellow", "black", "cyan"))

maps <- ggplot(dd, aes(x, y)) +
      geom_point(data=e, color="red", size=10) +
      geom_raster(fill=dd$color) +
      facet_wrap(~ syndrome, nrow=1, scales="free") +
      theme_void()
ggsave("figures/windsheds/syndrome_archetypes.png", maps, width=8, height=2, units="in")

