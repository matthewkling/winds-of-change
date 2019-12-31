

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
library(patchwork)

select <- dplyr::select




################################# windshed model setup #####################################

# adapted from windshed_03_examples.r

stop("get functions and input datasets from windshed_02_analysis script")



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


dg <- function(lat) distGeo(c(0, lat), c(1, lat))/1000
areas <- f %>% select(y) %>% distinct()
areas$area <- sapply(areas$y, dg)

df <- list()
for(drn in c("inbound", "outbound")){
   
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
      arrange(wind) %>%
      mutate(windrnk = cumsum(area) / sum(area)) %>%
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
   
   df[[drn]] <- d
}



################################# maps for syndrome exemplars ######################################

d <- df$outbound

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
ggsave("figures/windsheds/syndromes/archetype_candidates_v2.png", maps, width=10, height=5, units="in")



########### final four ###########

e <- filter(d, cell_id %in% c(148858, 104183, 111686, 178733))
#e$syndrome <- sub("-", "-\n", e$syndrome)

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


archetypes <- ggplot(dd, aes(x, y)) +
   geom_point(data=e, color="gray75", size=10) +
   geom_raster(fill=dd$color) +
   facet_wrap(~ syndrome, nrow=1, scales="free") +
   theme_void() +
   theme(strip.text=element_text(size=12))
ggsave("figures/windsheds/syndromes/syndrome_archetypes.png", 
       archetypes, width=8, height=2, units="in")








### 3d climate - anisotropy - hindrance classification

for(drn in c("inbound", "outbound")){
   
   d <- df[[drn]]
   
   palette <- c("white", "red", "gold", "black")
   
   hist <- d %>%
      mutate(y = round(y)) %>%
      group_by(y, syndrome) %>%
      summarize(n = sum(area)) %>%
      ggplot(aes(x=y, y=n, fill=syndrome)) +
      geom_histogram(stat="identity", position="fill", width=1) +
      coord_flip() +
      scale_fill_manual(values=palette) +
      #scale_y_continuous(expand=c(0,0)) +
      scale_x_continuous(expand=c(0,0), limits=c(-90, 90), 
                         breaks=seq(-90, 90, 30)) +
      theme(legend.position="none",
            panel.grid=element_blank(),
            panel.background = element_rect(fill="gray75")) +
      labs(y = "Proportion of area", x = "°N")
   
   d <- d %>% mutate(color = colors2d(cbind(.$andivrnk, .$windfillrnk),
                                      palette[c(4,3,2,4)]))
   
   whiten <- function(hex, scalar){
      #hex <- d$color[1:5]; scalar <- seq(0, 1, length.out=5)
      col2rgb(hex) %>% 
         rbind(scalar) %>%
         apply(2, function(x) 255 - ((255 - x[1:3]) * (1 - x[4]))) %>%
         t() %>% rgb(maxColorValue = 255)
   }
   d$color2 <- whiten(d$color, (1-d$climrnk)^1.5)
   
   ds <- d %>% sample_n(nrow(.))
   dummy <- data.frame(clim=max(ds$clim[ds$syndrome=="climate-limited"]),
                       windfill=min(ds$windfill[ds$syndrome=="wind-facilitated"]))
   scatter <- ds %>% 
      ggplot(aes(clim, windfill)) + 
      
      # a horrid hack to place the legend for the archetype panel
      geom_point(data=dummy, size=6, shape=15, aes(fill="climate similarity")) +
      geom_point(data=dummy, size=6, shape=15, aes(fill="wind accessibility")) +
      geom_point(data=dummy, size=6, shape=15, aes(fill="wind-climate overlap")) +
      scale_fill_manual(values=c("yellow", "green", "cyan"), 
                        guide=guide_legend(override.aes=list(shape=22, size=8))) +
      
      geom_point(color=ds$color2, size=.5) +
      annotate(geom="segment", linetype=2, 
               x=max(ds$clim[ds$syndrome=="climate-limited"]), 
               xend=max(ds$clim[ds$syndrome=="climate-limited"]), 
               y=min(ds$windfill), yend=max(ds$windfill)) +
      annotate(geom="segment", linetype=2, 
               x=max(ds$clim[ds$syndrome=="climate-limited"]), 
               xend=max(ds$clim), 
               y=min(ds$windfill[ds$syndrome=="wind-facilitated"]), 
               yend=min(ds$windfill[ds$syndrome=="wind-facilitated"])) +
      theme(legend.position=c(3, -.1), legend.direction="horizontal",
            legend.text=element_text(size=11),
            legend.key = element_rect(color="white", fill="white"),
            panel.background = element_rect(fill="gray75"),
            panel.grid=element_blank()) +
      labs(x="Analog climate availability",
           y="Wind facilitation (1/h)",
           color=NULL)
   
   map <- ggplot(d, aes(x, y, fill=syndrome)) +
      geom_point(data=d %>% filter(abs(x)<80, abs(y)<80) %>% sample_n(100), color="gray75", size=.01) +
      geom_raster(fill=d$color2) +
      scale_fill_manual(values=palette, 
                        guide=guide_legend(override.aes=list(shape=22, size=8))) +
      scale_alpha_continuous(range=0:1, guide="none") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
      theme_void() +
      theme(legend.position=c(.5, -.065), legend.direction="horizontal",
            legend.text = element_text(color="black", size=11),
            panel.background = element_rect(fill="gray75")) +
      labs(fill = NULL)
   
   
   p <- scatter + archetypes + hist + map +  
      plot_layout(nrow=2, heights=c(1.5, 3), widths=c(1, 4))
   
   
   source("E:/edges/range-edges/code/utilities.r")
   
   outfile <- paste0("figures/windsheds/syndromes/syndrome_multi_", drn, ".png")
   
   ggs(outfile, p,
       width=10, height=7, units="in",
       add = grid.text(letters[1:4], 
                       x=c(.02, .275, .02, .275), 
                       y=c(.97, .97, .58, .58),
                       gp=gpar(fontsize=20, fontface="bold", col="black")))
   
   file.copy(outfile, paste0("figures/manuscript/", basename(outfile)), overwrite = T)
   
   
   #########################################
   
   next()
   
   if(drn == "inbound"){
      map_in <- map
      hist_in <- hist
   }
   if(drn == "outbound"){
      stop()
      hist <- hist + theme(axis.text.x=element_blank(), 
                           axis.title.x=element_blank(),
                           axis.ticks.x=element_blank())
      patch <- map + hist + map_in + hist_in + 
         plot_layout(nrow=2, widths=c(4, 1), heights=c(1,1))
      ggs(paste0("figures/windsheds/syndromes/syndrome_cont_hist.png"), 
          patch, width=10, height=9, units="in",
          add = grid.text(letters[1:4], 
                          x=c(.03, .93, .03, .93), 
                          y=c(.96, .96, .49, .49),
                          gp=gpar(fontsize=25, fontface="bold", col="black")))
   }
   
   # a discrete version
   map <- ggplot(d, aes(x, y, fill=syndrome)) +
      geom_raster() +
      scale_fill_manual(values=palette) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
      theme_void() +
      theme(legend.position=c(.12, .35),
            legend.text = element_text(color="black"),
            panel.background = element_rect(fill="gray75")) +
      labs(fill = paste0(drn, "\nsyndrome"))
   ggsave(paste0("figures/windsheds/syndromes/syndrome_", drn, "2.png"), 
          width=8, height=4, units="in")
   
   
   #####
   
   
   
   ds <- d %>% arrange(clim)
   scatter_alt <- ds %>% 
      ggplot(aes(andiv, windfill)) + 
      
      # a horrid hack to place the legend for the archetype panel
      # geom_point(data=dummy, size=6, shape=15, aes(fill="climate similarity")) +
      # geom_point(data=dummy, size=6, shape=15, aes(fill="wind accessibility")) +
      # geom_point(data=dummy, size=6, shape=15, aes(fill="wind-climate overlap")) +
      scale_fill_manual(values=c("yellow", "green", "cyan"), 
                        guide=guide_legend(override.aes=list(shape=22, size=8))) +
      
      geom_point(color=ds$color2, size=.5) +
      # annotate(geom="segment", linetype=2, 
      #          x=max(ds$clim[ds$syndrome=="climate-limited"]), 
      #          xend=max(ds$clim[ds$syndrome=="climate-limited"]), 
      #          y=min(ds$windfill), yend=max(ds$windfill)) +
      # annotate(geom="segment", linetype=2, 
      #          x=max(ds$clim[ds$syndrome=="climate-limited"]), 
      #          xend=max(ds$clim), 
      #          y=min(ds$windfill[ds$syndrome=="wind-facilitated"]), 
      #          yend=min(ds$windfill[ds$syndrome=="wind-facilitated"])) +
      scale_x_sqrt() +
      theme(legend.position=c(3, -.1), legend.direction="horizontal",
            legend.text=element_text(size=11),
            legend.key = element_rect(color="white", fill="white"),
            panel.background = element_rect(fill="gray75"),
            panel.grid=element_blank()) +
      labs(x="Directional divergence",
           y="Wind facilitation (1/h)",
           color=NULL)
   
   #dr <- sample_n(d, 5000)
   #plot3d(dr$clim, dr$andiv, dr$overlap, 
   #       col=dr$color2, size=10,
   #       xlab="clim", ylab="divergence", zlab="overlap")
   
   
   
   d <- d %>% mutate(color = colors2d(cbind(.$climrnk, .$windfillrnk),
                                      c("green", "yellow", "red", "blue")))
   dr <- sample_n(d, 5000)
   plot3d(dr$clim, dr$andiv, dr$overlap, 
          col=dr$color, size=10,
          xlab="clim", ylab="divergence", zlab="overlap")
   
   
   d <- d %>% mutate(color = colors3d(cbind(.$climrnk, .$windfillrnk, .$andivrnk)))
   dr <- sample_n(d, 5000)
   plot3d(dr$clim, dr$andiv, dr$overlap, 
          col=dr$color, size=10,
          xlab="clim", ylab="divergence", zlab="overlap")
   
   d <- d %>% mutate(color = colors3d(cbind(.$climrnk, .$windrnk, .$overlaprnk)))
   dr <- sample_n(d, 5000)
   plot3d(dr$clim, dr$wind, dr$overlap, 
          col=dr$color, size=10,
          xlab="clim", ylab="divergence", zlab="overlap")
   
   
   
   
   
   d <- d %>% mutate(color = colors2d(cbind(.$andivrnk, .$windfillrnk),
                                        c("darkblue", "red", "gold", "forestgreen")))
   d$color2 <- whiten(d$color, (1-d$climrnk)^4)
   
   dr <- sample_n(d, 1000)
   plot3d(dr$clim, dr$andiv, dr$windfill, 
          col=dr$color2, size=10,
          xlab="clim", ylab="divergence", zlab="overlap")
   
   
}

