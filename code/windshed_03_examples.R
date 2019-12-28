


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


stop("get functions and input datasets from windshed_02_analysis script")


f <- read_csv("data/windshed/p1_30y_250km_inv.csv") %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_")



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

#### some random examples ####

em <- data.frame(x=c(98.12386, 143.43630,  -66.25091, -38.12595, 
                     28.74896, 104.3739, -61.56342, 118.7488),
                 y=c(26.38193, -32.31396,  -2.341591, 71.34048, 
                     -12.64459, 27.31857, 2.966016, 54.16881))
em <- em[c(3, 7, 5, 8, 1, 6, 2, 4),]
#em <- em[1,]

n <- 8
e <- f %>%
      select(x, y) %>%
      filter(abs(y)<75) %>%
      sample_n(n-nrow(em)) %>%
      bind_rows(em, .) %>%
      mutate(id = 1:nrow(.))

time_conv <- identity

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=time_conv, radius=250,
          output="rasters")


d <- r[1:n] %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)) d[[i]]$id <- i
d <- bind_rows(d)
d <- filter(d, is.finite(wind_fwd))

#d <- mutate(d, wind_fwd = 1/(sqrt(wind_fwd)))
d <- mutate(d, wind_fwd = ifelse(wind_fwd > quantile(wind_fwd, .95),
                                 quantile(wind_fwd, .95),
                                 wind_fwd))
dd <- e
m <- 1

p <- ggplot() +
      geom_raster(data=d, aes(x, y, fill=wind_fwd)) +
      facet_wrap(~ id, nrow=2, scales="free") +
      geom_contour(data=d, aes(x, y, z=wind_fwd), binwidth=max(d$wind_fwd)/50,
                   color="#290000", alpha=.5, size=1) +
      geom_vline(data=e[1:n,], aes(xintercept=x), color="white", size=.5) +
      geom_hline(data=e[1:n,], aes(yintercept=y), color="white", size=.5) +
      scale_fill_gradientn(colors=c("cyan", "purple", "darkred", "#290000")) +
      guides(fill=guide_colorbar(barwidth=25)) +
      theme_void() + 
      theme(strip.text = element_blank(),
            plot.background=element_rect(fill="black", color="black"),
            text = element_text(size=25, color="white"),
            plot.title=element_text(size=10),
            legend.position="top") +
      labs(fill="wind hours from origin")
ggsave("figures/windsheds/examples/examples_wind.png", p,
       width=16, height=9, units="in")




p <- ggplot() +
      geom_raster(data=d %>% mutate(clim_fwd=ifelse(clim_fwd != 0, clim_fwd, NA)), 
                  aes(x, y, fill=clim_fwd)) +
      facet_wrap(~ id, nrow=4, scales="free") +
      scale_fill_gradientn(colors=c("yellow", "darkred", "black") %>%
                                 rev(),
                           na.value="lightblue", limits=0:1) +
      geom_point(data=e[1:n,], aes(x, y), color="green", size=2.5) +
      guides(fill=guide_colorbar(barwidth=25)) +
      theme_void() + 
      theme(strip.text = element_blank(),
            text = element_text(size=25),
            plot.title=element_text(size=10),
            legend.position="top") +
      labs(fill="climatic similarity")
ggsave("figures/windsheds/examples/examples_analogs.png", p,
       width=16, height=9, units="in")


p <- ggplot() +
      geom_raster(data=d, aes(x, y, fill=-log10(wind_fwd) * (1-clim_fwd) + 4)) +
      facet_wrap(~ id, nrow=4, scales="free") +
      scale_fill_gradientn(colors=c("black", "darkred", "yellow"),
                           na.value="lightblue") +
      #geom_point(data=e[1:n,], aes(x, y), color="cyan") +
      guides(fill=guide_colorbar(barwidth=25)) +
      theme_void() + 
      theme(strip.text = element_blank(),
            text = element_text(size=25),
            plot.title=element_text(size=10),
            legend.position="top")
ggsave("figures/windsheds/examples/examples_overlap.png", p,
       width=16, height=9, units="in")


d$color <- colors2d(select(d, clim_fwd, wind_fwd),
                    c("green", "gold", "black", "cyan"))

maps <- ggplot(d, aes(x, y)) +
      geom_raster(fill=d$color) +
      facet_wrap(~ id, ncol=8, scales="free") +
      geom_point(data=e, aes(x, y), color="red") +
      theme_void() + 
      theme(strip.text = element_blank(),
            plot.title=element_text(size=10))

scatter <- ggplot(d, aes(clim_fwd, wind_fwd)) +
      geom_point(color=d$color) +
      theme_minimal() +
      theme(plot.background = element_rect(fill="white", color="white")) +
      labs(x = "climatic analogy",
           y = "wind connectivity")
png("figures/windsheds/examples/examples.png", width=1000, height=750)
plot(maps)
plot(scatter, vp=viewport(x=0, y=1, width=.25, height=.333,
                          just = c("left", "top")))
dev.off()



### a single example of different landscape fields

ex <- filter(d, id==2) %>%
      select(-color, -id) %>%
      gather(stat, value, -x, -y) %>%
      separate(stat, c("property", "direction")) %>%
      filter(property != "overlap") %>%
      group_by(property, direction) %>%
      mutate(value = scales::rescale(value)) %>%
      spread(property, value) %>%
      mutate(overlap = wind * clim) %>%
      gather(property, value, wind, clim, overlap) %>%
      ungroup() %>%
      mutate(property = factor(property,
                               levels = c("wind", "clim", "overlap"),
                               labels = c("wind accessibility", "climate analogy",
                                          "wind-climate overlap")),
             direction = factor(direction,
                                levels=c("fwd", "rev"),
                                labels=c("forward (emigration)", "reverse (immigration)")))
saveRDS(ex, "figures/windsheds/ex_dat.rds")

center <- ex %>%
      filter(property == "wind accessibility", 
             direction == "forward (emigration)") %>%
      filter(value == max(value)) %>%
      select(x, y)

p <- ggplot() +
      facet_grid(direction ~ property) +
      geom_raster(data=ex, aes(x, y, fill=value)) +
      geom_point(data=center, aes(x, y), color="red", size=1) +
      scale_fill_gradientn(colors=c("gray90", "black")) +
      theme_void() +
      theme(strip.text.y = element_text(angle=-90),
            legend.position="none") +
      coord_fixed()
ggsave("figures/windsheds/examples/example_landscape.png", 
       p, width=6, height=4, units="in")




#### examples based on summary stats ####

# classify by isotropy vs anisotropy, hindrance vs facilitation
# and grab some prospective examples of each type

d <- f %>%
      filter(direction == "fwd",
             moment == "windshed",
             stat == "size") %>%
      spread(property, value) %>%
      mutate(windfill = overlap / clim)
d <- f %>%
      filter(direction == "fwd",
             moment == "windshed",
             property == "wind",
             stat == "isotropy") %>%
      rename(isotropy = value) %>%
      select(x, y, isotropy) %>%
      left_join(d, .)

e <- d %>%
      filter(clim > .05, clim < .1) %>% # control for suitable area
      mutate(hf = case_when(windfill < .1 ~ "hindrance",
                            windfill > .4 ~ "facilitation",
                            TRUE ~ NA_character_),
             ia = case_when(isotropy > .9 ~ "isotropic",
                            isotropy > .5 ~ "anisotropic",
                            TRUE ~ NA_character_)) %>%
      na.omit() %>%
      arrange(ia, hf) %>%
      group_by(ia, hf) %>%
      sample_n(10) %>%
      ungroup() %>%
      mutate(id = 1:nrow(.))

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=time_conv, radius=500,
          output="rasters")


d <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)) d[[i]]$id <- e$id[i]
d <- bind_rows(d) %>%
   filter(is.finite(clim_fwd), is.finite(wind_fwd))

d$color <- colors2d(select(d, clim_fwd, wind_fwd),
                    c("green", "gold", "black", "cyan"))

maps <- ggplot(d, aes(x, y)) +
      geom_raster(fill=d$color) +
      facet_wrap(~ id, nrow=4, scales="free") +
      theme_void()
ggsave("figures/windsheds/archetype_candidates_v2.png", maps, width=10, height=5, units="in")


# identify archetypes, get data, plot

#archetype_af <- e[5,] # amazing example
#archetype_ah <- e[18,] # ok example
#archetype_ih <- e[40,] # great example
#archetype_if <- e[27,] # ok example

archetypes <- bind_rows(archetype_af, archetype_ah, archetype_if, archetype_ih) %>%
      mutate(id = 1:nrow(.))
saveRDS(archetypes, "figures/windsheds/archetypes.rds")

r <- map2(archetypes$x, archetypes$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, cost_to_flow=cost_to_flow, radius=500,
          output="rasters")

d <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)) d[[i]]$id <- archetypes$id[i]
d <- bind_rows(d)

d <- archetypes %>%
      rename(cx = x, cy = y) %>%
      left_join(d)

d$color <- colors2d(select(d, clim_fwd, wind_fwd),
                    c("green", "gold", "black", "cyan"))

landscapes <- ggplot(d, aes(x, y)) +
      geom_raster(fill=d$color) +
      facet_wrap(~ ia + hf, scales="free", nrow=2) +
      geom_point(aes(cx, cy), color="red") +
      theme_void() +
      theme(strip.text=element_text(size=10))

scatter <- ggplot(d, aes(clim_fwd, wind_fwd)) +
      geom_point(color=d$color, size=2) +
      theme_minimal() +
      labs(x="climatic similarity",
           y="wind accessibility")

ld <- land %>% rasterToPoints() %>% as.data.frame()
map <- ggplot() +
      geom_raster(data=ld, aes(x, y), fill="gray") +
      geom_raster(data=d, aes(x, y), fill="black") +
      geom_point(data=d, aes(cx, cy), color="red", size=1) +
      theme_void()

p <- arrangeGrob(scatter, map, ncol=1, heights=c(1.5, 1))
p <- arrangeGrob(landscapes, p, nrow=1, widths=c(1.5, 1))
ggsave("figures/windsheds/archetypes.png", p, width=9, height=6, units="in")

