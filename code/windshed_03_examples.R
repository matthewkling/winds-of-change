


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

f <- read_csv("data/windshed/force_500km_v2.csv") %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_")



   

#### some random examples ####

e <- f %>%
      select(x, y) %>%
      sample_n(48) %>%
      mutate(id = 1:nrow(.))

cost_to_flow <- function(cost) (1/(cost)) ^ (1/3)

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, cost_to_flow=cost_to_flow, radius=500,
          output="rasters")

d <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)) d[[i]]$id <- i
d <- bind_rows(d)

d <- filter(d, is.finite(wind_fwd))

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
png("figures/windsheds/examples.png", width=1000, height=750)
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
ggsave("figures/windsheds/example_landscape.png", 
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
          windrose=rose, climate=climate, cost_to_flow=cost_to_flow, radius=500,
          output="rasters")

d <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)) d[[i]]$id <- e$id[i]
d <- bind_rows(d)

d$color <- colors2d(select(d, clim_fwd, wind_fwd),
                    c("green", "gold", "black", "cyan"))

maps <- ggplot(d, aes(x, y)) +
      geom_raster(fill=d$color) +
      facet_wrap(~ id, nrow=4, scales="free") +
      theme_void()
ggsave("figures/windsheds/archetype_candidates.png", maps, width=10, height=5, units="in")


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

