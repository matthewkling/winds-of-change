


# examples of landscapes falling in 4 quadrants:
# hindrance vs facilitation, isotropic vs anisotropic


library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(colormap)
library(grid)
library(gridExtra)

select <- dplyr::select


stop("get functions from windshed_02_analysis script")


#### find examples based on summary stats ####

d <- read_csv("data/windshed/force_500km.csv") %>%
      mutate(windfill_fwd = overlap_fwd / clim_fwd)

e <- d %>%
      filter(clim_fwd > .05, clim_fwd < 1) %>% # control for suitable area
      # mutate(hf = ifelse(log10(windfill_fwd) > -1, "facilitation",
      #                    ifelse(log10(windfill_fwd) < -3, "hindrance", NA)),
      #        ia = ifelse(abs(log10(overlap_fwd / overlap_rev)) > 1, 
      #                    "anisotropic", "isotropic")) %>%
      # filter(!is.na(hf)) %>%
      # group_by(hf, ia) %>%
      # 
      mutate(class = case_when(
            log10(windfill_fwd) < -3.75 & # low filling 
                  log10(wind_fwd) > -2.25 & log10(wind_fwd) < -2 & # moderate wind
                  abs(log10(overlap_fwd / overlap_rev)) > 1 ~ # directionality (ish)
                  "anisotropic hindrance",
            log10(windfill_fwd) < -3.75 & # low filling 
                  log10(wind_fwd) > -2.8 & log10(wind_fwd) < -2.6 & # low wind
                  abs(log10(overlap_fwd / overlap_rev)) < .25 ~ # non-directionality (ish)
                  "isotropic hindrance",
            log10(windfill_fwd) > -1.25 & # high filling 
                  log10(wind_fwd) > -2.25 & log10(wind_fwd) < -1.5 & # moderate wind
                  abs(log10(overlap_fwd / overlap_rev)) > .75 ~ # directionality (ish)
                  "anisotropic facilitation",
            log10(windfill_fwd) > -1 & # high filling 
                  log10(wind_fwd) > -1.5 & # high wind
                  abs(log10(overlap_fwd / overlap_rev)) < .5 ~ # non-directionality (ish)
                  "isotropic facilitation",
            
            TRUE ~ NA_character_)
      ) %>%
      filter(!is.na(class)) %>%
      group_by(class) %>%
      
      sample_n(10) %>%
      mutate(group_id = 1:length(wind_fwd)) %>%
      ungroup() %>%
      mutate(id = 1:nrow(.))
nrow(e)




#### get raster data ####

#transform <- function(x) .5 ^ (x * 1e9) # .05
transform <- function(x) x

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=wr, climate=climate, radius=500,
          output="rasters", transform=transform)

d <- r %>%
      lapply(function(x) x %>%
                   subset(c("wind_fwd", "clim_fwd")) %>%
                   rasterToPoints() %>%
                   as.data.frame())
for(i in 1:length(d)) d[[i]]$id <- e$id[i]
d <- d %>% bind_rows()

d <- e %>%
      mutate(px = x, py = y) %>%
      #select(hf, ia, group_id, id, px, py) %>%
      select(class, group_id, id, px, py) %>%
      left_join(d)


#### plot ####

# wind_fwd is a cost metric, invert/convert to a diffusion metric 
d <- mutate(d, wind= (1 - ecdf(wind_fwd)(wind_fwd)) ^ 2)


#d$wind <- .005 ^ (d$wind_fwd * 1000)

centers <- d %>%
      select(px, py, class, group_id, id) %>%
      distinct()

p <- ggplot(d, aes(x, y, fill=wind)) +
      geom_raster() +
      facet_wrap(~ group_id + class, ncol=4, scales="free") +
      geom_point(aes(px, py), color="red") +
      scale_fill_viridis_c() +
      theme_void() +
      theme(strip.text = element_blank(),
            plot.title=element_text(size=10)) +
      labs(title = "anisotropic facilitation -- anisotropic hidrance -- isotropic facilitation -- isotropic hindrance")



d$color <- colors2d(select(d, clim_fwd, wind),
                    c("green", "gold", "black", "cyan4"))

p <- ggplot(d, aes(x, y, fill=wind)) +
      geom_raster(fill=d$color) +
      facet_wrap(~ group_id + class, ncol=4, scales="free") +
      geom_point(aes(px, py), color="red") +
      theme_void() + # note: removing the legend or the fill aes breaks the plot
      theme(strip.text = element_blank(),
            plot.title=element_text(size=10)) +
      labs(title = "anisotropic facilitation -- anisotropic hidrance -- isotropic facilitation -- isotropic hindrance")
p
ggsave("figures/windsheds/examples.png", p, width=6, height=12, units="in")






ggplot(d, aes(x, y, 
              #fill=wind
              fill=clim_fwd
)) +
      #facet_wrap(~ group_id + hf + ia, ncol=4, scales="free") +
      facet_wrap(~ id, scales="free") +
      geom_raster() +
      geom_point(aes(px, py), color="red") +
      #geom_contour() +
      scale_fill_gradientn(colors=c("gray90", "darkblue")) +
      theme_void() +
      theme(legend.position="bottom")

pd <- d %>% filter(id == 34)
ggplot(pd, aes(x, y, 
               z=wind, fill=clim_fwd
)) +
      geom_raster() +
      stat_contour(geom='polygon', breaks=c(.1, .5, .9), 
                   color="blue", fill="blue", alpha=.25) +
      geom_point(aes(px, py)) +
      scale_fill_gradientn(colors=c("gray90", "darkred")) +
      theme_void() +
      theme(legend.position="bottom")



pd <- d %>% filter(id == 34)
ggplot(pd, aes(x, y)) +
      geom_raster(fill=pd$color) +
      geom_point(aes(px, py)) +
      theme_void()


######## make an expository plot of nice examples #######

examples <- e %>%
      filter(class=="isotropic facilitation" & group_id %in% c(4)) %>%
      select(x, y)

examples <- data.frame(x=c(9723770.1, -736629.9, 6439370, 10303370), 
                       y=c(6380287, 1725787, 4727287, 5597287),
                       class=c("anisotropic hindrance", "isotropic hindrance", 
                               "anisotropic facilitation", "isotropic facilitation"),
                       id=1:4)

r <- map2(examples$x, examples$y, possibly(woc, NULL), 
          windrose=wr, climate=climate, radius=500,
          output="rasters", transform=transform)

d <- r %>%
      lapply(function(x) x %>%
                   rasterToPoints() %>%
                   as.data.frame())
for(i in 1:length(d)) d[[i]]$id <- e$id[i]
d <- d %>% bind_rows()

d <- examples %>%
      mutate(px = x, py = y) %>%
      select(id, class, px, py) %>%
      left_join(d)

d <- mutate(d, wind= (1 - ecdf(wind_fwd)(wind_fwd)) ^ 2)

centers <- d %>%
      select(px, py, class, id) %>%
      distinct()

d$color <- colors2d(select(d, clim_fwd, wind),
                    c("green", "red", "black", "cyan3"))

landscapes <- ggplot(d, aes(x, y, fill=wind)) +
      geom_raster(fill=d$color) +
      facet_wrap(~ class, scales="free") +
      geom_point(aes(px, py), color="yellow") +
      theme_void() # note: removing the legend or the fill aes breaks the plot :(
      
scatter <- ggplot(d, aes(clim_fwd, wind)) +
      geom_point(color=d$color, size=2) +
      theme_minimal() +
      labs(x="climatic similarity",
           y="wind accessibility")

ld <- land %>% rasterToPoints() %>% as.data.frame()
map <- ggplot() +
      geom_raster(data=ld, aes(x, y), fill="gray") +
      geom_raster(data=d, aes(x, y), fill="black") +
      geom_point(data=centers, aes(px, py), color="red", size=1) +
      theme_void()

p <- arrangeGrob(scatter, map, ncol=1, heights=c(1.5, 1))
p <- arrangeGrob(landscapes, p, nrow=1, widths=c(1.5, 1))
ggsave("figures/windsheds/landscape_examples.png", p, width=9, height=6, units="in")
