



library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)
library(colormap)
library(grid)
library(gridExtra)
library(ecoclim)
library(scales)

select <- dplyr::select

#####################


f <- read_csv("data/windshed/p2_500km.csv") %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_")



### correlations among windshed stats

v <- c("centroid_bearing", "centroid_distance", 
       "windshed_bearing", "windshed_distance", 
       "windshed_isotropy", "windshed_size")
d <- f %>%
      filter(property == "wind",
             direction == "fwd") %>%
      mutate(stat = paste0(moment, "_", stat)) %>%
      filter(stat %in% v) %>%
      select(x, y, stat, value) %>%
      spread(stat, value) %>%
      sample_n(5000) %>%
      pairsData(v, c("x", "y"), mirror=T) %>%
      mutate(x_var = factor(x_var, levels=v),
             y_var = factor(y_var, levels=v))
p <- ggplot(d, aes(x_value, y_value)) +
      facet_grid(y_var ~ x_var, scales="free") +
      geom_point(size=.25) +
      labs(title = "Relationships among forward windshed summary statistics") +
      theme(axis.title=element_blank())
ggsave("figures/windsheds/windshed_stat_correlations.png", 
       width=8, height=8, units="in")



### speed-isotropy-direction space

d <- f %>%
      filter(property == "wind",
             direction == "fwd",
             moment == "windshed") %>%
      spread(stat, value)







### maps of basic windshed size properties, fwd and rev

d <- f %>%
      filter(moment == "windshed",
             stat == "size") %>%
      spread(property, value) %>%
      mutate(windfill = overlap / clim) %>%
      gather(property, value, clim, wind, overlap, windfill) %>%
      mutate(direction = ifelse(direction == "fwd", "forward", "reverse"))


for(var in unique(d$property)){
      v <- filter(d, property==var)
      trans = switch(var, 
                     clim="identity", wind="identity", 
                     overlap="log10", windfill="log10")
      map <- ggplot(v, aes(x, y, fill=value)) +
            facet_grid(direction ~ .) +
            geom_raster() +
            scale_fill_viridis_c(trans=trans) +
            theme_void() +
            theme(text=element_text(size=20),
                  strip.text=element_text(size=20, angle=-90),
                  legend.position=c(.1, .7)) +
            guides(fill=guide_colorbar(barheight=15)) +
            labs(fill = var)
      ggsave(paste0("figures/windsheds/", var, ".png"), 
             width=12, height=12, units="in")
      
      # forward vs reverse
      
      v <- v %>% 
            spread(direction, value) %>%
            mutate(color = colors2d(cbind(rank(.$forward), rank(.$reverse)),
                                    #cbind(.$forward, .$reverse),
                                    c("black", "cyan", "gray85", "magenta")))
      map <- ggplot(v, aes(x, y)) +
            geom_raster(fill = v$color) +
            theme_void()
      
      legend <- ggplot(v, aes(forward, reverse)) +
            geom_point(color=v$color, size=.1) +
            scale_y_continuous(trans=trans) +
            scale_x_continuous(trans=trans) +
            theme_minimal() +
            theme(text=element_text(size = 45)) +
            labs(x = paste("forward", var),
                 y = paste("reverse", var))
      png(paste0("figures/windsheds/", var, "_FR.png"), 
          width=3000, height=1500)
      plot(map)
      plot(legend, vp=viewport(x=.15, y=.33, width=.2, height=.4))
      dev.off()
}



## 2d map clim vs windfilling

dd <- d %>%
      select(-moment, -stat) %>%
      spread(property, value) %>%
      mutate(color = colors2d(cbind(rank(.$clim), rank(.$windfill)),
                              c("forestgreen", "yellow", "red", "black")))

txt <- data.frame(direction = c("forward", "reverse"),
                  text = c("forward\n(emigration)", "reverse\n(immigration)"),
                  x = min(dd$x), y=min(dd$y) + .3 *(diff(range(dd$y))))

map <- ggplot(dd, aes(x, y)) +
      facet_grid(direction ~ .) +
      geom_raster(fill=dd$color) +
      geom_text(data=txt, aes(x, y, label=text), 
                hjust=0, size=6, lineheight=.75) +
      theme_void() +
      #theme(strip.text=element_text(size=50, angle=-90)) +
      theme(strip.text=element_blank())
legend <- ggplot(dd, aes(clim, windfill)) +
      geom_point(color=dd$color, size=.2) +
      #scale_y_log10(breaks=c(.001, .003, .01, .03, .1, .3, 1)) +
      xlim(0, .75) +
      theme_minimal() +
      theme(#text=element_text(size = 45),
            strip.text=element_blank()) +
      labs(x = "area of analog climate",
           y = "proportion windfilling")

scat <- d %>%
      select(-moment, -stat) %>%
      filter(property %in% c("clim", "wind", "overlap")) %>%
      mutate(property = factor(property, levels=c("clim", "wind", "overlap"),
                               labels=c("area of\nanalog climate",
                                        "area of\naccessible windshed",
                                        "area of overlap\n(analog and accessible)"))) %>%
      spread(direction, value) %>%
      sample_n(30000) %>%
      ggplot(aes(forward, reverse)) +
      facet_wrap(~property, ncol=1, scales="free") +
      geom_point(size=.25, alpha=.05) +
      #coord_fixed() +
      theme_minimal() +
      labs(x="forward (emigration)",
           y="reverse (immigration)")

p <- arrangeGrob(legend, scat, ncol=1, heights=c(1, 3))
p <- arrangeGrob(map, p, ncol=2, widths=c(4, 1))
ggsave("figures/manuscript/fig_4.png", p, width=10, height=8, units="in")
ggsave("figures/windsheds/windfill_clim.png", p, width=10, height=8, units="in")

# png("figures/windsheds/windfill_clim.png", width=3000, height=2000)
# grid.draw(p)
# dev.off()





### map isotropy vs hindrance

for(drn in c("fwd", "rev")){
      
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
            left_join(d, .)
      
      d <- d %>% mutate(color = colors2d(cbind(log(.$isotropy), rank(.$windfill)),
                                         c("cyan", "magenta", "darkred", "darkblue")))
      
      map <- ggplot(d, aes(x, y)) +
            geom_raster(fill=d$color) +
            theme_void() +
            theme(strip.text=element_text(size=50, angle=-90))
      legend <- ggplot(d, aes(isotropy, windfill)) +
            geom_point(color=d$color, size=.2) +
            #scale_y_log10(breaks=c(.001, .003, .01, .03, .1, .3, 1)) +
            xlim(.25, NA) +
            theme_minimal() +
            theme(text=element_text(size = 45),
                  strip.text=element_blank()) +
            labs(x = paste(drn, "windshed isotropy"),
                 y = paste(drn, "proportion windfilling"))
      png(paste0("figures/windsheds/isotropy_windfilling_", drn, ".png"), 
          width=3000, height=1500)
      plot(map)
      plot(legend, vp=viewport(x=.15, y=.33, width=.25, height=.4))
      dev.off()
}



### climate analog curve

sigma <- 2
kernel <- function(x) exp(-.5*(x/sigma)^2)
p <- ggplot(data.frame(temp_diff = seq(-10, 10, .1)) %>%
                  mutate(similarity = kernel(temp_diff)),
            aes(temp_diff, similarity)) +
      geom_line(color="darkred", size=1) +
      theme_minimal() +
      scale_x_continuous(breaks=-10:10, limits=c(-6, 6)) +
      labs(x = "spatiotemporal temperature difference (°C)",
           y = "similarity score")
ggsave("figures/windsheds/temp_kernel.png", 
       width=6, height=4, units="in")





