



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



infile <- "data/windshed/p1_500km.csv"



truncate <- function(x, q=.005, sides=c("high")){
      q <- quantile(x, c(q, 1-q), na.rm=T)
      if("low" %in% sides) x[x<q[1]] <- q[1]
      if("high" %in% sides) x[x>q[2]] <- q[2]
      x
}

### correlations between geography and windshed stats

# geographic and windshed attributes
xy <- read_csv(infile) %>% select(x, y)
coordinates(xy) <- c("x", "y")
geo <- list.files("data/geographic/processed", full.names=T)
geo <- stack(geo[!grepl("temperature", geo)])
geo <- extract(geo, xy) %>% cbind(latitude = coordinates(xy)[,"y"])
geo <- as.data.frame(geo)
geo <- rename(geo, continentality = coastal_distance)

f <- read_csv(infile) %>%
      select(-runtime) 

d <- geo %>% 
      select(elevation, latitude, continentality) %>%
      mutate(latitude = abs(latitude)) %>%
      cbind(f) %>%
      filter(latitude < 85) %>%
      gather(var, value, -x, -y, -elevation, -latitude, -continentality) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_") %>%
      mutate(direction = ifelse(direction=="fwd", "outbound", "inbound"))

d <- d %>%
      filter(moment == "windshed",
             stat == "size"#,
             #direction == "outbound"
      ) %>%
      gather(attr, coord, elevation, latitude, continentality)

d <- spread(d, property, value) %>%
      mutate(windfill = overlap / clim) %>%
      filter(wind > 0) %>%
      gather(property, value, clim, wind, overlap, windfill)

d <- mutate(d, property = factor(property, levels=c("clim", "wind", "overlap", "windfill"), 
                                 labels=c("climate\narea", "windshed\narea", "climate-wind\noverlap", "windfill")),
            attr = factor(attr, levels=c("latitude", "elevation", "continentality")))

p <- ggplot(d %>% sample_n(200000) %>% filter(property != "windfill"), 
            aes(coord, value)) +
      facet_grid(property ~ attr, scales="free") +
      geom_point(size=.2) +
      geom_smooth(se=F, aes(color=direction)) +
      scale_color_manual(values=c("dodgerblue", "red")) +
      theme_minimal() +
      labs(title = "Geographic relationships with forward windshed summary statistics") +
      theme(axis.title=element_blank(),
            text=element_text(size=25))
ggsave("figures/windsheds/global/scatter_windshed_geography.png", 
       width=12, height=8, units="in")


p <- ggplot(d %>% sample_n(200000) %>% filter(property != "windfill") %>% 
                  mutate(coord = case_when(attr == "elevation" ~ log10(coord),
                                           attr == "continentality" ~ log10(coord),
                                           TRUE ~ coord)), 
            aes(coord, value)) +
      facet_grid(property ~ attr, scales="free") +
      geom_point(size=.2, color="white") +
      geom_smooth(se=F, aes(color=direction)) +
      scale_color_manual(values=c("dodgerblue", "red")) +
      theme_minimal() +
      labs(title = "Geographic relationships with forward windshed summary statistics") +
      theme(axis.title=element_blank(),
            strip.text=element_text(size=25, color="white"),
            axis.text=element_text(size=15, color="white"),
            plot.background = element_rect(fill="black"))
ggsave("figures/windsheds/global/scatter_windshed_geography_bw.png", 
       width=12, height=8, units="in")


for(prp in unique(d$property)){
      
      p <- ggplot(d %>% sample_n(200000) %>% filter(property == prp) %>% 
                        mutate(coord = case_when(attr == "elevation" ~ log10(coord),
                                                 attr == "continentality" ~ log10(coord),
                                                 TRUE ~ coord),
                               coord = truncate(coord, .001, "low"),
                               value = truncate(value, .005, "high")), 
                  aes(coord, value)) +
            facet_grid(. ~ attr, scales="free") +
            geom_point(size=.2, color="white") +
            geom_smooth(se=F, aes(color=paste(direction, "     ")), size=2) +
            scale_color_manual(values=c("red", "dodgerblue")) +
            theme_minimal() +
            labs(title = "Geographic relationships with forward windshed summary statistics",
                y = sub("\\\n", " ", prp)) +
            theme(axis.title.x=element_blank(),
                  axis.title.y=element_text(color="white", size=25),
                  legend.position="bottom",
                  strip.text=element_text(size=25, color="white"),
                  legend.text=element_text(size=25, color="white"),
                  axis.text=element_text(size=15, color="white"),
                  legend.title=element_blank(),
                  plot.background = element_rect(fill="black"))
      ggsave(paste0("figures/windsheds/global/scatter_windshed_geography_bw_", 
                    sub("\\\n", "", prp), ".png"), 
             width=12, height=7.25, units="in")
}






#####################

f <- read_csv(infile) %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_") %>%
      mutate(direction = ifelse(direction=="fwd", "outbound", "inbound"))



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
ggsave("figures/windsheds/global/windshed_stat_correlations.png", 
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
      gather(property, value, clim, wind, overlap, windfill)


for(var in unique(d$property)){
      
      v <- filter(d, property==var) %>% mutate(value = truncate(value))
      
      trans <- switch(var, 
                      clim="identity", wind="identity", 
                      overlap="identity", windfill="identity")
      
      label <- switch(var, 
                      clim="analog\narea", wind="windshed\nsize", 
                      overlap="windshed-\nanalog\noverlap", 
                      windfill="analog\naccessibility")
      
      map <- ggplot(v, aes(x, y, fill=value)) +
            facet_grid(direction ~ .) +
            geom_raster() +
            #scale_fill_viridis_c(trans=trans) +
            scale_fill_gradientn(colors=c("darkblue", "red", "yellow"),
                                 trans=trans) +
            ylim(-90, 90) +
            theme_void() +
            theme(text=element_text(size=20, color="white"),
                  strip.text=element_text(size=20, angle=-90),
                  plot.background = element_rect(fill="black", color="black"),
                  legend.position=c(.1, .73)) +
            guides(fill=guide_colorbar(barheight=15)) +
            labs(fill = label)
      ggsave(paste0("figures/windsheds/global/", var, ".png"), 
             width=12, height=12, units="in")
      
      for(drn in unique(v$direction)){
            map <- ggplot(v %>% filter(direction == drn), 
                          aes(x, y, fill=value)) +
                  geom_raster() +
                  scale_fill_gradientn(colors=c("darkblue", "red", "yellow"),
                                       trans=trans) +
                  scale_x_continuous(expand=c(0,0)) +
                  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
                  theme_void() +
                  theme(text=element_text(size=20, color="white"),
                        strip.text=element_text(size=20, angle=-90),
                        plot.background = element_rect(fill="black", color="black"),
                        legend.position=c(.08, .45)) +
                  guides(fill=guide_colorbar(barheight=15)) +
                  labs(fill = label)
            ggsave(paste0("figures/windsheds/global/", var, "_", drn, ".png"), 
                   width=12, height=6, units="in")
      }
      
      # forward vs reverse
      
      v <- filter(d, property==var) %>% mutate(value = truncate(value))
      v <- v %>% 
            spread(direction, value) %>%
            mutate(color = colors2d(cbind(rank(.$outbound), rank(.$inbound)),
                                    c("white", "cyan", "gray10", "magenta")))
      #c("black", "cyan", "gray85", "magenta")))
      map <- ggplot(v, aes(x, y)) +
            geom_raster(fill = v$color) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
            theme_void() +
            theme(text=element_text(size = 45, color="white"),
                  plot.background = element_rect(fill="black", color="black"))
      
      legend <- ggplot(v, aes(outbound, inbound)) +
            geom_point(color=v$color, size=.1) +
            scale_y_continuous(trans=trans) +
            scale_x_continuous(trans=trans) +
            theme_minimal() +
            theme(text=element_text(size = 55, color="white"),
                  axis.text=element_text(color="white", size=30)) +
            labs(x = paste("outbound"),
                 y = paste("inbound"))
      png(paste0("figures/windsheds/global/", var, "_FR.png"), 
          width=3000, height=1500)
      plot(map)
      plot(legend, vp=viewport(x=.12, y=.35, width=.22, height=.44))
      dev.off()
}



## 2d map clim vs windfilling

dd <- d %>%
      select(-moment, -stat) %>%
      spread(property, value) %>%
      mutate(color = colors2d(cbind(rank(.$clim), rank(.$windfill)),
                              c("forestgreen", "yellow", "red", "black")))

txt <- data.frame(direction = c("outbound", "inbound"),
                  text = c("outbound", "inbound"),
                  x = min(dd$x), y=min(dd$y) + .3 *(diff(range(dd$y))))

map <- ggplot(dd, aes(x, y)) +
      facet_grid(direction ~ .) +
      geom_raster(fill=dd$color) +
      geom_text(data=txt, aes(x, y, label=text), 
                hjust=0, size=6, lineheight=.75) +
      ylim(-90, 90) +
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
      #sample_n(30000) %>%
      ggplot(aes(outbound, inbound)) +
      facet_wrap(~property, ncol=1, scales="free") +
      geom_point(size=.25, alpha=.05) +
      #coord_fixed() +
      theme_minimal() +
      labs(x="outbound",
           y="inbound")

p <- arrangeGrob(legend, scat, ncol=1, heights=c(1, 3))
p <- arrangeGrob(map, p, ncol=2, widths=c(4, 1))
ggsave("figures/manuscript/fig_4.png", p, width=10, height=8, units="in")
ggsave("figures/windsheds/global/windfill_clim.png", p, width=10, height=8, units="in")

# png("figures/windsheds/windfill_clim.png", width=3000, height=2000)
# grid.draw(p)
# dev.off()



### climate vs windfilling

for(drn in c("inbound", "outbound")){
      
      dd <- d %>%
            filter(direction==drn) %>%
            select(-moment, -stat) %>%
            spread(property, value) %>%
            mutate(windfill = overlap / clim) %>%
            filter(is.finite(clim), is.finite(windfill)) %>%
            mutate(clim = truncate(clim),
                   windfill = truncate(windfill)) %>%
            #mutate(color = colors2d(cbind((.$clim), (.$windfill)),
            #                        c("white", "red", "black", "forestgreen"))) %>%
            mutate(color = colorwheel2d(cbind((.$clim), (.$windfill)),
                                        c("gray",
                                          
                                          "pink", "white", "lightgreen", "green",
                                          "darkgreen", "black", "darkred", "red")))
      
      map <- ggplot(dd, aes(x, y)) +
            geom_raster(fill=dd$color) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
            theme_void() +
            theme(text=element_text(color="white"),
                  plot.background = element_rect(fill="black", color="black"))
      legend <- ggplot(dd, aes(clim, windfill)) +
            geom_point(color=dd$color, size=.2) +
            theme_minimal() +
            theme(text=element_text(size = 45, color="white"),
                  axis.text=element_text(color="white"),
                  strip.text=element_blank()) +
            labs(x = "area of analog climate",
                 y = "proportion windfilling")
      png(paste0("figures/windsheds/global/windfill_clim_", drn, ".png"), 
          width=3000, height=1500)
      plot(map)
      plot(legend, vp=viewport(x=.12, y=.35, width=.22, height=.44))
      dev.off()
      
}


### map isotropy vs hindrance

for(drn in c("inbound", "outbound")){
      
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
      
      d <- filter(d, is.finite(windfill),
                  is.finite(isotropy)) %>%
            mutate(isotropy=truncate(isotropy))
      
      d <- d %>% mutate(color = colors2d(cbind(rank(.$isotropy), rank(.$windfill)),
                                         c("cyan", "magenta", "darkred", "darkblue")))
      
      map <- ggplot(d, aes(x, y)) +
            geom_raster(fill=d$color) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
            theme_void() +
            theme(#strip.text=element_text(size=50, angle=-90),
                  text=element_text(color="white"),
                  plot.background = element_rect(fill="black", color="black"))
      legend <- ggplot(d, aes(isotropy, windfill)) +
            geom_point(color=d$color, size=.5) +
            #scale_y_log10(breaks=c(.001, .003, .01, .03, .1, .3, 1)) +
            xlim(.25, NA) +
            theme_minimal() +
            theme(text=element_text(size = 45, color="white"),
                  axis.text=element_text(color="white"),
                  strip.text=element_blank()) +
            labs(x = paste(drn, "windshed isotropy"),
                 y = paste(drn, "proportion windfilling"))
      png(paste0("figures/windsheds/global/isotropy_windfilling_", drn, ".png"), 
          width=3000, height=1500)
      plot(map)
      plot(legend, vp=viewport(x=.12, y=.35, width=.22, height=.44))
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





