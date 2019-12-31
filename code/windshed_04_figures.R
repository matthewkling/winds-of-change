



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
library(rgl)
library(patchwork)

select <- dplyr::select



##infile <- "data/windshed/p1_30y_250km_exp995.csv"
##infile <- "data/windshed/p1_30y_250km_exp99.csv"
infile <- "data/windshed/p1_30y_250km_inv.csv"

source("E:/edges/range-edges/code/utilities.r")

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
geo <- stack(geo[!grepl("temperature|uvma", geo)])
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
                                 labels=c("Climate\navailability", 
                                          "Wind\naccessibility", 
                                          "Wind-climate\noverlap", 
                                          "Wind\nfacilitation")),
            attr = factor(attr, levels=c("latitude", "elevation", "continentality"), 
                          labels=c("Latitude (deg)", "Elevation (m)", "Continentality (km)")))

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

if(F){
  library(mgcv)
  png("figures/windsheds/global/gams.png", width=2000, height=1000)
  par(mfrow=c(2, 4))
  for(drn in c("inbound", "outbound")){
    for(prp in unique(d$property)){
      pd <- d %>%# sample_n(200000) %>% 
        filter(property == prp,
               direction == drn) %>%
        #filter(x < -20, y > 20) %>%
        spread(attr, coord) %>%
        mutate(elevation = as.vector(scale(log(elevation + 1))),
               continentality = as.vector(scale(log(continentality + 1))),
               value=rank(value))
      
      fit <- gam(value ~ s(elevation) + s(continentality) + s(x) + s(y), data=pd)
      #fit <- gam(value ~ elevation + continentality + s(x) + s(y), data=pd)
      vis.gam(fit, theta=220, main=paste(drn, sub("\n", " ", prp)),
              zlim=range(pd$value, na.rm=T))
    }
  }
  dev.off()
  
  png("figures/windsheds/global/gams_lat_elev.png", width=2000, height=1000)
  par(mfrow=c(2, 4))
  for(drn in c("inbound", "outbound")){
    for(prp in unique(d$property)){
      pd <- d %>%# sample_n(200000) %>% 
        filter(property == prp,
               direction == drn) %>%
        #filter(x < -20, y > 20) %>%
        spread(attr, coord) %>%
        mutate(elevation = as.vector(scale(log(elevation + 1))),
               continentality = as.vector(scale(log(continentality + 1))),
               value=rank(value))
      
      fit <- gam(value ~ s(elevation) + s(continentality) + s(x) + s(y), data=pd)
      #fit <- gam(value ~ elevation + continentality + s(x) + s(y), data=pd)
      vis.gam(fit, c("y", "elevation"), theta=130, phi=25, 
              main=paste(drn, sub("\n", " ", prp)),
              zlim=range(pd$value, na.rm=T))
      
      fit <- gam(value ~ te(x, y) + s(elevation, continentality), 
                 data=pd %>% sample_n(5000))
      vis.gam(fit, c("continentality", "elevation"), theta=130)
    }
  }
  dev.off()
}

for(prp in unique(d$property)){
  pd <- d %>% sample_n(200000) %>% filter(property == prp) %>% 
    mutate(coord = case_when(attr == "Elevation (m)" ~ log10(coord),
                             attr == "Continentality (km)" ~ log10(coord),
                             TRUE ~ coord),
           coord = truncate(coord, .001, "low"),
           value = truncate(value, .005, "high"))
  
  title <- case_when(prp == "Climate\navailability" ~ "Climate analog availability",
                     prp == "Wind\naccessibility" ~ "Wind accessibility (1/h)",   
                     prp == "Wind-climate\noverlap" ~ "Wind-climate overlap (1/h)",
                     prp == "Wind\nfacilitation" ~ "Wind facilitation (1/h)")
  
  style <- theme(axis.title=element_text(color="black", size=25),
                 legend.position="bottom",
                 strip.text=element_text(size=25, color="black"),
                 legend.text=element_text(size=25, color="black"),
                 axis.text=element_text(size=15, color="black"),
                 legend.title=element_blank())
  
  plat <- pd %>% filter(attr == "Latitude (deg)") %>%
    ggplot(aes(coord, value)) +
    geom_point(size=.2, color="black") +
    geom_smooth(se=F, aes(color=paste(direction, "     ")), size=2) +
    scale_color_manual(values=c("red", "dodgerblue")) +
    theme_minimal() + style +
    labs(y = title,
         x = "Latitude (deg)") +
    ylim(0, NA)
    
  
  pelev <- pd %>% filter(attr == "Elevation (m)") %>%
    ggplot(aes(coord, value)) +
    geom_point(size=.2, color="black") +
    geom_smooth(se=F, aes(color=paste(direction, "     ")), size=2) +
    scale_color_manual(values=c("red", "dodgerblue")) +
    theme_minimal() + style +
    labs(y = gsub("\\\n", " ", prp),
         x = "Elevation (m)") +
    scale_x_continuous(labels = function(x) ifelse(x < 4, 10^x, x)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ylim(0, NA)
  
  pcont <- pd %>% filter(attr == "Continentality (km)") %>%
    ggplot(aes(coord, value)) +
    geom_point(size=.2, color="black") +
    geom_smooth(se=F, aes(color=paste(direction, "     ")), size=2) +
    scale_color_manual(values=c("red", "dodgerblue")) +
    theme_minimal() + style +
    labs(y = gsub("\\\n", " ", prp),
         x = "Continentality (km)") +
    scale_x_continuous(labels = function(x) ifelse(x < 4, 10^x, x)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ylim(0, NA)
  
  
  p <- (plat + pelev + pcont) / guide_area() + 
    plot_layout(guides = "collect", heights=c(1, .1))
  
  ggsave(paste0("figures/windsheds/global/scatter_windshed_geography_bw_",
                gsub("\\\n", "", prp), ".png"), p,
         bg="black", width=15, height=6, units="in")
}



######################

## correlations between wind regimes, local alignment, and windshed stats

# wind regime data
land <- raster("f:/cfsr/land.tif")
w <- stack("data/wind_regime/regimes.tif") %>%
  mask(land) %>%
  rotate()
names(w) <- c("speed", "direction", "anisotropy", "u", "v")
w <- as.data.frame(rasterToPoints(w))

# local alignment data
a <- stack("data/alignment/alignment_s3.tif")
names(a) <- c("wx", "wy", "wa", "wm", "ta", "land")
a <- a %>% rasterToPoints %>% 
  as.data.frame() %>%
  mutate(theta = abs(wa - ta),
         theta = ifelse(theta > 180, 360-theta, theta)) %>%
  filter(is.finite(theta)) %>%
  select(x, y, theta)

dw <- w %>% 
  select(x, y, speed, anisotropy) %>% mutate(x=round(x, 2), y=round(y, 2)) %>%
  left_join(f %>% mutate(x=round(x, 2), y=round(y, 2))) %>%
  left_join(a %>% mutate(x=round(x, 2), y=round(y, 2))) %>%
  mutate(asr = anisotropy / speed) %>%
  gather(var, value, -x, -y, -speed, -anisotropy, -asr, -theta) %>%
  separate(var, c("property", "direction", "moment", "stat"), sep="_") %>%
  mutate(direction = ifelse(direction=="fwd", "outbound", "inbound")) %>%
  filter(moment == "windshed", stat == "size") %>%
  gather(attr, coord, speed, anisotropy, asr, theta)

d <- spread(dw, property, value) %>%
  mutate(windfill = overlap / clim) %>%
  filter(wind > 0) %>%
  gather(property, value, clim, wind, overlap, windfill)

# d <- mutate(d, property = factor(property, levels=c("clim", "wind", "overlap", "windfill"), 
#                                  labels=c("climate\nsimilarity", "wind\naccessibility", 
#                                           "climate-wind\noverlap", "wind\nfacilitation")))
#             #attr = factor(attr, levels=c("latitude", "elevation", "continentality")))

p <- ggplot(d %>% sample_n(200000) %>% #sample_n(200000) %>% 
              filter(property %in% c("wind", "windfill")), 
            aes(coord, value)) +
  facet_grid(property ~ attr, scales="free") +
  geom_point(size=.2) +
  geom_smooth(se=F, aes(color=direction)) +
  scale_color_manual(values=c("dodgerblue", "red")) +
  scale_x_log10() + scale_y_log10() +
  theme_minimal() +
  labs(title = "Wind regime relationships with forward windshed summary statistics") +
  theme(axis.title=element_blank(),
        text=element_text(size=25))
ggsave("figures/windsheds/global/scatter_windshed_regime.png", 
       width=12, height=8, units="in")



dw <- w %>% 
  select(x, y, speed, anisotropy, v) %>% 
  mutate(x=round(x, 2), y=round(y, 2)) %>%
  left_join(f %>% 
              mutate(x=round(x, 2), y=round(y, 2),
                     windfill = overlap_fwd_windshed_size/clim_fwd_windshed_size) %>%
              select(x, y, windfill, wind_fwd_windshed_size)) %>%
  left_join(a %>% mutate(x=round(x, 2), y=round(y, 2))) %>%
  rename(wind = wind_fwd_windshed_size) %>%
  mutate(v = v * sign(y)) %>%
  mutate(y = abs(y)) %>%
  mutate(asr = anisotropy/speed) %>%
  mutate(speed=log(speed), anisotropy=-log(anisotropy), asr=log(asr),
         windfill=log(windfill)) %>%
  mutate_all(scale) %>%
  na.omit()

lm(windfill ~ speed + speed^2, data=dw) %>% summary()
lm(windfill ~ speed + speed^2 + asr + asr ^ 2 + speed*asr, data=dw) %>% summary()
lm(windfill ~ speed + speed^2 + asr + asr ^ 2 + speed*asr + v + v^2, data=dw) %>% summary()

lm(windfill ~ speed*anisotropy, data=dw) %>% summary()
lm(windfill ~ speed*anisotropy*theta, data=dw) %>% summary()



ggplot(dw %>% sample_n(2000), aes(speed, asr)) + geom_point()
ggplot(dw %>% sample_n(2000), aes(speed, resid)) + geom_point()
ggplot(dw %>% sample_n(2000), aes(asr, resid)) + geom_point()

library(mgcv)
fit <- gam(windfill ~ s(speed) + s(asr) + s(v), data=dw %>% sample_n(10000))
vis.gam(fit, theta=-45, view=c("asr", "v"))


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
         direction == "outbound") %>%
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

# d <- f %>%
#       filter(property == "wind",
#              direction == "outbound",
#              moment == "windshed") %>%
#       spread(stat, value)







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
                  clim="Climate\nanalog\navailability", 
                  wind="Wind\naccessibility\n(1/h)", 
                  overlap="Climate-wind\noverlap\n(1/h)", 
                  windfill="Wind\nfacilitation\n(1/h)")
  
  map0 <- ggplot(v, aes(x, y, fill=value)) +
    facet_grid(direction ~ .) +
    geom_raster() +
    scale_fill_gradientn(colors=c("darkblue", "red", "yellow"),
                         trans=trans) +
    theme_minimal() +
    theme(text=element_text(size=20, color="black"),
          strip.text=element_text(size=20, angle=-90),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          #plot.background = element_rect(fill="black", color="black"),
          legend.position=c(.07, .73)) +
    guides(fill=guide_colorbar(barheight=12)) +
    scale_x_continuous(expan=c(0,0)) +
    scale_y_continuous(expan=c(0,0), limits=c(-90, 90)) +
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
  
  # assemble multi-panel figure for SI
  lbl <- switch(var, 
                clim="Climateavailability", 
                wind="Windaccessibility", 
                overlap="Wind-climateoverlap", 
                windfill="Windfacilitation")
  library(png)
  img <- paste0("figures/windsheds/global/scatter_windshed_geography_bw_",
                lbl, ".png") %>%
    readPNG() %>%
    rasterGrob(interpolate=TRUE)
  p <- arrangeGrob(map0, img, nrow=2, heights=c(12, 5))
  
  ggs(paste0("figures/manuscript/SI_fig_", var, ".png"),
      p, width=12, height=17, units="in",
      add = grid.text(letters[1:5], 
                      x=c(.03, .03, .09, .41, .73), 
                      y=c(.97, .63, .26, .26, .26),
                      gp=gpar(fontsize=30, fontface="bold", col="black")))
}



## 2d map clim vs windfilling

dd <- d %>%
  select(-moment, -stat) %>%
  spread(property, value) %>%
  mutate(windfill = truncate(windfill)) %>%
  filter(is.finite(windfill), is.finite(clim)) %>%
  mutate(color = colors2d(cbind(rank(.$clim), rank(.$windfill)),
                          c("forestgreen", "yellow", "red", "darkblue"))) %>%
  mutate(direction = factor(direction, levels=c("outbound", "inbound")))

txt <- data.frame(direction = c("outbound", "inbound"),
                  text = c("outbound", "inbound"),
                  x = min(dd$x) + 25, y=min(dd$y) + .2 *(diff(range(dd$y))))

map <- ggplot(dd, aes(x, y)) +
  facet_grid(direction ~ .) +
  geom_raster(fill=dd$color) +
  geom_text(data=txt, aes(x, y, label=text), 
            hjust=0, size=6, lineheight=.75) +
  ylim(-90, 90) +
  theme_void() +
  theme(strip.text=element_blank())
lgnd <- ggplot(dd, aes(clim, windfill)) +
  geom_point(color=dd$color, size=.2) +
  xlim(0, .75) +
  theme_minimal() +
  theme(strip.text=element_blank()) +
  labs(x = "Climate availability",
       y = "Wind facilitation (1/h)")

sd <- d %>%
  select(-moment, -stat) %>%
  filter(property %in% c("clim", "wind", "overlap")) %>%
  mutate(property = factor(property, levels=c("clim", "wind", "overlap"),
                           labels=c("Climate  \navailability  ",
                                    "Wind  \naccessibility  \n(1/h)  ",
                                    "Wind-climate  \noverlap  \n(1/h)  "))) %>%
  group_by(property) %>%
  mutate(value = ifelse(grepl("overlap", property), truncate(value), value)) %>%
  ungroup() %>%
  spread(direction, value)


tsd <- sd %>%
  group_by(property) %>%
  summarize(r2pears = cor(inbound, outbound, use="pairwise.complete.obs") ^2,
            r2spear = cor(inbound, outbound, use="pairwise.complete.obs", method="spearman") ^2,
            inbound = max(inbound),
            outbound = max(range(outbound)))

scat <- sd %>% #sample_n(500) %>%
  ggplot(aes(outbound, inbound)) +
  facet_wrap(~property, ncol=1, scales="free") +
  geom_point(size=.25, alpha=.05) +
  geom_text(data=tsd, aes(label=property), 
            hjust=1, vjust=1, lineheight=.7, fontface="bold", size=4) +
  geom_text(data=tsd, aes(label=paste("\n\n\n\n\n\n\n\n\n\nr2 =", round(r2pears, 2))), color="red",
            hjust=1, vjust=1, lineheight=.7, fontface="bold", size=4) +
  theme_minimal() +
  theme(strip.text=element_blank()) +
  labs(x="outbound",
       y="inbound") +
  scale_x_continuous(labels=function(x) ifelse(x %in% c(0.0025, 0.0075, .0125), "", x),
                     limits = c(0, NA)) +
  scale_y_continuous(labels=function(x) ifelse(x %in% c(0.0025, 0.0075, .0125), "", x),
                     limits = c(0, NA))

p <- arrangeGrob(lgnd, scat, ncol=1, heights=c(1, 3))
p <- arrangeGrob(map, p, ncol=2, widths=c(4, 1))

source("E:/edges/range-edges/code/utilities.r")
ggs("figures/windsheds/global/windfill_clim.png", 
    p, width=10, height=8, units="in",
    add = grid.text(letters[1:6], 
                    x=c(.05, .05, .82, .82, .82, .82), 
                    y=c(.62, .115, .78, .55, .30, .05),
                    gp=gpar(fontsize=20, fontface="bold", col="black")))
file.copy("figures/windsheds/global/windfill_clim.png",
          "figures/manuscript/SI_fig_windfillclim.png", overwrite = T)



dg <- function(lat) distGeo(c(0, lat), c(1, lat))/1000
areas <- d %>% select(y) %>% distinct()
areas$area <- sapply(areas$y, dg)





### climate vs windfilling

for(drn in c("inbound", "outbound")){
  
  dd <- d %>%
    filter(direction==drn) %>%
    select(-moment, -stat) %>%
    spread(property, value) %>%
    mutate(windfill = overlap / clim) %>%
    filter(is.finite(clim), is.finite(windfill)) %>%
    mutate(clim = truncate(clim, .001),
           windfill = truncate(windfill, .001, c("high", "low"))) %>%
    
    left_join(areas) %>%
    arrange(clim) %>%
    mutate(climrnk = cumsum(area)/ sum(area)) %>%
    arrange(windfill) %>%
    mutate(windfillrnk = cumsum(area) / sum(area)) %>%
    
    mutate(color2 = colors2d(cbind(.$climrnk, .$windfillrnk),
                             c("forestgreen", "yellow", "red", "darkblue"))) %>%
    
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
  
  
  map <- ggplot(dd, aes(x, y)) +
    geom_raster(fill=dd$color2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
    theme_void()
  legend <- ggplot(dd, aes(clim, windfill)) +
    geom_point(color=dd$color2, size=.3) +
    theme_minimal() +
    theme(text=element_text(size = 45)) +
    labs(x = "area of analog climate",
         y = "proportion wind-accessible")
  png(paste0("figures/windsheds/global/windfill_clim_", drn, "_white.png"), 
      width=3000, height=1500)
  plot(map)
  plot(legend, vp=viewport(x=.12, y=.35, width=.22, height=.44))
  dev.off()
  
  updown <- ifelse(drn=="outbound", "down", "up")
  
  legend <- ggplot(dd, aes(clim, windfill)) +
    geom_point(color=dd$color2, size=.3) +
    geom_vline(xintercept=max(dd$clim[dd$syndrome=="FALSE FALSE"]), linetype=2, color="white") +
    geom_hline(yintercept=max(dd$windfill[dd$syndrome=="FALSE FALSE"]), linetype=2, color="white") +
    scale_x_continuous(expand=c(0,0), limits=c(0, NA)) +
    scale_y_continuous(expand=c(0,0), trans="log10") +
    theme_minimal() +
    theme(text=element_text(size = 17),
          axis.text.y=element_text(angle=90),
          legend.position="none") +
    labs(x = paste0("Analog climate availability\n"),
         y = paste0("\nWind facilitation (1/h)"))
  
  
  
  dd <- dd %>%
    mutate(syndrome = paste(climrnk > .5, windfillrnk > .5),
           syndrome = factor(syndrome, 
                             levels = c("TRUE TRUE", "FALSE TRUE", 
                                        "FALSE FALSE", "TRUE FALSE")))
  
  palette <- c("forestgreen", "darkblue", "red", "yellow")
  hist <- dd %>%
    mutate(y = round(y)) %>%
    group_by(y, syndrome) %>%
    summarize(n = sum(area)) %>%
    ggplot(aes(x=y, y=n, fill=syndrome)) +
    geom_histogram(stat="identity", position="fill", width=1) +
    coord_flip() +
    scale_fill_manual(values=palette) +
    scale_y_continuous(expand=c(0,0), 
                       breaks=c(0, .5, 1)) +
    scale_x_continuous(expand=c(0,0), limits=c(-90, 90), 
                       breaks=seq(-90, 90, 30)) +
    theme(legend.position="none",
          text=element_text(size = 17),
          panel.grid=element_blank(),
          panel.background = element_rect(fill="white")) +
    labs(y = "Proportion of area", x = "°N")
  
  library(patchwork)
  
  
  library(png)
  img <- readPNG(paste0("figures/windsheds/rgl/climwind_", drn, ".png"))
  img <- apply(img, c(1, 2), function(x){if(mean(x)>.75){return(c(1,1,1))}; return(x)})
  img <- aperm(img, c(2, 3, 1))
  img <- rasterGrob(img, interpolate=TRUE)
  
  p <- (legend + img) / 
    (hist + map +  plot_layout(widths=c(1, 7)))
    
  
  source("E:/edges/range-edges/code/utilities.r")
  ggs(paste0("figures/windsheds/global/windfill_clim_", drn, "_array.png"),
      p, width=12, height=12, units="in", dpi=1000,
      add = list(grid.text(letters[1:3], 
                           x=c(.03, .03, .23), 
                           y=c(.95, .47, .47),
                           gp=gpar(fontsize=30, fontface="bold", col="black")),
                 grid.text(letters[4], 
                           x=c(.56), 
                           y=c(.95),
                           gp=gpar(fontsize=30, fontface="bold", col="white"))))
  
  if(drn=="outbound") file.copy(paste0("figures/windsheds/global/windfill_clim_", drn, "_array.png"),
                                "figures/manuscript/fig_4.png",
                                overwrite=T)
}


### map anisotropy vs hindrance

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
    mutate(isotropy=truncate(isotropy),
           windfill = truncate(windfill, .001, c("high", "low")))
  
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
  legend <- ggplot(d, aes(1 - isotropy, windfill)) +
    geom_point(color=d$color, size=.5) +
    scale_y_sqrt() +
    xlim(.25, NA) +
    theme_minimal() +
    theme(text=element_text(size = 45, color="white"),
          axis.text=element_text(color="white"),
          strip.text=element_blank()) +
    labs(x = paste(drn, "windshed anisotropy"),
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
ggsave("figures/manuscript/SI_fig_tempkernel.png", 
       width=6, height=4, units="in")







