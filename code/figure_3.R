
# grab functions from windshed_02_analysis.r before running


f <- read_csv("data/windshed/p1_30y_250km_inv.csv") %>%
  select(-runtime) %>%
  gather(var, value, -x, -y) %>%
  separate(var, c("property", "direction", "moment", "stat"), sep="_")

ff <- filter(f, abs(y) < 45)






i <- sample(1:nrow(ff), 1)
#saveRDS(i, "data/windshed/example_landscape_index.rds")
#i <- readRDS("data/windshed/example_landscape_index.rds")

# 183380, 2083858
i <- 2083858




# load current/future climate data
climate <- stack("data/geographic/processed/temperature.tif") %>% unwrap(180)

# load land data
land <- raster("f:/cfsr/land.tif") %>% 
  rotate() %>% unwrap(180)

# load windrose data
rose <- stack("data/windrose/windrose_p1_wnd10m.tif") %>%
  rotate() %>% unwrap(180)

# downweight conductance over water
water=.1
rose <- land %>%
  reclassify(c(NA, NA, water)) %>% # weight=.1 makes it possible to cross narrow waterways
  "*"(rose)

access=list(name = "inv", fx = function(x){1/x}, form = "1/x")
time_conv=identity

radius = 250
r <- woc(ff$x[i], ff$y[i], windrose=rose, climate=climate,
         radius = radius, time_conv=time_conv,
         sigma = 2,
         output = "rasters")


d <-  r %>%
  rasterToPoints() %>%
  as.data.frame() %>%
  filter(is.finite(wind_fwd)) %>%
  mutate(windr_fwd = 1/wind_fwd,
         windr_rev = 1/wind_rev,
         overlap_fwd = clim_fwd * windr_fwd,
         overlap_rev = clim_rev * windr_rev) %>%
  gather(stat, value, -x, -y) %>%
  separate(stat, c("stat", "direction")) %>%
  spread(stat, value)

pointsize <- 2

keys <- guides(fill = guide_colourbar(barwidth=5, barheight=.5, 
                                      title.position="top", title.hjust = 0.5))
style <- theme(legend.position="top",
               strip.text.y=element_text(angle=-90, size=10),
               strip.text.x=element_blank())
style <- theme(legend.position="top",
               strip.text=element_blank())
origins <- geom_point(data=ff[i,c("x", "y")], aes(x, y), color="red", size=pointsize)

clim <- d %>%
  ggplot() +
  geom_raster(aes(x, y, fill=clim)) +
  facet_grid(direction ~ .) + 
  theme_void() + style + keys + origins +
  scale_fill_gradientn(colours=c("black", "yellow"),
                       limits=0:1, breaks=c(0,.5,1)) +
  labs(fill = "Climate similarity")

time <- d %>% 
  ggplot() +
  geom_raster(aes(x, y, fill=wind)) +
  facet_grid(direction ~ .) + 
  theme_void() + style + keys + origins +
  scale_fill_gradientn(colours=c("cyan", "black"),
                       breaks=c(0, 400, 800)) +
  labs(fill = "Wind time (h)")

access <- d %>% 
  ggplot() +
  geom_raster(aes(x, y, fill=windr)) +
  facet_grid(direction ~ .) + 
  theme_void() + style + keys + origins +
  scale_fill_gradientn(colours=c("black", "cyan"), trans="sqrt", 
                       limits=c(0, NA), breaks=c(0, .01, .04)) +
  labs(fill = "Accessibility (1/h)")

overlap <- d %>% mutate(direction = factor(direction, levels=c("fwd", "rev"),
                                           labels=c("Outbound", "Inbound"))) %>%
  ggplot() +
  geom_raster(aes(x, y, fill=overlap)) +
  facet_grid(direction ~ .) + 
  theme_void() + style + keys + origins +
  scale_fill_gradientn(colours=c("black", "green"), trans="sqrt",
                       limits=c(0, NA), breaks=c(0, .01, .04)) +
  labs(fill = "Overlap (1/h)") +
  theme(strip.text = element_text(size=14, angle=-90))




md <- map_data("world")

origin <- matrix(c(ff$x[i], ff$y[i]), ncol=2)
coords_ll <- SpatialPoints(origin, crs(climate)) %>%
  spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  coordinates()
circle <- geo_circle(coords_ll, width = radius * 1000)

key <- ggplot() +
  geom_polygon(data=md, aes(long, lat, group=group),
               fill="gray30") +
  geom_polygon(data=fortify(circle), aes(long, lat, group=group), 
               fill="white") +
  geom_point(data=ff[i,c("x", "y")], 
             aes(x, y), color="red", size=pointsize/2) +
  coord_map(projection="ortho", 
            orientation=c(ff[i, "y"], ff[i, "x"], 0)) +
  scale_x_continuous(breaks=seq(-180, 180, 20)) +
  scale_y_continuous(breaks=seq(-90, 90, 15)) +
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_line(color="gray40"))

td <- climate[[1]] %>% 
  crop(circle) %>% mask(circle) %>%
  rasterToPoints() %>% as.data.frame()
temp <- ggplot() +
  geom_raster(data=td, aes(x, y, fill=layer.1)) +
  geom_point(data=ff[i,c("x", "y")], 
             aes(x, y), color="red", size=pointsize) +
  scale_fill_gradientn(colors=c("cyan", "dodgerblue", "blue", "purple", "red", "orange", "gold", "yellow"),
                       breaks=c(10, 14, 18)) +
  theme_void() +
  theme_void() + theme(legend.position="top", legend.title = element_text(size=10)) + 
  guides(fill = guide_colourbar(barwidth=5, barheight=.5, title.position="top", title.hjust = 0.5)) +
  labs(fill="Temperature (°C)")


library(patchwork)

p <- (temp / key) | clim | time | access | overlap 

source("E:/edges/range-edges/code/utilities.r")

library(grid)

ggs("figures/manuscript/fig_3.pdf", 
    p, width=8.5, height=4.3, units="in", 
    add=list(grid.text(letters[1:10], 
                  x=rep(seq(.04, .8, length.out=5), each=2), 
                  y=rep(c(.8, .39), 5),
                  gp=gpar(fontsize=12, fontface="bold", col="black")),
             grid.rect(x=.608, y=c(.24, .65), width=.78, height=.39,
                      gp=gpar(col="black", fill=rgb(1,1,1,0), lwd=.25))
             ))

