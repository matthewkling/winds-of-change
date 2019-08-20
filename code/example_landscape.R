



f <- read_csv("data/windshed/p1_30y_250km.csv") %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_")

ff <- filter(f, abs(y) < 45)

i <- sample(1:nrow(ff), 1)

r <- woc(ff$x[i], ff$y[i], windrose=rose, climate=climate,
         radius = 1000, time_conv=time_conv,
         sigma = 2,
         output = "rasters")

base <- .995


d0 <-  r %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(is.finite(wind_fwd)) %>%
      mutate(windr_fwd = 1/wind_fwd,
             windr_rev = 1/wind_rev,
             wind_fwd = base ^ (1/wind_fwd),
             wind_rev = base ^ (1/wind_rev),
             overlap_fwd = clim_fwd * wind_fwd,
             overlap_rev = clim_rev * wind_rev) 
d <- d0 %>%
      gather(stat, value, -x, -y) %>%
      separate(stat, c("stat", "direction")) %>%
      mutate(stat2 = case_when(stat=="windr" ~ "windr",
                               stat=="clim" ~ "clim",
                               TRUE ~ "wind")) %>%
      group_by(stat2) %>%
      #group_by(stat) %>%
      mutate(#value = ifelse(stat %in% c("clim", "windr"), value, log(value)),
             value = scales::rescale(value)) %>%
      ungroup() %>%
      mutate(stat = factor(stat, 
                           levels=c("windr", "wind", "clim", "overlap"),
                           labels=c("wind time", "wind accessibility", 
                                    "climate similarity", "wind-climate overlap")),
             direction = factor(direction,
                                levels = c("fwd", "rev"),
                                labels = c("outbound direction", "inbound direction")))


maps <- ggplot() +
      geom_raster(data=d, aes(x, y, fill=value)) +
      geom_point(data=ff[i,c("x", "y")], 
                 aes(x, y), color="red", size=1) +
      facet_grid(direction ~ stat#, switch="y"
                 ) + 
      theme_void() +
      scale_fill_viridis_c() +
      theme(legend.position="none",
            strip.text.y=element_text(angle=-90, size=10),
            strip.text.x=element_blank())# +
      #guides(fill = guide_colourbar(barwidth=5, title.position="top", title.hjust = 0.5))



library(lemon)
legend1 <- d0 %>% 
      ggplot(aes(x, y, fill=windr_fwd)) + geom_raster() +
      theme_void() + theme(legend.position="bottom", legend.title = element_text(size=10)) + 
      scale_fill_viridis_c(limits=c(0, max(c(d0$windr_fwd, d0$windr_rev))), breaks=c(0, 500, 1000)) +
      guides(fill = guide_colourbar(barwidth=5, barheight=.5, title.position="top", title.hjust = 0.5)) +
      labs(fill = "wind time (h)")
legend1 <- g_legend(legend1)

legend2 <- d0 %>% 
      ggplot(aes(x, y, fill=wind_fwd)) + geom_raster() +
      theme_void() + theme(legend.position="bottom", legend.title = element_text(size=10)) + 
      scale_fill_viridis_c(limits=0:1, breaks=c(0, .5, 1)) +
      guides(fill = guide_colourbar(barwidth=5, barheight=.5, title.position="top", title.hjust = 0.5)) +
      labs(fill = "wind accessibility")
legend2 <- g_legend(legend2)

legend3 <- d0 %>% 
      ggplot(aes(x, y, fill=clim_fwd)) + geom_raster() +
      theme_void() + theme(legend.position="bottom", legend.title = element_text(size=10)) + 
      scale_fill_viridis_c(limits=0:1, breaks=c(0, .5, 1)) +
      guides(fill = guide_colourbar(barwidth=5, barheight=.5, title.position="top", title.hjust = 0.5)) +
      labs(fill = "climate similarity")
legend3 <- g_legend(legend3)

legend4 <- d0 %>% 
      ggplot(aes(x, y, fill=overlap_fwd)) + geom_raster() +
      theme_void() + theme(legend.position="bottom", legend.title = element_text(size=10)) + 
      scale_fill_viridis_c(limits=0:1, breaks=c(0, .5, 1)) +
      guides(fill = guide_colourbar(barwidth=5, barheight=.5, title.position="top", title.hjust = 0.5)) +
      labs(fill = "wind-climate overlap")
legend4 <- g_legend(legend4)


md <- map_data("world")

origin <- matrix(c(ff$x[i], ff$y[i]), ncol=2)
coords_ll <- SpatialPoints(origin, crs(climate)) %>%
      spTransform(CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
      coordinates()
circle <- geo_circle(coords_ll, width = 1000 * 1000)

key <- ggplot() +
      geom_polygon(data=md, aes(long, lat, group=group),
                   fill="gray30") +
      geom_polygon(data=fortify(circle), aes(long, lat, group=group), fill="dodgerblue") +
      geom_point(data=ff[i,c("x", "y")], 
                 aes(x, y), color="red", size=1) +
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
                 aes(x, y), color="red", size=1) +
      scale_fill_gradientn(colors=c("dodgerblue", "blue", "purple", "red", "yellow")) +
      theme_void() +
      theme_void() + theme(legend.position="top", legend.title = element_text(size=10)) + 
      guides(fill = guide_colourbar(barwidth=5, barheight=.5, title.position="top", title.hjust = 0.5)) +
      labs(fill="temperature (°C)")
temp_legend <- g_legend(temp)
temp <- temp + theme(legend.position="none")

blank <- r <- rectGrob(gp=gpar(fill="white", col="white"))

legends <- arrangeGrob(legend1, legend2, legend3, legend4, blank, nrow=1, widths=c(1, 1, 1, 1, .01))
main <- arrangeGrob(legends, maps, nrow=2, heights=c(1.25, 6))
extras <- arrangeGrob(temp, key, nrow=2, heights=c(1, 1.05))
extras <- arrangeGrob(temp_legend, extras, nrow=2, heights=c(1.25, 6))
p <- arrangeGrob(extras, main, nrow=1, widths=c(1.05, 4.5))


source("E:/edges/range-edges/code/utilities.r")
ggs("figures/manuscript/fig_3.png", p, width=8.5, height=4, units="in", 
    add=grid.text(letters[1:10], 
                  x=rep(seq(.02, .82, length.out=5), each=2), 
                  y=rep(c(.8, .37), 5),
                  gp=gpar(fontsize=12, fontface="bold", col="black")))

###########################

fs <- sample_n(f, 1000)

r <- map2(fs$x, fs$y, woc, windrose=rose, climate=climate,
         radius = 250, time_conv=time_conv,
         sigma = 2,
         output = "rasters")

mwh <- lapply(r, function(x) 1 / x$wind_fwd) %>%
      lapply(values) %>%
      sapply(max,na.rm=T)
max(mwh)
plot(r[[which(!is.finite(mwh))[1]]])

