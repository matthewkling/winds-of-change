

dtc <- raster("data/geographic/processed/coastal_distance.tif")

f <- read_csv("data/windshed/p1_500km.csv") %>% select(x, y)
coordinates(f) <- c("x", "y")
f$dtc <- raster::extract(dtc, f)
f <- as.data.frame(f)


n <- 100
e <- f %>%
      filter(dtc > 500, abs(y)<75) %>%
      select(x, y) %>%
      sample_n(n) %>%
      mutate(id = 1:nrow(.))

time_conv <- function(x) .995 ^ x
time_conv <- function(x) .99 ^ x
time_conv <- function(x) sqrt(1 / x)
time_conv <- function(x) (1 / x)

r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=time_conv, radius=250,
          output="rasters")

d <- r[1:n] %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)){
      d[[i]]$id <- i
      #d[[i]]$size <- rmd[[i]]["wind_fwd_windshed_size"]
      #d[[i]]$iso <- rmd[[i]]["wind_fwd_windshed_isotropy"]
}
d <- bind_rows(d)
d <- filter(d, is.finite(wind_fwd))


#x <- sample(d$wind_fwd, 10000)
#plot(x, .99 ^ x)

#d <- mutate(d, wind_fwd = 1/(sqrt(wind_fwd)))

#d <- mutate(d, wind_fwd = ifelse(wind_fwd > quantile(wind_fwd, .95),
#                                 quantile(wind_fwd, .95),
#                                 wind_fwd))
#d <- mutate(d, wind_fwd = .995 ^ wind_fwd)


d <- d %>% group_by(id) %>% mutate(size = sum(wind_fwd, na.rm=T)) %>% ungroup() %>%
     mutate(size=as.integer(factor(size)))

# dp <- d %>%
#       group_by(id) %>% sample_n(1)
# library(FNN)
# points <- select(dp, size, iso) %>%
#       ungroup() %>%
#       mutate(px = scales::rescale(rank(size), c(1, 10)),
#              py = scales::rescale(rank(iso), c(1, 10)))
# ranks <- expand.grid(px=1:10, py=1:10)
# 
# for(i in 1:nrow(points)){
#       dists <- get.knnx(select(points, px, py), select(ranks, px, py), k=1)
#       target <- which.max(dists$nn.dist)
#       mover <- dists$nn.index[target]
#       points$px[mover] <- ranks$px[target]
#       points$py[mover] <- ranks$py[target]
# }
# 
# points$loc <- paste(letters[points$px], letters[points$py])
# 
# dd <- left_join(d, points)



p <- ggplot() +
      geom_raster(data=d, aes(x, y, fill=wind_fwd)) +
      facet_wrap(~ size, nrow=10, scales="free") +
      #geom_contour(data=d, aes(x, y, z=wind_fwd), binwidth=max(d$wind_fwd)/50,
      #             color="#290000", alpha=.5, size=1) +
      #geom_vline(data=e[1:n,], aes(xintercept=x), color="white", size=.5) +
      #geom_hline(data=e[1:n,], aes(yintercept=y), color="white", size=.5) +
      # scale_fill_gradientn(colors=c("cyan", "purple", "darkred", "#290000") %>% rev(),
      #                      limits=0:1) +
      scale_fill_viridis_c() +
      guides(fill=guide_colorbar(barwidth=25)) +
      theme_void() + 
      theme(strip.text = element_blank(),
            plot.background=element_rect(fill="black", color="black"),
            text = element_text(size=25, color="white"),
            plot.title=element_text(size=10),
            legend.position="top") +
      labs(fill="windsheds")
ggsave("figures/windsheds/tests/examples_timeconv.png", p,
       width=20, height=20, units="in")


d %>% group_by(id) %>% 
      arrange(desc(wind_fwd)) %>% 
      mutate(wind_cum = cumsum(wind_fwd),
             wind_cum = wind_cum/max(wind_cum),
             lat = abs(mean(y)),
             pix=1:length(wind_cum),
             pix=pix/max(pix))%>%
      ggplot(aes(pix, wind_cum, group=id, color=lat)) +
      geom_line() +
      scale_color_viridis_c()


x=1:1000
xinv=1/x
plot(x, cumsum(xinv))
