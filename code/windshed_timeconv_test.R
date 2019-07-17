


f <- read_csv("data/windshed/p2_500km.csv") %>%
      select(-runtime) %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("property", "direction", "moment", "stat"), sep="_")



n <- 100
e <- f %>%
      select(x, y) %>%
      filter(abs(y)<75) %>%
      sample_n(n) %>%
      mutate(id = 1:nrow(.))


r <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=time_conv, radius=500,
          output="rasters")


d <- r[1:n] %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(d)) d[[i]]$id <- i
d <- bind_rows(d)
d <- filter(d, is.finite(wind_fwd))


#x <- sample(d$wind_fwd, 10000)
#plot(x, .99 ^ x)

#d <- mutate(d, wind_fwd = 1/(sqrt(wind_fwd)))

#d <- mutate(d, wind_fwd = ifelse(wind_fwd > quantile(wind_fwd, .95),
#                                 quantile(wind_fwd, .95),
#                                 wind_fwd))
#d <- mutate(d, wind_fwd = .995 ^ wind_fwd)
hist(d$wind_fwd)

dd <- e
m <- 1

d <- d %>% group_by(id) %>% mutate(size = sum(wind_fwd, na.rm=T)) %>% ungroup() %>%
      mutate(size=as.integer(factor(size)))
      

p <- ggplot() +
      geom_raster(data=d, aes(x, y, fill=wind_fwd)) +
      facet_wrap(~ size, nrow=10, scales="free") +
      #geom_contour(data=d, aes(x, y, z=wind_fwd), binwidth=max(d$wind_fwd)/50,
      #             color="#290000", alpha=.5, size=1) +
      #geom_vline(data=e[1:n,], aes(xintercept=x), color="white", size=.5) +
      #geom_hline(data=e[1:n,], aes(yintercept=y), color="white", size=.5) +
      scale_fill_gradientn(colors=c("cyan", "purple", "darkred", "#290000") %>% rev()) +
      guides(fill=guide_colorbar(barwidth=25)) +
      theme_void() + 
      theme(strip.text = element_blank(),
            plot.background=element_rect(fill="black", color="black"),
            text = element_text(size=25, color="white"),
            plot.title=element_text(size=10),
            legend.position="top") +
      labs(fill="windsheds")
ggsave("figures/windsheds/examples/examples_timeconv.png", p,
       width=20, height=20, units="in")
