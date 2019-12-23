

library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)

select <- dplyr::select




# calculate the spatial distribution of climate analogs for a given location
analogs <- function(climate, # raster stack of first, second timesteps
                    coords, # lon-lat vector
                    reverse = FALSE, # should reverse analogs be calculated, instead of forward
                    method = "gaussian", # c("gaussian", "triangular", "threshold)
                    sigma = 1 # width of similarity kernel, in deg C; details vary by method
                    # gaussian = standard deviation, triangular = delta of half suitability,
                    # threshold = cutoff
){
  browser()
  if(reverse) climate <- climate[[2:1]]
  coords <- matrix(coords, ncol=2)
  target <- raster::extract(climate[[1]], coords)
  if(method == "gaussian") kernel <- function(x) exp(-.5*(x/sigma)^2)
  if(method == "triangular") kernel <- function(x) pmax(0, 1 - (abs(x)/sigma/2))
  if(method == "threshold") kernel <- function(x) as.integer(abs(x) <= sigma)
  calc(climate[[2]] - target, kernel)
}


unwrap <- function(x, width=20){
  x <- x %>% crop(extent(180-width, 180, -90, 90)) %>% shift(-360) %>% merge(x)
  x <- x %>% crop(extent(-180, -180+width, -90, 90)) %>% shift(360) %>% merge(x)
  x
}


geo_circle <- function(pts, width, thresh=100) {
  #pts = coords_ll, width = 500000
  # https://stackoverflow.com/questions/25411251/buffer-geospatial-points-in-r-with-gbuffer
  angles <- seq(from = 0, to = 360, by = 5)
  vertices <- geosphere::destPoint(p = pts, b = angles, d = width)
  
  if(pts[1,1] > thresh) vertices[,1] <- ifelse(vertices[,1] > 0, vertices[,1], vertices[,1] + 360)
  if(pts[1,1] < -thresh) vertices[,1] <- ifelse(vertices[,1] < 0, vertices[,1], vertices[,1] - 360)
  
  poly <- Polygon(vertices, hole=F)
  poly <- Polygons(list(poly), ID=NA)
  poly <- SpatialPolygons(Srl = list(poly), 
                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  poly
}

add_coords <- function(windrose){
  rows <- cols <- windrose[[1]]
  rows[] <- rep(1:nrow(rows), each=ncol(rows))
  cols[] <- rep(1:ncol(rows), nrow(rows))
  windrose <- stack(windrose, rows, cols)
  names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
  return(windrose)
}


woc <- function(x, y, windrose, climate,
                radius = 1000, time_conv=identity,
                method = "gaussian", sigma = 2,
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
  s <- list(wind_fwd = accCost(downtrans, origin) %>% "/"(3600) %>% calc(time_conv),
            wind_rev = accCost(uptrans, origin) %>% "/"(3600) %>% calc(time_conv),
            clim_fwd = analogs(clim, coords, method=method, sigma=sigma),
            clim_rev = analogs(clim, coords, method=method, sigma=sigma, reverse = T)) %>%
    stack()
  
  # climate similarity over water is 0
  n_cells <- length(na.omit(values(s$clim_fwd)))
  s[is.na(values(s))] <- 0
  s <- mask(s, circle)
  
  # overlap between wind and climate
  s$overlap_fwd <- s$wind_fwd * s$clim_fwd
  s$overlap_rev <- s$wind_rev * s$clim_rev
  
  if(output == "rasters") return(s)
  
  
  ### summary statistics of wind and climate surfaces
  
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



### test ###
# time_conv <- function(x) 1 / x
# 
# coords <- c(179, 67)
# x <- woc(coords[1], coords[2], rose, climate, time_conv=time_conv)
# 
# x <- land %>% rasterToPoints() %>% as.data.frame() %>%
#   filter(abs(x) <= 180, y > 65) %>%
#   sample_n(1)
# woc(x$x, x$y, rose, climate, time_conv=time_conv)
# x <- woc(x$x, x$y, windrose=rose, climate=climate, 
#          time_conv=function(x) 1/x, radius=250,
#          output="rasters")
# plot(x, colNA="black")




windscapes <- function(data="data/windrose/windrose_p1_wnd10m.tif", 
                       tag="p1_wnd10m",
                       access=list(name = "inv", fx = function(x){1/x}, form = "1/x"),
                       method="gaussian", # c("gaussian", "triangular", "threshold")
                       sigma=2, # climate sigma
                       radius=250, # km
                       water=.1, # accessibility multiplier
                       subsample=NULL, # list(n=1000, seed=12345)
                       output="csv", # c("csv", "df")
                       overwrite=F, 
                       ncores = 7){
  
  # admin
  if(output == "csv"){
    outfile <- paste0("data/windshed/", tag, "_", radius, "km_", access$name, ".csv")
    if(!overwrite & file.exists(outfile)){
      message(paste(access$name, radius, "output already exists -- aborting"))
      return(NULL)
    }
  }
  
  # load current/future climate data
  climate <- stack("data/geographic/processed/temperature.tif") %>% unwrap(180)
  
  # load land data
  land <- raster("f:/cfsr/land.tif") %>% 
    rotate() %>% unwrap(180)
  
  # load windrose data
  rose <- stack(data) %>%
    rotate() %>% unwrap(180)
  
  # downweight conductance over water
  rose <- land %>%
    reclassify(c(NA, NA, water)) %>% # weight=.1 makes it possible to cross narrow waterways
    "*"(rose)
  
  # data frame of locations to model
  pixels <- land %>%
    rasterToPoints() %>%
    as.data.frame() %>%
    filter(abs(x) <= 180)
  if(!is.null(sample)){
    set.seed(subsample$seed)
    pixels <- sample_n(pixels, subsample$n)
  } 
  pixels <- pixels %>% 
    mutate(batch = sample(1:ncores, nrow(.), replace=T)) %>%
    split(.$batch)
  
  # models for each location
  require(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  d <- foreach(pts = pixels,
               .combine="rbind",
               .export=c("woc", "rose", "climate", "geo_circle", "sigma", 
                         "method", "add_coords", "analogs", "access"),
               .packages=(.packages())) %dopar% {
                 map2(pts$x, pts$y, possibly(woc, NULL), 
                      windrose=rose, climate=climate, radius=radius, 
                      method=method, sigma=sigma,
                      time_conv=access$fx) %>%
                   do.call("rbind", .) %>%
                   as.data.frame()
               }
  stopCluster(cl)
  
  if(output == "csv") write_csv(d, outfile)
  if(output == "df") return(mutate(d, tag=tag))
}



##################################################################################



####### main results (defaults for all params) ##########

windscapes()



############ sensitivity analyses ##############

sensitivity_plot <- function(x, maintag, title, outfile, height, 
                             log=F, method="spearman",
                             comptag=NULL){
  
  maintag <- sub(" = ", "_", maintag)
  x$tag <- sub(" = ", "_", x$tag)
  xx <- x
  
  d <- x %>%
    mutate(var = overlap_fwd_windshed_size / clim_fwd_windshed_size) %>%
    select(x, y, tag, var) %>%
    spread(tag, var) %>%
    mutate_(main =  maintag) %>%
    gather(tag, comp, -main, -x, -y) %>%
    mutate(tag=factor(tag, levels=unique(xx$tag))) %>%
    mutate(weight = cos(y/180*pi))
  
  wcor <- function(x, y, weights) cov.wt(cbind(x, y), wt = weights, cor = TRUE)$cor[1,2]
  
  if(log){
    cm <- d %>%
      filter(is.finite(main), is.finite(comp)) %>%
      mutate_at(vars(main, comp), log10) %>%
      group_by(tag) %>%
      summarize(r = wcor(main, comp, weights=weight),
                main = max(10 ^ main, na.rm=T),
                x = 10^((max(comp, na.rm=T)+min(comp, na.rm=T))/2))
  } else{
    cm <- d %>%
      filter(is.finite(main), is.finite(comp)) %>%
      group_by(tag) %>%
      summarize(r = wcor(main, comp, weights=weight),
                main = max(main, na.rm=T),
                x = (max(comp, na.rm=T)+min(comp, na.rm=T))/2)
  }
  
  p <- ggplot(d, aes(comp, main, color = tag==maintag)) +
    facet_grid(. ~ tag, scales="free") +
    geom_point(alpha = 0.2, shape = 16, size = 0.5) +
    geom_text(data=cm, aes(x=x, label=paste("r =", round(r, 2))), 
              color="red", hjust=.5) +
    scale_color_manual(values=c("black", "red")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x="Facilitation ratio",
         y=paste0("Facilitation ratio for\n", maintag))
  
  if(log) p <- p + scale_x_log10() + scale_y_log10()
  ggsave(outfile, p, width=8, height=height, units="in")
  
  
  
  ########
  
  require(patchwork)
  
  if(is.null(comptag)) comptag <- setdiff(unique(d$tag), maintag)[1]
  comptag <- sub(" = ", "_", comptag)
  
  d <- x %>%
    mutate(var = overlap_fwd_windshed_size / clim_fwd_windshed_size) %>%
    select(x, y, tag, var) %>%
    spread(tag, var) %>%
    mutate(diff=.[,maintag] - .[,comptag])
  
  if(log) d <- d %>% mutate(diff=log10(.[,maintag]) - log10(.[,comptag]))
  
  scatter <- ggplot(d, aes_string(maintag, comptag, color="diff")) + 
    geom_abline(slope=1, intercept=0) +
    geom_point()
  if(log) scatter <- scatter + scale_x_log10() + scale_y_log10()
  
  world <- map_data("world")
  map <- ggplot() + 
    geom_polygon(data=world, aes(long, lat, group=group)) +
    geom_point(data=d, aes(x, y, color=diff))
  
  p <- map + scatter + 
    plot_layout(guides="collect") & 
    theme_minimal() &
    scale_color_gradientn(colours=c("blue", "gray", "red"))
  ggsave(outfile %>% str_replace("\\.png", "diffmap.png"), 
         p, width=8, height=4, units="in")
  
}




# p exponents
ws <- list(data=c("data/windrose/windrose_p0_wnd10m.tif",
                  "data/windrose/windrose_p1_wnd10m.tif",
                  "data/windrose/windrose_p2_wnd10m.tif",
                  "data/windrose/windrose_p3_wnd10m.tif"),
           tag=c("p0", "p1", "p2", "p3")) %>%
  pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="p1", comptag="p2", 
                 title="Sensitivity to conductance function",
                 outfile="figures/windsheds/global/sensitivity_scatters_p.png",
                 height=3)
sensitivity_plot(ws, maintag="p1", comptag="p2", 
                 title="Sensitivity to conductance function", log=T,
                 outfile="figures/windsheds/global/sensitivity_scatters_p_log.png",
                 height=3)



# seasonality
ws <- list(data=c("data/windrose/windrose_p1_wnd10m.tif",
                  "data/windrose/windrose_p1_wnd10m_DJF.tif",
                  "data/windrose/windrose_p1_wnd10m_MAM.tif",
                  "data/windrose/windrose_p1_wnd10m_JJA.tif",
                  "data/windrose/windrose_p1_wnd10m_SON.tif"),
           tag=c("annual", "DJF", "MAM", "JJA", "SON")) %>%
  pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="annual", comptag="DJF",
                 title="Sensitivity to season of year",
                 outfile="figures/windsheds/global/sensitivity_scatters_seasons.png",
                 height=2)
sensitivity_plot(ws, maintag="JJA", comptag="DJF",
                 title="Sensitivity to season of year",
                 outfile="figures/windsheds/global/sensitivity_scatters_seasons_v2.png",
                 height=2)
sensitivity_plot(ws, maintag="annual", comptag="DJF",
                 title="Sensitivity to season of year", log=T,
                 outfile="figures/windsheds/global/sensitivity_scatters_seasons_log.png",
                 height=2)


# atmospheric layers
ws <- list(data=c("data/windrose/windrose_p1_wnd10m.tif",
                  "data/windrose/windrose_p1_wnd1000.tif",
                  "data/windrose/windrose_p1_wnd850.tif",
                  "data/windrose/windrose_p1_wnd700.tif"),
           tag=c("z = 10m", "z = 1000hpa", "z = 850hpa", "z = 700hpa")) %>%
  pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="z = 10m", comptag="z = 850hpa", title="Sensitivity to altitude", 
                 outfile="figures/windsheds/global/sensitivity_scatters_altitude.png",
                 height=2.5)
sensitivity_plot(ws, maintag="z = 10m", comptag="z = 850hpa", title="Sensitivity to altitude", log=T, 
                 outfile="figures/windsheds/global/sensitivity_scatters_altitude_log.png",
                 height=2.5)



# accessibility functions
fx <- list(access=list(list(name = "inv", fx = function(x){1/x}, form = "1/x"),
                       list(name = "sqrtinv", fx = function(x) {sqrt(1/x)}, form = "sqrt(1/x)"),
                       list(name = "exp995", fx = function(x) {.995 ^ x}, form = ".995 ^ x"),
                       list(name = "exp99", fx = function(x) {.99 ^ x}, form = ".99 ^ x"),
                       list(name = "exp95", fx = function(x) {.95 ^ x}, form = ".95 ^ x"),
                       list(name = "invlog", fx = function(x) {1/log(x)}, form = "1/log(x)")),
           tag=c("inv", "sqrtinv", "exp995", "exp99", "exp95", "invlog"))
ws <- pmap_df(fx, windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="inv", comptag="exp95",
                 title="Sensitivity to accessibility function",
                 outfile="figures/windsheds/global/sensitivity_scatters_access.png",
                 height=2)
sensitivity_plot(ws %>% na.omit(), maintag="inv", comptag="exp95",
                 title="Sensitivity to accessibility function", log=T,
                 outfile="figures/windsheds/global/sensitivity_scatters_access_log.png",
                 height=2)
curves <- map_df(fx$access, 
                 function(fun) data.frame(name=fun$name, fun=fun$form, 
                                          x = 1:1000, y = fun$fx(1:1000))) %>%
  ggplot(aes(x, y, color=paste0("\"", name, "\": ", fun))) + 
  geom_line() +
  theme_minimal() +
  ylim(0, 1) +
  theme(legend.position=c(.7, .7),
        legend.title = element_blank()) +
  labs(x="wind hours",
       y="wind accessibility") +
  coord_fixed(ratio=1000)
ggsave("figures/windsheds/global/sensitivity_scatters_access_curves.png", curves,
       width=4, height=4, units="in")



# analog function
ws <- expand.grid(method=c("gaussian", "triangular", "threshold"),
                  sigma=c(.5, 1, 2, 3, 4)) %>%
  mutate(tag=paste0(method, "_", sigma)) %>%
  pmap_df(windscapes, data="data/windrose/windrose_p1_wnd10m.tif", 
          output="df", subsample=list(n=1000, seed=12345))
maintag="gaussian_2"
title="Sensitivity to climate analog specificity"
outfile="figures/windsheds/global/sensitivity_scatters_climate.png"
height=4
d <- ws %>%
  mutate(var = overlap_fwd_windshed_size / clim_fwd_windshed_size) %>%
  select(x, y, tag, var) %>%
  spread(tag, var) %>%
  select(-x, -y) %>%
  mutate_(main =  maintag) %>%
  gather(tag, comp, -main) %>%
  separate(tag, c("form", "breadth"), sep="_", remove=F)
cm <- d %>%
  mutate_at(vars(main, comp), log10) %>%
  group_by(breadth) %>%
  mutate(x = 10^((max(comp, na.rm=T)+min(comp, na.rm=T))/2)) %>%
  group_by(breadth, form) %>%
  summarize(r = cor(main, comp, use = "pairwise.complete.obs", method="spearman"),
            main = max(10 ^ main, na.rm=T),
            x=mean(x))
p <- ggplot(d, aes(comp, main, color = tag==maintag)) +
  facet_grid(form ~ breadth, scales="free", labeller = label_both) +
  geom_point(alpha = 0.2, shape = 16, size = 0.5) +
  geom_text(data=cm, aes(x=x, label=round(r, 2)), 
            color="red", hjust=.5, vjust=1) +
  scale_color_manual(values=c("black", "red")) +
  scale_x_log10() + scale_y_log10() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x="Facilitation ratio",
       y=paste0("Facilitation ratio for Gaussian with breadth = 2"))
ggsave(outfile, p, width=8, height=height, units="in")

afun <- function(method = "gaussian", sigma = 1){
  if(method == "gaussian") kernel <- function(x) exp(-.5*(x/sigma)^2)
  if(method == "triangular") kernel <- function(x) pmax(0, 1 - (abs(x)/sigma/2))
  if(method == "threshold") kernel <- function(x) as.integer(abs(x) <= sigma)
  tibble(x=seq(-10, 10, .001),
         y=kernel(x),
         form=method, breadth=sigma)
}
af <- d %>%
  select(form, breadth) %>%
  distinct() %>%
  mutate(sigma = as.numeric(breadth),
         method = form) %>%
  select(method, sigma) %>%
  pmap_df(afun) %>%
  ggplot(aes(x, y, color = form=="gaussian" & breadth==2)) +
  facet_grid(form ~ breadth, labeller = label_both) +
  geom_line() +
  scale_color_manual(values=c("black", "red")) +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x="Climatic difference (degrees C)",
       y="Similarity")
ggsave(outfile %>% sub("\\.png", "_curves.png", .), 
       af, width=8, height=height, units="in")

