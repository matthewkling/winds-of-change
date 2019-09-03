

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
                    sigma = 1 # standard deviation of similarity kernel, in degrees Celsius
){
      if(reverse) climate <- climate[[2:1]]
      coords <- matrix(coords, ncol=2)
      target <- raster::extract(climate[[1]], coords)
      kernel <- function(x) exp(-.5*(x/sigma)^2)
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
                sigma = 2,
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
                clim_fwd = analogs(clim, coords, sigma=sigma),
                clim_rev = analogs(clim, coords, sigma=sigma, reverse = T)) %>%
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


############


land <- raster("f:/cfsr/land.tif") %>% 
      rotate() %>% unwrap(180)

# load windrose data
rose <- stack("data/windrose/windrose_p1_1980_2009.tif") %>%
      rotate() %>% unwrap(180)

# downweight conductance over water
rose <- land %>%
      reclassify(c(NA, NA, .1)) %>% # weight=.1 makes it possible to cross narrow waterways
      "*"(rose)

# load current/future climate data
climate <- stack("data/geographic/processed/temperature.tif") %>% unwrap(180)


#time_conv <- function(x) .995 ^ x
time_conv <- function(x) 1 / x

# test
coords <- c(179, 67)
x <- woc(coords[1], coords[2], rose, climate, time_conv=time_conv)
#profvis({ x <- woc(coords[1], coords[2], rose, climate, time_conv=identity) })

x <- land %>% rasterToPoints() %>% as.data.frame() %>%
      filter(abs(x) <= 180, y > 65) %>%
      sample_n(1)
woc(x$x, x$y, rose, climate, time_conv=time_conv)
x <- woc(x$x, x$y, windrose=rose, climate=climate, 
     time_conv=function(x) 1/x, radius=250,
     output="rasters")
plot(x, colNA="black")


windscapes <- function(fun, ncores = 7, radius=250, overwrite=F){
      
      outfile <- paste0("data/windshed/p1_30y_", radius, "km_", fun$name, ".csv")
      if(!overwrite & file.exists(outfile)){
            message(paste(fun$name, radius, "output already exists -- aborting"))
            return(NULL)
      }
      
      pixels <- land %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            filter(abs(x) <= 180) %>%
            #sample_n(10) %>%
            mutate(batch = sample(1:ncores, nrow(.), replace=T)) %>%
            split(.$batch)
      
      require(doParallel)
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      d <- foreach(pts = pixels,
                   .combine="rbind",
                   .export=c("woc", "rose", "climate", "geo_circle", 
                             "add_coords", "analogs", "fun"),
                   .packages=(.packages())) %dopar% {
                         map2(pts$x, pts$y, possibly(woc, NULL), 
                              windrose=rose, climate=climate, radius=radius,
                              time_conv=fun$fx) %>%
                               do.call("rbind", .) %>%
                               as.data.frame()
                   }
      stopCluster(cl)
      write_csv(d, outfile)
}

functions <- list(list(name = "inv", fx = function(x){1/x}, form = "1/x"),
                  list(name = "sqrtinv", fx = function(x) {sqrt(1/x)}, form = "sqrt(1/x)"),
                  list(name = "exp995", fx = function(x) {.995 ^ x}, form = ".995 ^ x"),
                  list(name = "exp99", fx = function(x) {.99 ^ x}, form = ".99 ^ x"),
                  list(name = "exp98", fx = function(x) {.95 ^ x}, form = ".95 ^ x"),
                  list(name = "invlog", fx = function(x) {1/log(x)}, form = "1/log(x)"))

if(F) lapply(functions[1], windscapes, overwrite=T)



####### sensitivity to accessibility function #########

library(ggforce)
library(gridExtra)
library(grid)
source("E:/edges/range-edges/code/utilities.r")

ncores <- 7
pixels <- land %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(abs(x) <= 180) %>%
      sample_n(1000) %>%
      mutate(batch = sample(1:ncores, nrow(.), replace=T)) %>%
      split(.$batch)


windscape_sample <- function(fun, pixels, ncores = 7, radius=250){
      
      require(doParallel)
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      d <- foreach(pts = pixels,
                   .combine="rbind",
                   .export=c("woc", "rose", "climate", "geo_circle", 
                             "add_coords", "analogs"),
                   .packages=(.packages())) %dopar% {
                         map2(pts$x, pts$y, woc, # possibly(woc, NULL), 
                              windrose=rose, climate=climate, radius=radius,
                              time_conv=fun$fx) %>%
                               do.call("rbind", .) %>%
                               as.data.frame()
                   }
      stopCluster(cl)
      d$fun <- fun$name
      return(d)
}

sd <- lapply(functions, windscape_sample, pixels=pixels)

d <- do.call("rbind", sd) %>%
      select(overlap_fwd_windshed_size, fun, x, y)

fd <- data.frame()
for(i in 1:length(functions)){
      fun <- functions[[i]]
      d$form[d$fun==fun$name] <- fun$form
      fdi <- data.frame(name=fun$name, fun=fun$form,
                        x = 1:1000, y = fun$fx(1:1000))
      fd <- rbind(fd, fdi)
}

d <- d %>%
      arrange(form) %>%
      select(-form) %>%
      spread(fun, overlap_fwd_windshed_size) %>%
      select(-x, -y) %>%
      mutate_all(rank)

colnames(d) <- sapply(functions, function(x)x$form)
scatters <- ggplot(d, aes(x = .panel_x, y = .panel_y)) + 
      geom_point(alpha = 0.2, shape = 16, size = 0.5) + 
      facet_matrix(vars(everything())) +
      theme_minimal() +
      theme(axis.text=element_blank())

curves <- ggplot(fd, aes(x, y, color=fun)) + geom_line() +
      theme_minimal() +
      ylim(0, 1) +
      theme(legend.position=c(.7, .7),
            legend.title = element_blank()) +
      labs(x="wind hours",
           y="wind accessibility") +
      coord_fixed(ratio=1000)

p <- arrangeGrob(curves, scatters, ncol=2, widths=c(1, 2))

corr <- as.vector(cor(d, method="spearman")) %>% round(2)

ggs("figures/windsheds/global/sensitivity.png", p, width=9, height=6, units="in",
    add = list(grid.text(letters[1:2], 
                         x=c(.03, .34), 
                         y=c(.75, .95),
                         gp=gpar(fontsize=20, fontface="bold", col="black")),
               grid.text(corr, 
                         x=rep(seq(.37, .89, length.out=6), 6), 
                         y=rep(seq(.92, .14, length.out=6), each=6),
                         gp=gpar(fontsize=8, col="red"))
    )
    )
file.copy("figures/windsheds/global/sensitivity.png",
          "figures/manuscript/SI_fig_sensitivity.png")







### define a quantile function based on the global frequency of wind hours

n <- 100
e <- pixels %>%
      bind_rows() %>%
      select(x, y) %>%
      sample_n(n) %>%
      mutate(id = 1:nrow(.))
e <- map2(e$x, e$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, 
          time_conv=identity, radius=250,
          output="rasters")
quants <- map(e, function(x) x$wind_fwd) %>%
      map(values) %>%
      do.call("c", .) %>%
      na.omit() %>%
      quantile(probs=seq(0, 1, .001))
# functions[[length(functions) + 1]] <- list(name = "quant", 
#                                            fx = function(x) {ecdf(quants)(x)}, 
#                                            form = "1-ecdf(x)")


functions <- list(list(name = "inv", fx = function(x){1/x}, form = "1/x"),
                  list(name = "sqrtinv", fx = function(x) {sqrt(1/x)}, form = "sqrt(1/x)"),
                  list(name = "exp995", fx = function(x) {.995 ^ x}, form = ".995 ^ x"),
                  list(name = "exp990", fx = function(x) {.990 ^ x}, form = ".99 ^ x"),
                  list(name = "exp975", fx = function(x) {.975 ^ x}, form = ".975 ^ x"),
                  list(name = "invlog", fx = function(x) {1/log(x)}, form = "1/log(x)"))


fd <- data.frame()
for(i in 1:length(functions)){
      fun <- functions[[i]]
      d$form[d$fun==fun$name] <- fun$form
      fdi <- data.frame(name=fun$name, fun=fun$form,
                        x = 1:1000, y = fun$fx(1:1000))
      fd <- rbind(fd, fdi)
}

qd <- data.frame(y=seq(0, 1, length.out=length(quants)), x=quants) %>%
      mutate(#y = cumsum(y),
             y = y/max(y)) %>%
      filter(x <= 1000)
# qd <- data.frame(q=round(quants, -1)) %>%
#       count(q) %>%
#       mutate(n = n / max(n)) %>%
#       filter(q <= 1000)

ggplot() + 
      geom_line(data=qd, aes(x, y)) +
      geom_line(data=fd, aes(x, y, color=fun)) +
      theme_minimal() +
      ylim(0, 1) +
      xlim(0, max(fd$x)) +
      theme(legend.position=c(.7, .7),
            legend.title = element_blank()) +
      labs(x="wind hours",
           y="wind accessibility") +
      coord_fixed(ratio=1000)




stop("wootwoot")


### some exploratory plots ###

pts <- pixels[[1]][1:49,]
r <- map2(pts$x, pts$y, possibly(woc, NULL), 
          windrose=rose, climate=climate, cost_to_flow=cost_to_flow, radius=500,
          output="rasters") %>%
      lapply(function(x) x$wind_fwd)
po <- map2(pts$x, pts$y, function(x, y) matrix(c(x, y), ncol=2))
ri <- map2(r, po, ws_summarize) %>%
      map_dbl("isotropy") %>%
      round(3)
rd <- map2(r, po, ws_summarize) %>%
      map_dbl("mean_distance") %>%
      round(3)

r <- r %>% lapply(rasterToPoints) %>% lapply(as.data.frame)
for(i in 1:length(r)) r[[i]]$id <- i
for(i in 1:length(r)) r[[i]]$isotropy <- ri[i]
for(i in 1:length(r)) r[[i]]$distance <- rd[i]
r <- bind_rows(r)

ggplot(r %>% mutate(dir = isotropy ^ distance), 
       aes(x, y, fill=wind_fwd)) +
      facet_wrap(~dir, scales="free") +
      geom_raster() +
      scale_fill_viridis_c() +
      theme_void()




