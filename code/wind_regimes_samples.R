
library(windscape)
library(tidyverse)
library(raster)



samples <- function(infile, outdir, points, ncores=1){
      
      # file admin
      message(infile)
      outfile <- paste0(outdir, "/", basename(infile))
      outfile <- sub("grb2", "rds", outfile)
      if(any(file.exists(c(outfile,
                           sub("\\.rds", paste0("_", ncores, ".rds"), outfile))))){
            return("skipping")}
      
      # load and organize hourly wind data
      inv <- paste0(infile, ".inv") %>%
            readLines()
      process <- which(!grepl("6 hour fcst", inv))
      w <- infile %>%
            brick() %>%
            subset(process)
      even <- function(x) x %% 2 == 0
      w <- w[[c(which(!even(1:nlayers(w))),
                which(even(1:nlayers(w))))]]
      
      p <- extract(w, points)
      
      uv <- p %>%
            split(1:nrow(.)) %>%
            lapply(function(x) matrix(x, ncol=2, byrow=F))
      for(i in 1:length(uv)){
            x <- as.data.frame(uv[[i]])
            names(x) <- c("u", "v")
            x$x <- points$x[i]
            x$y <- points$y[i]
            uv[[i]] <- x
      }
      
      saveRDS(uv, outfile)
      
}


# hourly input data
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]

set.seed(123)
land <- raster("f:/cfsr/land.tif") %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(y > -60) %>%
      sample_n(1000) %>%
      dplyr::select(x, y)

map(f, samples, outdir="data/wind_regime/samples", points=land, ncores=6)




###



d <- list.files("data/wind_regime/samples", full.names=T) %>%
   lapply(readRDS) %>% 
   map(bind_rows) %>% bind_rows() %>%
   mutate(id = as.integer(factor(paste(x, y))))




regime <- function(u, v, return="speed"){
   uv <- cbind(u, v)
   colnames(uv) <- NULL
   u <- mean(uv[,1])
   v <- mean(uv[,2])
   
   speeds <- sqrt(uv[,1]^2 + uv[,2]^2)
   speed <- mean(speeds)
   if(return=="speed") return(speed)
   
   dir <- do.call(atan2, as.list(colMeans(uv))) / pi * 180
   
   dirs <- atan2(uv[,1], uv[,2])
   x <- sum(sin(dirs) * speeds) / sum(speeds)
   y <- sum(cos(dirs) * speeds) / sum(speeds)
   
   rbar <- sqrt(x^2 + y^2) # mean resultant vector
   iso <- sqrt(1-rbar) # circular standard deviation
   
   if(return=="anisotropy") return(1 - iso)
   
   #c(speed = speed, dir = dir, 
   #  rbar = rbar, aniso = 1 - iso,
   #  u = u, v = v, x = x, y = y)
}


r <- d %>%
   group_by(id) %>%
   summarize(speed = regime(u, v, "speed"),
             aniso = regime(u, v, "anisotropy"),
             x = x[1],
             y = y[1])






# weak, aniso: 601 (or 85, maybe)
# weak iso: 562
# strong, aniso: 668
# strong, iso: 818


ex <- c(601, 668, 818, 562)
ed <- d %>%
   filter(id %in% ex) %>%
   mutate(id = factor(id, levels = ex))
centroids <- ed %>%
   group_by(id) %>%
   summarize(u = mean(u), v = mean(v))


rad <- 5
circle <- data.frame(angle=seq(0, pi*2, .01)) %>%
   mutate(x=cos(angle) * rad, 
          y=sin(angle) * rad)
acol <- "black"
lim <- 7.5


p <- ed %>% filter(abs(u) < lim, abs(v) < lim) %>%
   ggplot(aes(u, v)) +
   facet_grid(.~id, space="free", scales="free") +
   geom_point(aes(color=id), alpha=.1, shape=16, size=1) +
   geom_path(data=circle, aes(x, y), color=acol) +
   geom_path(data=circle, aes(x/2, y/2), color=acol) +
   annotate(geom="segment", x=0, xend=0, y=-rad*.1, yend=rad * 1.1, color=acol) +
   annotate(geom="segment", x=-rad*.1, xend=rad*.1, y=0, yend=0, color=acol) +
   annotate(geom="text", x=0, y=rad * 1.2, label="N", color=acol, fontface="bold") +
   annotate(geom="text", x=0, y=-rad * 1.1, label=paste0(rad, " m/s"), color=acol) +
   annotate(geom="text", x=0, y=-rad * .6, label=paste0(rad/2, " m/s"), color=acol) +
   geom_point(data=centroids, shape=21, size=3, color="black", fill="white") +
   scale_color_manual(values=c("blue", "red", "gold", "cyan")) +
   theme_void() +
   theme(strip.text = element_blank(), legend.position="none")
ggsave(paste0("figures/regimes/examples_final.png"), 
       p, width=16, height=4.65, units="in", dpi=1200)



########


stop("slow examples below")

s <- r %>%
   mutate(spd = scales::rescale(sqrt(speed)),
          ani = scales::rescale(log10(aniso)),
          pos = scales::rescale(spd + ani),
          neg = scales::rescale(spd - ani))

for(var in c("pos", "neg")){
   for(ext in c("high", "low")){
      if(ext=="high") si <- s %>% mutate(temp = s[[var]]) %>% arrange(temp)  %>% slice(1:25)
      if(ext=="low") si <- s %>% mutate(temp = s[[var]]) %>% arrange(desc(temp)) %>% slice(1:25)
      
      p <- d %>%
         filter(id %in% si$id) %>%
         ggplot(aes(u, v)) +
         facet_wrap(~id, scales="free") +
         geom_hline(yintercept = 0) +
         geom_vline(xintercept = 0) +
         geom_point(alpha=.1, shape=16, size=.2) +
         geom_path(data=circle, aes(x, y), color="red") +
         #coord_fixed() +
         theme_minimal()
      ggsave(paste0("figures/regimes/examples_", var, "_", ext, ".png"), 
             p, width=12, height=12, units="in")
   }
}




rad <- 3
circle <- data.frame(angle=seq(0, pi*2, .01)) %>%
   mutate(x=cos(angle) * rad, 
          y=sin(angle) * rad)

lim <- 7

for(i in seq(0, 900, 100)){
   p <- d %>%
      filter(id %in% (i + 1:100)) %>%
      ggplot(aes(u, v)) +
      facet_wrap(~id) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(alpha=.1, shape=16, size=.1) +
      geom_path(data=circle, aes(x, y), color="red") +
      coord_fixed() +
      theme_minimal() +
      xlim(-lim, lim) + 
      ylim(-lim, lim)
   ggsave(paste0("figures/regimes/examples_", i, ".png"), 
          p, width=12, height=12, units="in")
}



