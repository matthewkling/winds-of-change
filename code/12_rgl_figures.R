


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

select <- dplyr::select

infile <- "data/windshed/p1_30y_250km_inv.csv"

source("E:/edges/range-edges/code/utilities.r")

truncate <- function(x, q=.005, sides=c("high")){
      q <- quantile(x, c(q, 1-q), na.rm=T)
      if("low" %in% sides) x[x<q[1]] <- q[1]
      if("high" %in% sides) x[x>q[2]] <- q[2]
      x
}


### 3d surfaces with RGL 
rglplots <- function(){
      library(rgl)
      
      ext <- extent(c(-160, -70, 0, 80))
      
      elev <- raster("F:/chelsa/elevation/mn30_grd/w001001.adf") %>% 
            crop(ext) %>%
            aggregate(3) %>%
            focal(matrix(1, 5, 5), mean, na.rm=T)
      
      f <- read_csv(infile) %>%
            select(-runtime) %>%
            gather(var, value, -x, -y) %>%
            separate(var, c("property", "direction", "moment", "stat"), sep="_") %>%
            mutate(direction = ifelse(direction=="fwd", "outbound", "inbound"))
      
      dem <- stack("data/geographic/processed/elevation.tif") %>% crop(ext)
      
      for(drn in unique(f$direction)){
            
            d <- f %>%
                  filter(direction == drn,
                         moment == "windshed",
                         stat == "size") %>%
                  spread(property, value) %>%
                  mutate(windfill = overlap / clim)
            
            
            library(geosphere)
            dg <- function(lat) distGeo(c(0, lat), c(1, lat))/1000
            areas <- d %>% select(y) %>% distinct()
            areas$area <- sapply(areas$y, dg)
            
            
            dd <- f %>%
                  filter(direction == drn,
                         moment == "windshed",
                         property == "wind",
                         stat == "isotropy") %>%
                  rename(isotropy = value) %>%
                  select(x, y, isotropy) %>%
                  left_join(d, .) %>%
                  
                  filter(is.finite(windfill),
                         is.finite(clim),
                         is.finite(isotropy)) %>%
                  mutate(isotropy=truncate(isotropy),
                         windfill = truncate(windfill, .001, c("high", "low"))) %>%
                  
                  left_join(areas) %>%
                  arrange(clim) %>%
                  mutate(climrnk = cumsum(area)) %>%
                  arrange(windfill) %>%
                  mutate(windfillrnk = cumsum(area)) %>%
                  arrange(isotropy) %>%
                  mutate(isotropyrnk = cumsum(area))
            
            
            for(stat in c("windfill", "syndrome", "climwind")){
                  
                  if(stat=="windfill"){
                        ramp <- colorRampPalette(c("red", "yellow", "limegreen"))
                        dd <- dd %>% arrange(windfillrnk) %>%
                              mutate(color = ramp(nrow(.)))
                  }
                  
                  if(stat=="syndrome"){
                        palette <- c("white", "red", "gold", "black")
                        dd <- dd %>% mutate(syndrome = case_when(climrnk < mean(climrnk) ~ "climate-limited",
                                                                 windfillrnk > mean(windfillrnk) ~ "wind-facilitated",
                                                                 isotropyrnk > mean(isotropyrnk) ~ "speed-hindered",
                                                                 TRUE ~ "direction-hindered")) %>% 
                              mutate(color = colors2d(cbind(.$isotropyrnk, .$windfillrnk),
                                                      palette[c(4,3,2,4)]))
                        
                        
                  }
                  
                  if(stat=="climwind"){
                        pal <- c("forestgreen", "yellow", "red", "darkblue")
                        dd <- dd %>% mutate(color = colors2d(cbind(.$climrnk, .$windfillrnk), pal))
                        
                  }
                  
                  cst <- dd %>%
                        filter(x >= ext@xmin, x <= ext@xmax, y >= ext@ymin, y <= ext@ymax)
                  cst <- expand.grid(x=unique(cst$x), y=unique(cst$y)) %>%
                        left_join(cst) %>%
                        mutate(color=ifelse(is.na(color), "white", color)) %>%
                        arrange(-y, x)
                  
                  crgb <- col2rgb(cst$color)
                  if(stat=="syndrome"){
                        crgb <- rbind(crgb, cst$climrnk/max(cst$climrnk, na.rm=T)) %>%
                              apply(2, function(x){if(any(is.na(x))) return(rep(255, 3))
                                    return(255 - ((255 - x[1:3]) * x[4]))})
                  }
                  
                  clr <- stack(dem, dem, dem)
                  clr[[1]][] <- crgb[1,]
                  clr[[2]][] <- crgb[2,]
                  clr[[3]][] <- crgb[3,]
                  clr <- resample(clr, elev)
                  clr <- trim(clr) %>% reclassify(c(-Inf, 0, 0, 255, Inf, 255))
                  elev <- crop(elev, clr)
                  
                  clr <- rgb(values(clr), maxColorValue = 255)
                  wind_color <- matrix(clr, nrow=nrow(elev), byrow=T) %>% t()
                  demm <- matrix(values(elev), nrow=nrow(elev), byrow=T) %>% t()
                  wind_color[demm<=0] <- "white"
                  
                  z <- demm
                  x <- -1000 * (1:nrow(z))
                  y <- 1000 * (1:ncol(z))
                  
                  open3d()
                  rgl.surface(x, -y, z*12, col=wind_color, shininess=128)
                  rgl.viewpoint(theta = 180, phi = 45, zoom=.17)
                  rgl.clear("lights")
                  light3d(phi=15, specular="gray60")
                  par3d("windowRect" = c(0,0,1000,1000))
                  snapshot3d(paste0("figures/windsheds/rgl/", stat, "_", drn, ".png"))
                  
            }
      }
      
      rgl.quit()
      
}
rglplots()
