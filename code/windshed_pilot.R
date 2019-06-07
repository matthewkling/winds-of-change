
library(gdistance)
library(ggplot2)
library(viridis)
library(tidyverse)

# TODO
# investigate possibility of using product rather than sum for accumulating cost distance


ww <- stack("E:/flow/dispersal/data/derived/04_08/windrose_force.tif")
#w <- stack("E:/flow/dispersal/data/derived/full_year/windrose.tif")

ext <- extent(3675098, 8808516, 2130760, 5394390) # us48
#pd <- data.frame(x=6510000, y=3500000)
#pt <- SpatialPoints(pd)

#ext <- extent(3059518, 4158675, 3036942, 3897151) # pacific ocean
#pd <- data.frame(x=3600000, y=3400000)
#pt <- SpatialPoints(pd)


# temperature data, transferred to wind grid
tmean <- raster("F:/chelsa/bio19/CHELSA_bio10_1.tif") %>%
      crop(extent(-126.9978,-110.282,32.37925,46.68744))
dtemp <- tmean %>%
      rasterToPoints() %>%
      as.data.frame()
coordinates(dtemp) <- c("x", "y")
crs(dtemp) <- crs(tmean)
dtemp <- spTransform(dtemp, crs(ww))

temp <- ww[[1]]
temp[] <- 1:ncell(temp)

temp <- rasterize(dtemp, temp,field="CHELSA_bio10_1", fun=mean, na.rm=T) %>%
      trim()
temp <- temp/10

#ext <- extent(4059518, 5158675, 3036942, 3897151) # cali nevada
pops <- list(site1=data.frame(x=4280000, y=3395000),
             site2=data.frame(x=4530000, y=3255000),
             site3=data.frame(x=4480000, y=3495000),
             site4=data.frame(x=4330000, y=3595000))

df <- data.frame()
for(pop in names(pops)){
      message(pop)
      
      pd <- pops[[pop]]
      pt <- SpatialPoints(pd)
      
      w <- crop(ww, ext)
      
      # add coordinates metadata as layers
      rows <- cols <- w[[1]]
      rows[] <- rep(1:nrow(rows), each=ncol(rows))
      cols[] <- rep(1:ncol(rows), nrow(rows))
      w <- stack(w, rows, cols)
      names(w) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
      
      # dispersal function to compute transition likelihood
      disperse <- function(x){
            p <- x[c(1,3,5,7,9,11,13,15)] # SW ...cw... S
            g <- x[17:20]
            if(g[1]<g[2] & g[4]<g[3]) return(p[1]/sqrt(2)/sqrt(2)) #SW
            if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
            if(g[2]<g[1] & g[4]<g[3]) return(p[3]/sqrt(2)/sqrt(2)) #NW
            if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
            if(g[2]<g[1] & g[3]<g[4]) return(p[5]/sqrt(2)/sqrt(2)) #NE
            if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
            if(g[1]<g[2] & g[3]<g[4]) return(p[7]/sqrt(2)/sqrt(2)) #SE
            if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
            # dividing by 2 root 2 prevents the geocorrection from distorting the probability field
      }
      
      
      transition_stack <- function(x, transitionFunction, directions, symm, ...){
            
            brk <- x
            x <- brk[[1]]
            
            tr <- new("TransitionLayer",
                      nrows=as.integer(nrow(x)),
                      ncols=as.integer(ncol(x)),
                      extent=extent(x),
                      crs=projection(x, asText=FALSE),
                      transitionMatrix = Matrix(0,ncell(x),ncell(x)),
                      transitionCells = 1:ncell(x))
            transitionMatr <- transitionMatrix(tr)
            Cells <- which(!is.na(getValues(x)))
            adj <- adjacent(x, cells=Cells, pairs=TRUE, target=Cells, directions=directions)
            
            ##### start modificiations #####
            # format raster data layers to feed to transitionFunction
            # col order is x[[1]][from, to] ... x[[n]][from,to]
            dataVals <- lapply(1:nlayers(brk), 
                               function(i) cbind(values(brk[[i]])[adj[,1]],
                                                 values(brk[[i]])[adj[,2]]))
            dataVals <- do.call("cbind", dataVals)
            ##### end modifications #####
            
            transition.values <- apply(dataVals, 1, transitionFunction, ...)
            
            #transition.values <- log(transition.values)
            #transition.values <- ecdf(transition.values)(transition.values)
            #plot(hist(transition.values[transition.values<10000]))
            transitionMatr[adj] <- as.vector(transition.values)
            transitionMatrix(tr) <- transitionMatr
            matrixValues(tr) <- "resistance"
            return(tr)
      }
      
      trans <- transition_stack(w, disperse, directions=8, symm=F)
      trans <- geoCorrection(trans, type="c")
      cost <- accCost(trans, pt)
      #plot(cost, colNA="black")
      #points(pt)
      
      cost_rev <- cost
      cost_rev[] <- costDistance(trans, coordinates(cost), pt)
      
      
      md <- map_data("state")
      coordinates(md) <- c("long", "lat")
      crs(md) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
      md <- spTransform(md, crs(w))
      md <- as.data.frame(md)
      
      cd <- as.data.frame(rasterToPoints(cost))
      wd <- as.data.frame(rasterToPoints(w))
      resolution <- mean(res(w))/5
      
      
      
      ### climate analogs ###
      
      for(direction in c("forward", "reverse")){
            
            cd <- cost
            delta <- 3
            
            if(direction=="reverse"){
                  cd <- cost_rev
                  delta <- -delta
            }
            
            target <- raster::extract(temp, pt) - delta
            diffs <- abs(temp - target)
            
            kernel <- function(x) exp(-(x^2))
            
            s <- stack(cd %>% crop(diffs) %>% mask(diffs), 
                       diffs %>% crop(cd)) %>%
                  trim() %>%
                  rasterToPoints() %>%
                  as.data.frame() %>%
                  rename(wind=layer.1,
                         suit=layer.2) %>%
                  filter(is.finite(suit), is.finite(wind),
                         suit != 0, wind != 0) %>%
                  mutate(suit = kernel(suit),
                         wind = .05^(wind*1e9),
                         overlap = suit * wind,
                         overlap_total = sum(suit * wind),
                         site = pop,
                         direction = direction)
            
            ggplot(s, aes(x, y, fill=overlap_total )) + 
                  geom_raster() +
                  scale_fill_gradientn(colors=c("black", "green"))
            
            df <- bind_rows(df, s)
      }
}

bbox <- c(xmin=min(s$x)-100000, xmax=5e6,
          ymin=min(s$y)+250000, ymax=max(s$y)-500000)

dfp <- df %>%
      gather(var, value, wind, suit, overlap) %>%
      mutate(var = factor(var, 
                          levels=c("wind", "suit", "overlap"),
                          labels=c("accessibility by wind", 
                                   "future suitability", 
                                   "suitable & accessible"))) %>%
      filter(is.finite(value),
             x>bbox[1], x<bbox[2], y>bbox[3], y<bbox[4],) %>%
      group_by(direction, var) %>%
      mutate(value=scales::rescale(value))

origins <- dfp %>%
      count(x, y)

p <- ggplot() +
      geom_raster(data=origins, aes(x, y), fill="red") +
      geom_raster(data=dfp, aes(x, y, fill=value)) +
      geom_path(data=md, aes(long, lat, group=group), 
                color="white", size=.5, alpha=.2) +
      facet_grid(var~direction+site) +
      scale_fill_gradientn(
            #colours=c("gray98", "darkred")
            #colours=c("black", "red", "yellow")
            colours=c("gray95", "black")
            #colours=c("black", "blue", "cyan")
      ) +
      theme(rect=element_blank(), line=element_blank(),
            axis.text=element_blank(), axis.title=element_blank(),
            legend.position="none") +
      coord_cartesian(xlim=range(dfp$x), ylim=range(dfp$y), 
                      expand=0)
ggsave("figures/windsheds/overlap.png", p, width=10, height=5, units="in")



do <- df %>%
      filter(x>bbox[1], x<bbox[2], y>bbox[3], y<bbox[4],) %>%
      select(direction, site, overlap_total) %>%
      distinct() %>%
      spread(direction, overlap_total)
p <- ggplot(do, aes(forward, reverse, label=site)) +
      geom_hline(yintercept=0) + geom_vline(xintercept=0) +
      geom_point() +
      geom_text(hjust=1, nudge_x=0) +
      theme_minimal() +
      labs(x = "emigration facilitation",
           y = "immigration facilitation")
ggsave("figures/windsheds/overlap_forward_reverse_scatter.png", p, width=6, height=6, units="in")



tempdf <- temp %>%
      rasterToPoints() %>%
      as.data.frame()
p <- ggplot() +
      geom_raster(data=tempdf, aes(x, y, fill=layer)) +
      geom_path(data=md, aes(long, lat, group=group), 
                color="white", size=.5) +
      scale_fill_gradientn(colours=c("black", "blue", "red", "yellow")) +
      coord_cartesian(xlim=range(s$x), ylim=range(s$y)+c(500000, -500000), 
                      expand=0) +
      theme_void() +
      labs(fill="temperature")
ggsave("figures/windsheds/overlap_temp.png", p, width=7, height=4, units="in")
