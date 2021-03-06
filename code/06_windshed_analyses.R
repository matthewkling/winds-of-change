

library(windscape)
library(tidyverse)
library(raster)
library(gdistance)
library(geosphere)

select <- dplyr::select


source("code/05_windshed_functions.R")

windscapes <- function(data="data/windrose/windrose_p1_wnd10m.tif", 
                       climdata="data/geographic/processed/temperature.tif", # vector of 1+ climate datasets
                       tag="p1_wnd10m",
                       access=list(name = "inv", fx = function(x){1/x}, form = "1/x"),
                       method="gaussian", # c("gaussian", "triangular", "threshold")
                       sigma=2, # vector of 1+ climate sigmas
                       radius=250, # km
                       water=.1, # accessibility multiplier
                       subsample=NULL, # list(n=1000, seed=12345)
                       output="csv", # c("csv", "df")
                       overwrite=F, 
                       ncores = 7){
      
      # admin
      message(tag)
      if(output == "csv"){
            outfile <- paste0("data/windshed/", tag, "_", radius, "km_", access$name, ".csv")
            if(!overwrite & file.exists(outfile)){
                  message(paste(access$name, radius, "output already exists -- aborting"))
                  return(NULL)
            }
      }
      
      # load current/future climate data
      climate <- stack(climdata) %>% unwrap(180)
      
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
      if(!is.null(subsample)){
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




####### main results ##########

windscapes(data = "data/windrose/windrose_p1_wnd10m.tif", 
           climdata = "data/geographic/processed/temperature.tif", 
           tag = "p1_wnd10m",
           access = list(name = "inv", fx = function(x){1/x}, form = "1/x"),
           method = "gaussian",
           sigma = 2,
           radius = 250,
           water = .1,
           subsample = NULL,
           output = "csv",
           overwrite = F, 
           ncores = 7)



############ sensitivity analyses ##############

sensitivity_plot <- function(x, maintag, title, outfile, height, 
                             tlabels = NULL, tlevels = NULL,
                             log=F, method="spearman", frows = 1,
                             comptag=NULL){
      
      xx <- x
      if(is.null(tlevels)) tlevels <- tlabels
      
      d <- x %>%
            mutate(var = overlap_fwd_windshed_size / clim_fwd_windshed_size) %>%
            select(x, y, tag, var) %>%
            spread(tag, var) %>%
            mutate_(main = maintag) %>%
            gather(tag, comp, -main, -x, -y) %>%
            mutate(tag=factor(tag, levels=unique(xx$tag))) %>%
            mutate(weight = cos(y/180*pi))
      
      wcor <- function(x, y, weights) cov.wt(cbind(x, y), wt = weights, cor = TRUE)$cor[1,2]
      
      if(!is.null(tlabels)) d <- mutate(d, tag = factor(tag, levels = tlevels, labels = tlabels))
      
      if(log){
            cm <- d %>%
                  filter(is.finite(main), is.finite(comp)) %>%
                  mutate_at(vars(main, comp), log10) %>%
                  group_by(tag) %>%
                  filter(is.finite(main + comp + weight)) %>%
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
      
      d <- d %>% mutate(mt = tag == tlabels[tlevels == maintag])
      
      p <- ggplot(d, aes(comp, main, color = mt )) +
            facet_grid(. ~ tag, scales="free") +
            geom_point(alpha = 0.2, shape = 16, size = 0.5) +
            geom_text(data=cm, aes(x=x, label=paste("r =", round(r, 2))), 
                      color="red", hjust=.5) +
            scale_color_manual(values=c("black", "red")) +
            theme_minimal() +
            theme(legend.position = "none",
                  axis.text.x=element_text(angle=45, hjust=1)) +
            labs(x="Wind facilitation (1/h)",
                 y=paste0("Wind facilitation\nfor ", tlabels[tlevels == maintag], " (1/h)"))
      
      if(log) p <- p + scale_x_log10() + scale_y_log10()
      if(frows > 1) p <- p + facet_wrap(~tag, nrow = frows)
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
      
      maintag <- sub("-", "", maintag)
      comptag <- sub("-", "", comptag)
      names(d) <- sub("-", "", names(d))
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





##### climate variables
temp <- "data/geographic/processed/temperature.tif"
ppt <- "data/geographic/processed/precipitation.tif"
bio5 <- "data/geographic/processed/bio5.tif"
bio6 <- "data/geographic/processed/bio6.tif"
bio13 <- "data/geographic/processed/bio13.tif"
bio14 <- "data/geographic/processed/bio14.tif"

# sigmas equivalent to temp sigma of 2:
get_sigma <- function(new_raster, temp_raster, temp_sigma){
      temp_sigma / sd(raster(temp_raster)[], na.rm = T) * sd(raster(new_raster)[], na.rm = T)
}
s1 <- 2
s12 <- get_sigma(ppt, temp, s1)
s5 <- get_sigma(bio5, temp, s1)
s6 <- get_sigma(bio6, temp, s1)
s13 <- get_sigma(bio13, temp, s1)
s14 <- get_sigma(bio14, temp, s1)

specs <- list(climdata = list(temp, bio5, bio6, ppt, bio13, bio14,
                              c(temp, ppt), c(bio5, bio6), c(bio13, bio14),
                              c(bio5, bio6, ppt), c(temp, bio13, bio14), c(bio5, bio6, bio13, bio14)),
              sigma = list(s1, s5, s6, s12, s13, s14,  
                           c(s1, s12), c(s5, s6), c(s13, s14),
                           c(s5, s6, s12), c(s1, s13, s14), c(s5, s6, s13, s14)),
              tag = c("temp", "b5", "b6", "prec", "b13", "b14",
                      "temp_prec", "b5_b6", "b13_b14",
                      "b5_b6_prec", "temp_b13_b14", "b5_b6_b13_b14"))
ws <- specs %>% pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="temp", comptag="temp_prec", frows = 2,
                 tlabels = specs$tag, tlevels = specs$tag,
                 title="Sensitivity to climate variables", log=F,
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_vars.png",
                 height=3)
sensitivity_plot(ws, maintag="temp", comptag="temp_prec", frows = 2,
                 tlabels = specs$tag, tlevels = specs$tag,
                 title="Sensitivity to climate variables", log=T,
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_vars_log.png",
                 height=3)
file.copy("figures/windsheds/sensitivity/sensitivity_scatters_vars_log.png",
          "figures/manuscript/SI_fig_sens_vars.png", overwrite=T)


##### conductance function
ws <- list(data=c("data/windrose/windrose_p0_wnd10m.tif",
                  "data/windrose/windrose_p1_wnd10m.tif",
                  "data/windrose/windrose_p2_wnd10m.tif",
                  "data/windrose/windrose_p3_wnd10m.tif"),
           tag=c("constant", "linear", "quadratic", "cubic")) %>%
      pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="linear", comptag="quadratic", 
                 tlevels = c("constant", "linear", "quadratic", "cubic"),
                 tlabels = c("constant", "linear", "quadratic", "cubic"),
                 title="Sensitivity to conductance function",
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_p.png",
                 height=3)
sensitivity_plot(ws, maintag="linear", comptag="quadratic", 
                 tlevels = c("constant", "linear", "quadratic", "cubic"),
                 tlabels = c("constant", "linear", "quadratic", "cubic"),
                 title="Sensitivity to conductance function", log=T,
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_p_log.png",
                 height=3)

curves <- tibble(x = seq(.1, 50, .1),
                 constant = 1,
                 linear = x,
                 quadratic = x^2,
                 cubic = x^3) %>%
      gather(stat, y, -x) %>%
      mutate(stat = factor(stat, levels = c("cubic", "quadratic", "linear", "constant"))) %>%
      
      ggplot(aes(x, y, color = stat)) +
      geom_line() +
      scale_x_log10(breaks = c(.01, .1, 1, 10), labels = c(.01, .1, 1, 10)) +
      scale_y_log10(breaks = c(.001, 1, 1000), labels = c(.001, 1, 1000)) +
      theme_minimal() +
      labs(x = "horizontal windspeed (m/s, log scale)",
           y = "conductance\n(log scale)",
           color = NULL)
library(grid)
img <- png::readPNG("figures/windsheds/sensitivity/sensitivity_scatters_p_log.png") %>%
      rasterGrob(interpolate=TRUE)

p <- (curves | plot_spacer() | plot_spacer()) / img + plot_layout(heights = c(1, 2))

source("E:/edges/range-edges/code/utilities.r")
ggs(paste0("figures/manuscript/SI_fig_sens_p.png"),
    p, width=8, height=5, units="in",
    add = grid.text(letters[1:2], 
                    x=c(.02, .02), 
                    y=c(.97, .53),
                    gp=gpar(fontsize=20, fontface="bold", col="black")))



##### season
ws <- list(data=c("data/windrose/windrose_p1_wnd10m.tif",
                  "data/windrose/windrose_p1_wnd10m_DJF.tif",
                  "data/windrose/windrose_p1_wnd10m_MAM.tif",
                  "data/windrose/windrose_p1_wnd10m_JJA.tif",
                  "data/windrose/windrose_p1_wnd10m_SON.tif"),
           tag=c("annual", "December-February", "March-May", "June-August", "September-November")) %>%
      pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="annual", comptag="December-February",
                 title="Sensitivity to season of year",
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_seasons.png",
                 height=2)
sensitivity_plot(ws, maintag="annual", comptag="December-February",
                 title="Sensitivity to season of year", log=T,
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_seasons_log.png",
                 height=2)
file.copy("figures/windsheds/sensitivity/sensitivity_scatters_seasons_log.png",
          "figures/manuscript/SI_fig_sens_seasons.png", overwrite=T)


##### landscape radius
tags <- c("d50km", "d100km", "d250km", "d500km", "d1000km", "d2500km")
labs <- c("50 km", "100 km", "250 km", "500 km", "1000 km", "2500 km")
ws <- list(data=c("data/windrose/windrose_p1_wnd10m.tif"),
           radius=c(50, 100, 250, 500, 1000, 2500),
           tag=c("d50km", "d100km", "d250km", "d500km", "d1000km", "d2500km")) %>%
      pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="d250km", comptag="d1000km",
                 tlevels = tags, tlabels = labs,
                 title="Sensitivity to landscape size",
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_radius.png",
                 height=2.5)
sensitivity_plot(ws, maintag="d250km", comptag="d1000km",
                 tlevels = tags, tlabels = labs,
                 title="Sensitivity to landscape size", log=T,
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_radius_log.png",
                 height=2.5)
file.copy("figures/windsheds/sensitivity/sensitivity_scatters_radius_log.png",
          "figures/manuscript/SI_fig_sens_radius.png", overwrite=T)


##### atmospheric layer
ws <- list(data=c("data/windrose/windrose_p1_wnd10m.tif",
                  "data/windrose/windrose_p1_wnd1000.tif",
                  "data/windrose/windrose_p1_wnd850.tif",
                  "data/windrose/windrose_p1_wnd700.tif"),
           tag=c("z10m", "z1000hpa", "z850hpa", "z700hpa")) %>%
      pmap_df(windscapes, output="df", subsample=list(n=1000, seed=12345))
sensitivity_plot(ws, maintag="z10m", comptag="z850hpa", title="Sensitivity to altitude", 
                 tlevels = c("z10m", "z1000hpa", "z850hpa", "z700hpa"),
                 tlabels = c("10 m", "1000 hpa", "850 hpa", "700 hpa"),
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_altitude.png",
                 height=2.5)
sensitivity_plot(ws, maintag="z10m", comptag="z850hpa", title="Sensitivity to altitude", log=T, 
                 tlevels = c("z10m", "z1000hpa", "z850hpa", "z700hpa"),
                 tlabels = c("10 m", "1000 hpa", "850 hpa", "700 hpa"),
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_altitude_log.png",
                 height=2.5)
file.copy("figures/windsheds/sensitivity/sensitivity_scatters_altitude_log.png",
          "figures/manuscript/SI_fig_sens_altitude.png", overwrite=T)


##### accessibility function
fx <- list(access=list(list(name = "inv", fx = function(x){1/x}, form = "1/h"),
                       list(name = "sqrtinv", fx = function(x) {sqrt(1/x)}, form = "sqrt(1/h)"),
                       list(name = "exp995", fx = function(x) {.995 ^ x}, form = ".995 ^ h"),
                       list(name = "exp99", fx = function(x) {.99 ^ x}, form = ".99 ^ h"),
                       list(name = "exp95", fx = function(x) {.95 ^ x}, form = ".95 ^ h"),
                       list(name = "invlog", fx = function(x) {1/log(x)}, form = "1/log(h)")),
           tag=c("inv", "sqrtinv", "exp995", "exp99", "exp95", "invlog"))
ws <- pmap_df(fx, windscapes, output="df", subsample=list(n=1000, seed=12345))
flabs <- map_chr(fx$access, "form")
flevs <- map_chr(fx$access, "name")
sensitivity_plot(ws, maintag="inv", comptag="exp95", tlevels = flevs, tlabels = flabs,
                 title="Sensitivity to accessibility function",
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_access.png",
                 height=2)
sensitivity_plot(ws %>% na.omit(), maintag="inv", comptag="exp95", tlevels = flevs, tlabels = flabs,
                 title="Sensitivity to accessibility function", log=T,
                 outfile="figures/windsheds/sensitivity/sensitivity_scatters_access_log.png",
                 height=2)
curves <- map_df(fx$access, 
                 function(fun) data.frame(name=fun$name, fun=fun$form, 
                                          x = 1:1000, y = fun$fx(1:1000))) %>%
      mutate(fun = factor(fun, levels = flabs)) %>%
      ggplot(aes(x, y, color=fun)) + 
      geom_line() +
      theme_minimal() +
      ylim(0, 1) +
      theme(legend.title = element_blank()) +
      labs(x="wind hours",
           y="wind accessibility") +
      coord_fixed(ratio=1000)
ggsave("figures/windsheds/sensitivity/sensitivity_scatters_access_curves.png", curves,
       width=4, height=4, units="in")

library(grid)
img <- png::readPNG("figures/windsheds/sensitivity/sensitivity_scatters_access_log.png") %>%
      rasterGrob(interpolate=TRUE)

p <- (curves | plot_spacer()) / img

source("E:/edges/range-edges/code/utilities.r")
ggs(paste0("figures/manuscript/SI_fig_sens_access.png"),
    p, width=8, height=5, units="in",
    add = grid.text(letters[1:2], 
                    x=c(.02, .02), 
                    y=c(.97, .43),
                    gp=gpar(fontsize=20, fontface="bold", col="black")))



##### analog function
ws <- expand.grid(method=c("gaussian", "triangular", "threshold"),
                  sigma=c(.5, 1, 2, 3, 4)) %>%
      mutate(tag=paste0(method, "_", sigma)) %>%
      pmap_df(windscapes, data="data/windrose/windrose_p1_wnd10m.tif", 
              output="df", subsample=list(n=1000, seed=12345))
maintag="gaussian_2"
title="Sensitivity to climate analog specificity"
outfile="figures/windsheds/sensitivity/sensitivity_scatters_climate.png"
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
      theme(legend.position = "none",
            axis.text.x=element_text(angle=45, hjust=1)) +
      labs(x="Wind facilitation (1/h)",
           y=paste0("Wind facilitation for Gaussian with breadth = 2 (1/h)"))
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
      labs(x="Climatic difference (�C)",
           y="Similarity")
ggsave(outfile %>% sub("\\.png", "_curves.png", .), 
       af, width=8, height=height, units="in")

pw <- af / p
ggs(paste0("figures/manuscript/SI_fig_sens_climate.png"),
    pw, width=8, height=10, units="in",
    add = grid.text(letters[1:2], 
                    x=c(.02, .02), 
                    y=c(.97, .47),
                    gp=gpar(fontsize=20, fontface="bold", col="black")))
