
library(tidyverse)
library(grid)
library(circular)

n <- 1000
mid <- pi/4
spread <- 1.5
speed_spread <- .5

circ_diff <- function(x, y){
      d <- abs(x - y)
      d[d > pi] <- 2*pi - d[d > pi]
      d
}

ref <- data.frame(angle = rvonmises(n, mid, spread) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid, speed_spread) * runif(n),
             speed = speed/mean(speed),
             panel = "ref")

stren <- data.frame(angle = rvonmises(n, mid, spread) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid, speed_spread) * runif(n),
             speed = speed/mean(speed) * .5,
             panel = "stren")

iso <- data.frame(angle = rvonmises(n, mid, spread * 4) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid, speed_spread * 4) * runif(n) / 2,
             speed = speed/mean(speed),
             panel = "iso")

dir <- data.frame(angle = rvonmises(n, mid + 1.5, spread) %>% as.numeric()) %>%
      mutate(speed = dvonmises(angle, mid + 1.5, speed_spread) * runif(n),
             speed = speed/mean(speed),
             panel = "dir")


d <- bind_rows(ref, dir, stren, iso) %>%
      mutate(panel = factor(panel, levels=c("ref", "dir", "iso", "stren"))) %>%
      mutate(y = cos(angle) * speed,
             x = sin(angle) * speed)


mx <- max(d$speed) * .5

circle <- data.frame(angle=seq(0, pi*2, .01)) %>%
      mutate(x=cos(angle) * mx, 
             y=sin(angle) * mx)

txt <- data.frame(panel = factor(c("ref", "dir", "iso", "stren"), 
                                 levels=c("ref", "dir", "iso", "stren")),
                  x = 0, y = -1.25 * mx,
                  label = c("a local\nwind regime", "direction\ndifference", 
                            "directionality\ndifference", "strength\ndifference"))

brk <- seq(-1, 1, .1)

eb <- element_blank()

p <- ggplot(d) +
      facet_grid(.~panel) +
      geom_path(data=circle, aes(x, y), color="red") +
      #annotate(geom="segment", x=0, xend=c(0,1,0,-1), y=0, yend=c(1,0,-1,0), color="red") +
      geom_segment(aes(x=x*.999, xend=x, y=y*.999, yend=y),
                   arrow = arrow(type="closed", angle=10, length=unit(.1, "in")),
                   alpha = .25) +
      annotate(geom="segment", x=0, xend=0, y=0, yend=.7 * mx, color="red") +
      annotate(geom="text", x=0, y=.85 * mx, label="N", color="red", fontface="bold") +
      annotate(geom="point", x=0, y=0, color="red", size=3) +
      geom_text(data=txt, aes(x, y, label=label), 
                color="red", fontface="bold", lineheight=.75) +
      scale_x_continuous(expand=c(0,0), breaks=brk) +
      scale_y_continuous(expand=c(0,0), breaks=brk) +
      coord_fixed() +
      theme(axis.title=eb, axis.text=eb, axis.ticks=eb,
            panel.background = eb, strip.text=eb, strip.background = eb)
saveRDS(p, "figures/tailwinds/regime_cartoon_plot.rds")
