

# see whether alignment is associated with treeline advance

library(tidyverse)


tls <- tl <- readxl::read_xlsx("data/treeline/treelines.xlsx")
coordinates(tls) <- c("Long °", "Lat °")

stop("run headwinds_02_analysis.r up midway through the first loop")

d <- cbind(tl, extract(s, tls)) %>%
      mutate(theta = abs(wa - ta),
             theta = ifelse(theta > 180, 360-theta, theta))

names(d) <- gsub(" ", "_", names(d))

ggplot(d %>% filter(Treeline_Type == "latitudinal"), 
       aes(Treeline_Form, cos(theta/180*pi), color=Advance)) +
      geom_boxplot()


md <- d %>% 
      #filter(Treeline_Type=="latitudinal") %>%
      filter(!is.na(Treeline_Form), Treeline_Form != "abrupt") %>%
      mutate(tailwind = cos(theta/180*pi),
             adv = as.integer(factor(Advance)))
lm(adv ~ tailwind + Treeline_Form, data=md) %>% summary()

ggplot(md, aes(tailwind, adv, color=paste(Treeline_Type, Treeline_Form))) + 
      geom_smooth() + coord_cartesian(ylim=1:2)


md <- d %>% 
      filter(Treeline_Type=="latitudinal") %>%
      filter(!is.na(Treeline_Form), Treeline_Form != "abrupt") %>%
      filter(Treeline_Form=="krummholz") %>%
      mutate(tailwind = cos(theta/180*pi),
             adv = as.integer(factor(Advance)))
lm(adv ~ tailwind, data=md) %>% summary()




ggplot(md, aes(Advance, tailwind, color=Advance)) + geom_boxplot()

ggplot(md, aes(tailwind, color=Advance)) + geom_density()




f <- filter(d, Treeline_Type=="latitudinal")
t.test(f$theta[f$Advance=="yes"], f$theta[f$Advance=="no"])


