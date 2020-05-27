
# reformat the raw dataset, which is too large for github, contains many extraneous variables, and has obtuse variable names

library(tidyverse)

d <- read_csv("results/figure_4/p1_wnd10m_250km_inv.csv") %>%
      select(x, y,
             wind_fwd_windshed_size, wind_rev_windshed_size,
             clim_fwd_windshed_size, clim_rev_windshed_size,
             overlap_fwd_windshed_size, overlap_rev_windshed_size) %>%
      mutate_at(vars(-x, -y), signif, digits = 5) %>%
      rename_at(vars(contains("windshed")), function(x) str_replace(x, "fwd_windshed_size", "outbound")) %>%
      rename_at(vars(contains("windshed")), function(x) str_replace(x, "rev_windshed_size", "inbound")) %>%
      write_csv("results/figure_4/windscape_data.csv")
