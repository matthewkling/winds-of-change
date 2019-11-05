

library(windscape)
library(tidyverse)
library(raster)


cfsr_rose <- function(infile, outdir, ncores, p){
      message(infile)
      
      # file admin
      outfile <- paste0(outdir, "/", basename(infile))
      outfile <- sub("grb2", "tif", outfile)
      if(any(file.exists(c(outfile,
                           sub("\\.tif", paste0("_", ncores, ".tif"), outfile))))){
            return("skipping")}
      
      # load inventory file and identify layers to process
      # (the "analysis" and "6-hour forecast" timesteps are redundant)
      inv <- paste0(infile, ".inv") %>%
            readLines()
      process <- which(!grepl("6 hour fcst", inv))
      
      # generate windrose
      x <- infile %>%
            brick() %>%
            subset(process) %>%
            windrose_rasters(ncores=ncores, 
                             p=p,
                             outfile=outfile)
      return("success")
}


wr_summarize <- function(fd, years=1980:2009, months=1:12, p, outfile){
      
      f <- list.files(fd, full.names=T)
      yr <- substr(f, nchar(f)-11, nchar(f)-8) %>% as.integer()
      mo <- substr(f, nchar(f)-7, nchar(f)-6) %>% as.integer()
      f <- f[yr %in% years & mo %in% months]
      
      f %>% 
            lapply(stack) %>%
            Reduce("+", .) %>%
            "/"(24 * 365/12 * length(months) * length(years)) %>% # number of hours
            "^"(1/p) %>%
            writeRaster(outfile, overwrite=T)
}


# hourly input data
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
#f <- f[grepl("gdas\\.199404", f)]

# velocity
intdir <- "data/windrose/monthly_p1"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=1)
wr_summarize(intdir, p=1, 
             outfile="data/windrose/windrose_p1_1980_2009.tif")

# seasonal versions
wr_summarize(intdir, months=c(12, 1, 2), p=1, 
             outfile="data/windrose/windrose_p1_DJF.tif")
wr_summarize(intdir, months=c(3:5), p=1, 
             outfile="data/windrose/windrose_p1_MAM.tif")
wr_summarize(intdir, months=c(6:8), p=1, 
             outfile="data/windrose/windrose_p1_JJA.tif")
wr_summarize(intdir, months=c(9:11), p=1, 
             outfile="data/windrose/windrose_p1_SON.tif")







# at 1000 hPa
f <- list.files("f:/CFSR/wnd1000", full.names=T)
f <- f[!grepl("inv", f)]
intdir <- "data/windrose/monthly_p1_wnd1000"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=1)
wr_summarize(intdir, p=1, outfile="data/windrose/windrose_p1_wnd1000.tif")

# at 850 hPa
f <- list.files("f:/CFSR/wnd850", full.names=T)
f <- f[!grepl("inv", f)]
intdir <- "data/windrose/monthly_p1_wnd850"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=1)
wr_summarize(intdir, p=1, outfile="data/windrose/windrose_p1_wnd850.tif")











# drag
intdir <- "data/windrose/monthly_p2"
map(f[1:120], possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=2)
wr_summarize(intdir, years=10, p=2, 
             outfile="data/windrose/windrose_p2_2000s.tif")

# force
intdir <- "data/windrose/monthly_p3"
map(f[1:120], possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=3)
wr_summarize(intdir, years=10, p=3, 
             outfile="data/windrose/windrose_p3_2000s.tif")






x2 <- substr(x, nchar(x[1])-11, nchar(x[1])-4)

