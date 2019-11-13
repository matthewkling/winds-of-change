

library(windscape)
library(tidyverse)
library(raster)


cfsr_rose <- function(infile, outdir, ncores, p){
   message(infile)
   
   # file admin
   if(!dir.exists(outdir)) dir.create(outdir)
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
   if(length(f) %% (length(years) * length(months)) != 0) stop("non-factorial data")
   
   f %>% 
      lapply(stack) %>%
      Reduce("+", .) %>%
      "/"(24 * 365/12 * length(months) * length(years)) %>% # number of hours
      "^"(1/p) %>%
      writeRaster(outfile, overwrite=T)
}

########################################################



# main analysis
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
intdir <- "data/windrose/monthly_p1_wnd10m"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=1)
wr_summarize(intdir, p=1, outfile="data/windrose/windrose_p1_wnd10m.tif")



# seasonal versions
wr_summarize(intdir, months=c(12, 1, 2), p=1, 
             outfile="data/windrose/windrose_p1_wnd10m_DJF.tif")
wr_summarize(intdir, months=c(3:5), p=1, 
             outfile="data/windrose/windrose_p1_wnd10m_MAM.tif")
wr_summarize(intdir, months=c(6:8), p=1, 
             outfile="data/windrose/windrose_p1_wnd10m_JJA.tif")
wr_summarize(intdir, months=c(9:11), p=1, 
             outfile="data/windrose/windrose_p1_wnd10m_SON.tif")

# p0
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
intdir <- "data/windrose/monthly_p0_wnd10m"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=0)
wr_summarize(intdir, p=0, outfile="data/windrose/windrose_p0_wnd10m.tif")

# p2
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
intdir <- "data/windrose/monthly_p2_wnd10m"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=2)
wr_summarize(intdir, p=2, outfile="data/windrose/windrose_p2_wnd10m.tif")

# p3
f <- list.files("f:/CFSR/wnd10m", full.names=T)
f <- f[!grepl("inv", f)]
intdir <- "data/windrose/monthly_p3_wnd10m"
map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=3)
wr_summarize(intdir, p=3, outfile="data/windrose/windrose_p3_wnd10m.tif")

# multiple atmospheric heights
for(height in c("wnd1000", "wnd850", "wnd700", "wnd500")){
   outfile <- paste0("data/windrose/windrose_p1_", height, ".tif")
   if(file.exists(outfile)) next()
   
   f <- list.files(paste0("f:/CFSR/", height), full.names=T)
   f <- f[!grepl("inv", f)]
   intdir <- paste0("data/windrose/monthly_p1_", height)
   dir.create(intdir)
   map(f, possibly(cfsr_rose, NULL), outdir=intdir, ncores=6, p=1)
   wr_summarize(intdir, p=1, outfile=outfile)
   
   stack(outfile) %>%
      resample(raster("data/windrose/windrose_p1_wnd10m.tif")) %>%
      writeRaster(outfile, overwrite=T)
}

