
library(tidyverse)

# wind at 10m
for(year in 1980:2009){
      for(month in str_pad(1:12, 2, "left", 0)){
            url <- paste0("https://nomads.ncdc.noaa.gov/data/cfsr/",
                          year, month, "/wnd10m.gdas.",
                          year, month, ".grb2")
            outfile <- paste0("f:/CFSR/wnd10m/", basename(url))
            if(file.exists(outfile)) next()
            message(paste(year, month))
            download.file(url, outfile, mode="wb")
            download.file(paste0(url, ".inv"), paste0(outfile, ".inv"))
      }
}
stop()

# air temp at 2m
for(year in 1980:1989){
      for(month in str_pad(1:12, 2, "left", 0)){
            url <- paste0("https://nomads.ncdc.noaa.gov/data/cfsr/",
                          year, month, "/tmp2m.gdas.",
                          year, month, ".grb2")
            outfile <- paste0("f:/CFSR/tmp2m/", basename(url))
            if(file.exists(outfile)) next()
            message(paste(year, month))
            download.file(url, outfile, mode="wb")
            download.file(paste0(url, ".inv"), paste0(outfile, ".inv"))
      }
}

# surface temperature