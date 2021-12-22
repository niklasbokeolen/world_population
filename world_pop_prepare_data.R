rm(list=ls())
library(raster)
library(rgdal)
library(SDMTools)


data_dir <- 'data/' 
out_dir <- 'output/'


pop2010 <- raster( paste(data_dir,'worldpop/ppp_2010_1km_Aggregated.tif',sep=''))


##Coast/water mask
water.mask <- raster(paste(data_dir,'coastlines/poly2raster_coast_arcmap_reclass.tif',sep='')) ## created with polygon 2 raster in arcmap. cell size set to match pop_2010. value set to FID.
water.mask <- crop(water.mask,pop_2010)

water.mask[!is.na(water.mask)] <- 1

plot(water.mask)
beginCluster(6)
water.mask <- resample(water.mask,pop_2010,method='ngb')
endCluster()
writeRaster(water.mask,paste(out_dir,'coast-mask.tif',sep=''),format='GTiff',overwrite=T,NAflag=-9999)


# worldpop 2010 - global one file

water.mask <- raster(paste(out_dir,'coast-mask.tif',sep=''))

pop2010 <- pop2010*water.mask

writeRaster(pop2010,paste(out_dir,'worldpop_2010_matched.grd',sep=''),format='raster',overwrite=T,datatype='FLT4S')


