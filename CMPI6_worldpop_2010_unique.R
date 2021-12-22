rm(list=ls())
library(raster)
library(rgdal)
library(gdalUtils)
library(foreach)
library(doMC)

cores <- 60
registerDoMC(cores)
source('world_pop_functions2.R')
source('SDMTools_functions.R')

## make sure run settings file is started at 1 =run
run_file <- "/media/niklas/data2/world_population/run_settings.txt"
rv <- read.table(run_file,head=F)
rv[,2] <- 1
write.table(rv,run_file,col.names = F,row.names = F)

##settings
tmp_folder <- '/media/data3/world_pop/tmp/'
data_dir <- '/media/data3/world_pop/data/' #DATA dir 
out_dir <- '/media/data3/world_pop/output/'
pixels <- 900 ## number of poppixels in each urb frac pixel
ssp_list <- c(  1  , 2, 3, 4, 5, 6)
rcp_list <- c('2_6','4_5','7','3_4','6','8_5')



## LOAD DATA ###

##
world <- readOGR(paste(data_dir,"WorldMap/TM_WORLD_BORDERS03_fixed.shp",sep=""),layer="TM_WORLD_BORDERS03_fixed") #needed to remove island from south africa and Eq Guinea
world.df <- data.frame(world)
world.df$id <- seq(1,length(world.df$FIPS))
pop2010 <- raster(paste(out_dir,'worldpop_2010_matched.grd',sep=''))
dist_urb <- raster(paste(data_dir,'urb_extnew/MGUP_2010_qgis_dist.tif',sep='')) #distance to urban created in QGIS from MGUP data

# country_code <- read.csv(paste(data_dir,'popssps/countryCodes.csv',sep=''))
dist_roads <- raster(paste(out_dir,'dist_to_road_arcHighres_matched.tif',sep='')) #distance to road created in QGIS from road data


scen <- 1
ssp <- ssp_list[scen] 
rcp <- rcp_list[scen]  





urb_frac <- raster(paste(data_dir,'LUH2_v2f_beta/hurtt_urbanfrac_cmip6_SSP',ssp,'_RCP',rcp,'.tif',sep=''))

pop <- pop2010


unique_test <- data.frame(id=NA,maxnum_same=NA,isna=NA)
cellstorun <- cellsFromExtent(urb_frac, extent(urb_frac), expand=FALSE)



rasterOptions(tmpdir = tmp_folder) 
print(paste('MAXcelltorun',max(cellstorun)))
# 1:max(cellstorun)

res.list <- NULL
res.list <- foreach(id=1:max(cellstorun),.inorder=F) %dopar% { 



# for(id in 1:max(cellstorun)){
  
  if(id %% 100 == 0) print(paste(id,'of',max(cellstorun),'--',Sys.time()))
  x <- xFromCell(urb_frac,id) - 0.125
  y <- yFromCell(urb_frac,id) - 0.125
  ext <- extent(x,x+0.25,y,y+0.25)
  a <- intersect(ext,extent(pop)) #check if this should be processed. only if within area chosen
  # pop_w <- NA
  if(!is.null(a)){
    if(a == ext){
      pop_cell <- crop(pop,ext)
      if(!is.null(pop_cell[!is.na(pop_cell)])){ ##only do if we have values
        
        b <- intersect(ext,extent(dist_roads))
        c <- intersect(ext,extent(dist_urb))
        
        if(!is.null(b) & !is.null(c)){
          
          
          
          dist_roads_cell <- crop(dist_roads,ext)
          dist_roads_cell <- 0.00001 + (max(dist_roads_cell[,])- dist_roads_cell) / (max(dist_roads_cell[,])- min(dist_roads_cell[,])) / 10000 #add just above zero-0.1 based on distance to road.
          dist_roads_cell <- resample(dist_roads_cell,pop_cell,method='ngb')  
          
          dist_urb_cell <- crop(dist_urb,ext)
          if(sum(dist_urb_cell[,],na.rm=T)!=0){
            dist_urb_cell <- 0.00001 + (max(dist_urb_cell[,])- dist_urb_cell) / (max(dist_urb_cell[,])- min(dist_urb_cell[,])) / 10000 #add just above zero-0.1 based on distance to road.
          }else{
            dist_urb_cell[,] <- 0
          }
          
          dist_urb_cell <- resample(dist_urb_cell,pop_cell,method='ngb')  
          
          if(sum(pop_cell[!is.na(pop_cell)])>0){
            
            ##center of gravity
            
            cog <- as.numeric(COGravity(pop_cell, y = NULL, z = NULL, wt = NULL))
            cog.spdf <- SpatialPointsDataFrame(coords = data.frame(latitude=cog[1],longitude=cog[3],data=1), data = data.frame(data=1), proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
            cog_dist <- pop_cell;  cog_dist[,] <- NA
            cog_dist[cellFromXY(cog_dist, cog.spdf) ] <- 1
            cog_dist <- raster::distance(cog_dist)
            cog_dist <- 0.00001 + (max(cog_dist[,])- cog_dist) / (max(cog_dist[,])- min(cog_dist[,])) / 10000
            
            
            pop_w <- pop_cell + dist_roads_cell + cog_dist + dist_urb_cell   #need to add otherwise we get pixels with zero that doesnt get any difference.
            
          }else if(!is.null(pop_cell[!is.na(pop_cell)])){ ##only zeros, add some values based on distance to road
            pop_w <- pop_cell + dist_roads_cell  + dist_urb_cell    #need to add otherwise we get pixels with zero that doesnt get any difference.
            
          }
          pop_w
          # res.list <- list(res.list,pop_w)
        }
      }
    }
    
    
    
  }
  
  
}    



raster.list <- list(); pp <- 1
for(i in 1:length(res.list)){ 
  val <- res.list[[i]]
  if(!is.null(val)){
    raster.list[[pp]] <- val  
    pp <- pp+1
  }
  
}


print('MOSAIC in progress')

ptm <- proc.time()
raster.list$fun <- mean
pop_unique <- do.call(merge,raster.list)
print(proc.time()-ptm)
print(length(raster.list))





print('WRITING to file')
writeRaster(pop_unique,paste(out_dir,'CMIP6_worldpop2010_cog_MODurb_road_unique.tif',sep=''),format='GTiff',overwrite=T,NAflag=-9999)
writeRaster(pop_unique,paste(out_dir,'CMIP6_worldpop2010_cog_MODurb_road_unique.grd',sep=''),format='raster',overwrite=T,datatype='FLT4S')
print('DONE')






