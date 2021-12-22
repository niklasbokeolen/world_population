rm(list=ls())
setwd('/home/niklas/Dropbox/doktorand/Matlabcodes/world_population/')
library(raster)
library(foreach)
library(doMC)
library(snow)
library(rgdal)
cores <- 1
registerDoMC(cores)

## FUNCTIONS
source('world_pop_functions3.R')
scen <- 1

data_dir <- '/media/niklas/data2/world_population/data/' #DATA dir 
out_dir <- '/media/niklas/data2/world_population/output/' #output dir


world <- readOGR(paste(data_dir,"WorldMap/",sep=""),layer="TM_WORLD_BORDERS03_fixed") #needed to remove island from south africa and Eq Guinea

world.df <- data.frame(world)
world.df$id <- seq(1,length(world.df$FIPS))



country_code <- read.csv(paste(data_dir,'popssps/countryCodes_newiso3.csv',sep=''))
ssp_list <- c(  1  , 2   , 3 , 4   , 4 , 5)
rcp_list <- c('2_6','4_5','7','3_4','6','8_5')

#scen <- 1
ssp <- ssp_list[scen] 
rcp <- rcp_list[scen]  
pop <- read.csv(paste(data_dir,'popssps/popSSP',ssp,'.csv',sep=''),header=TRUE,row.names=1)
colnames(pop) <- c(1:length(pop[1,]))
urban<- read.csv(paste(data_dir,'popssps/urbanShareSSP',ssp,'.csv',sep=''),header=TRUE,row.names=1)
colnames(urban) <- c(1:length(urban[1,]))



##creat output folders
folder <- paste(out_dir,'SSP',ssp,'_RCP',rcp,sep='')
if(!file.exists(folder)) {
  system(paste('mkdir',folder))
}




urb_frac.s <- stack(paste(data_dir,'LUH2_v2f_beta/hurtt_urbanfrac_cmip6_SSP',ssp,'_RCP',rcp,'.tif',sep=''))
pixels <- 900 ## number of pop pixels in each urb_frac


for (year in 2010:2100){
  print(paste(year,'---',Sys.time()))
  outfile_main <- paste(folder,'/','CMPI6_grid_pop_count',year,'_SSP',ssp,'_RCP',rcp,'.tif',sep='')
  

  if(!file.exists(outfile_main)){
    
    if(year==2010){
      pop_grid <-  raster(paste(out_dir,'CMIP6_worldpop2010_cog_MODurb_road_unique.grd',sep=''))
                        
    }else{
      pop_grid <- raster(paste(folder,'/','CMPI6_grid_pop_count',year-1,'_SSP',ssp,'_RCP',rcp,'.tif',sep=''))
    }
    
    
    ##Load Hurtt urban fraction 
    urb_frac <- urb_frac.s[[year-2009]]
    
    ##URBAN MASK #####
    urban_mask <- to_urban_mask_parallell(urb_frac,pop_grid,pixels,cores)  
   
    
    ##############
    
    
    c_pop_grid.list <- foreach(countryid=1:length(world.df$ISO3),.inorder=FALSE) %dopar% {   
     
      country <-  world[world.df$id==countryid,]
      a <- intersect(extent(urb_frac.s[[1]]),extent(country)) #check if country should be processed. only if within area chosen
      c <- FALSE
      if(!is.null(a)) c <- round(a) == round(extent(country)) #if the interesct make sure that the entire country is present in data
      nr <- as.numeric(country_code$nr[as.character(country_code$iso)==as.character(world.df$ISO3[countryid])])
      if(length(nr) > 0 && c){
        
        c_pop <- pop[year-1999,paste(nr)] #country population total from ssp file
        
        c_pop_urb <- c_pop*urban[year-1999,paste(nr)]/100 #number of people in urban
        c_pop_rur <- c_pop - c_pop_urb
        
        
        c_pop_grid <- mask(crop(pop_grid,country),country)
        c_urb_mask <- mask(crop(urban_mask,c_pop_grid),country)
        
        if(c_pop_urb==0) c_urb_mask[c_urb_mask==1] <- 0 # in the cases where urban_frac values are missing distirbute for all in same way
        
        
        c_pop_urb_grid <- overlay(c_pop_grid, c_urb_mask, fun=function(x,y){return(x*y)}) #c_pop_grid * c_urb_mask
        
        if(sum(c_pop_urb_grid[,],na.rm=T)>0){  #make sure we have urban pixels.
          
          
          #gridded urban pop
          
          c_pop_urb_grid <- setValues(c_pop_urb_grid,c_pop_urb_grid[,]* c_pop_urb / sum(c_pop_urb_grid[,],na.rm=T ))
          
          
        }else{ # end if we have urban pixels
          if(c_pop_urb>0) c_pop_rur <- c_pop_rur + c_pop_urb #Special case if we do not have urban pixels but have urban population put that into rural.
          
          c_pop_urb_grid <- c_pop_grid; c_pop_urb_grid[!is.na(c_pop_urb_grid)] <- 0
          
        }
        c_pop_rur_grid <-  overlay(c_pop_grid, c_urb_mask,fun=function(x,y){return(x*!y)}) #c_pop_grid * !c_urb_mask
        
        c_pop_rur_grid <- setValues(c_pop_rur_grid,c_pop_rur_grid[,]* c_pop_rur / sum(c_pop_rur_grid[,],na.rm=T ))
        #gridded rural pop
        
        c_pop_grid <- overlay(c_pop_urb_grid,c_pop_rur_grid,fun=function(x,y){return(x+y)}) #c_pop_urb_grid + c_pop_rur_grid #country combined pop
        
        c_pop_grid
        
      }else{


	  	country <-  world[world.df$id==countryid,]
      	a <- intersect(extent(pop_grid),extent(country)) #check if country should be processed. only if within area chosen
      	c <- FALSE
      	if(!is.null(a)) c <- round(a) == round(extent(country)) #if the interesct make sure that the entire country is present in data
      	if(c){
	 		c_pop_grid <- mask(crop(pop_grid,country),country)
			c_pop_grid
		} else{
			NA
		} 

		
      } # end if nr
      
      
      
    }  #end countryid.
    
    
    
    raster.list <- list(); pp <- 1
    for(i in 1:length(c_pop_grid.list)){ 
      val <- c_pop_grid.list[[i]]
      if(length(val)>1){
        raster.list[[pp]] <- val  
        pp <- pp+1
      }
      
    }
    
    
    
    
    
    raster.list$fun <- mean
    pop_result_grid <- do.call(merge,raster.list)
    
    
    writeRaster(pop_result_grid,outfile_main,format='GTiff',NAflag=-9,overwrite=T)
    
   
  }
  
}


  

