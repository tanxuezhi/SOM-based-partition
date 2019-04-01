
# module load R/3.3.1
# export R_LIBS_USER=/home/xtan/R/library

rm(list=ls())

library(ncdf4)
require(kohonen)
year <- 1979:2016
month <- c("01","02","03","04","05","06","07","08","09","10","11","12")
days_noleap = c(31,28,31,30,31,30,31,31,30,31,30,31)
days_cum <- cumsum(days_noleap)
days_leap   = c(31,29,31,30,31,30,31,31,30,31,30,31)
days_cum_leap <- cumsum(days_leap)
filename <- list.files()

region <-  list(WNA=c(73,105,3,21),
                CNA=c(97,121,3,21),
                ENA=c(107,137,3,21))  ############### region for Canada study ##############

region <-  list(WNA=c(85,99,13,25),
                CNA=c(99,109,13,25),
                ENA=c(109,127,13,25),
                WEU=c(139,144,1,11,13,25),
                CEA=c(11,27,13,25),
                CA=c(27,43,13,25),
                EA=c(43,59,13,25),
                SA=c(25,49,37,49),
                NSA=c(117,133,29,41),
                SSA=c(117,133,13,29),
                NAF=c(137,144,1,21,37,49),
                SAF=c(137,144,1,21,25,37),
                OCE=c(21,44,21,37))

for(i_region in 1:length(region)){
  
  hgt.summer <- hgt.winter <- hgt.spring <- hgt.fall <- NULL
  
  for(i in 2143:length(filename)){
    file <- nc_open(filename = filename[i])
    lat <- ncvar_get(file, varid="lat")
    lon <- ncvar_get(file, varid="lon")
    time <- ncvar_get(file, varid="time")
    if (length(region[[i_region]])>4){
      corner <- region[[i_region]]
      data <- ncvar_get(file, varid="HGT_L100",start=c(corner[1],corner[5],1),
                        count=c(corner[2]-corner[1]+1,corner[6]-corner[5]+1,-1))
      data1 <- ncvar_get(file, varid="HGT_L100",start=c(corner[3],corner[5],1),
                         count=c(corner[4]-corner[3]+1,corner[6]-corner[5]+1,-1))
      library(abind)
      data <- abind(data,data1,along = 1)
    }
    else{
      corner <- region[[i_region]]
      data <- ncvar_get(file, varid="HGT_L100",start=c(corner[1],corner[3],1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1))
    }
    nc_close(file)
    data1 <- array(0,dim=c(dim(data)[1],dim(data)[2],dim(data)[3]/4))
    for(iii in 1:(dim(data)[3]/4)){
      data1[,,iii] <- apply(X=data[,,((iii-1)*4+1):(iii*4)],MARGIN = c(1,2),FUN=mean,na.rm=T)
    }
    data.reshape <- matrix(NA,dim(data1)[3],dim(data1)[1]*dim(data1)[2])
    for (j in 1:dim(data1)[3]){
      data.reshape[j,] <- t(as.vector(data1[,,j]))
    }
    
    if((i-2142)%%72 %in% c(1:12,67:71,0)){
      hgt.winter <- rbind(hgt.winter,data.reshape)
    }
    else if((i-2142)%%72 %in% c(31:48)){
      hgt.summer <- rbind(hgt.summer,data.reshape)
    }
    else if((i-2142)%%72 %in% c(13:30)){
      hgt.spring <- rbind(hgt.spring,data.reshape) 
    }
    else{
      hgt.fall <- rbind(hgt.fall,data.reshape) 
    }
    
  }
    for(i in 1:2102){
      file <- nc_open(filename = filename[i])
      lat <- ncvar_get(file, varid="lat")
      lon <- ncvar_get(file, varid="lon")
      time <- ncvar_get(file, varid="time")
      if (length(region[[i_region]])>4){
        corner <- region[[i_region]]
        data <- ncvar_get(file, varid="HGT_L100",start=c(corner[1],corner[5],1),
                          count=c(corner[2]-corner[1]+1,corner[6]-corner[5]+1,-1))
        data1 <- ncvar_get(file, varid="HGT_L100",start=c(corner[3],corner[5],1),
                           count=c(corner[4]-corner[3]+1,corner[6]-corner[5]+1,-1))
        library(abind)
        data <- abind(data,data1,along = 1)
      }
      else{
        corner <- region[[i_region]]
        data <- ncvar_get(file, varid="HGT_L100",start=c(corner[1],corner[3],1),
                          count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1))
      }
      nc_close(file)
      data1 <- array(0,dim=c(dim(data)[1],dim(data)[2],dim(data)[3]/4))
      for(iii in 1:(dim(data)[3]/4)){
        data1[,,iii] <- apply(X=data[,,((iii-1)*4+1):(iii*4)],MARGIN = c(1,2),FUN=mean,na.rm=T)
      }
      data.reshape <- matrix(NA,dim(data1)[3],dim(data1)[1]*dim(data1)[2])
      for (j in 1:dim(data1)[3]){
        data.reshape[j,] <- t(as.vector(data1[,,j]))
      }
      
      if(substr(filename[i],start = 11,stop = 12) %in% c("12","01","02")){
        hgt.winter <- rbind(hgt.winter,data.reshape)
      }
      else if(substr(filename[i],start = 11,stop = 12) %in% c("06","07","08")){
        hgt.summer <- rbind(hgt.summer,data.reshape)
      }
      else if(substr(filename[i],start = 11,stop = 12) %in% c("03","04","05")){
        hgt.spring <- rbind(hgt.spring,data.reshape) 
      }
      else{
        hgt.fall <- rbind(hgt.fall,data.reshape) 
      }
      
  }
  hgt.spring1 <- hgt.summer1 <- matrix(0,nrow=nrow(hgt.summer),ncol = ncol(hgt.summer))
  for(days in 1:92){
    for(jjj in seq(days,nrow(hgt.summer),length(year))){
      # hgt.summer1[jjj,] <- hgt.summer[jjj,]- colMeans(hgt.summer[seq(days,nrow(hgt.summer),length(year)),],na.rm = T)
      # hgt.spring1[jjj,] <- hgt.spring[jjj,]- colMeans(hgt.spring[seq(days,nrow(hgt.summer),length(year)),],na.rm = T)
      hgt.summer1[jjj,] <- hgt.summer[jjj,]- mean(hgt.summer[jjj,],na.rm = T)
      hgt.spring1[jjj,] <- hgt.spring[jjj,]- mean(hgt.spring[jjj,],na.rm = T)
    }
  }
  hgt.fall1 <- matrix(0,nrow=nrow(hgt.fall),ncol = ncol(hgt.fall))
  for(days in 1:91){
    for(jjj in seq(days,nrow(hgt.fall),length(year))){
      # hgt.fall1[jjj,] <- hgt.fall[jjj,]- colMeans(hgt.fall[seq(days,nrow(hgt.fall),length(year)),],na.rm = T)
      hgt.fall1[jjj,] <- hgt.fall[jjj,]- mean(hgt.fall[jjj,],na.rm = T)
      
    }
  }
  
  days.winter <- NULL
  for(ii in 1:length(year) ){
    if(year[ii]%%4 !=0){
      days.winter <- c(days.winter,c(1:days_cum[2],(days_cum[11]+1):days_cum[12]))
    }
    else{
      days.winter <- c(days.winter,c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]))
    }
  }
  hgt.winter1 <- matrix(0,nrow=nrow(hgt.winter),ncol = ncol(hgt.winter))
  for(days in 1:nrow(hgt.winter)){
    # hgt.winter1[days,] <- hgt.winter[days,] - colMeans(hgt.winter[days.winter==days.winter[days],],na.rm = T)
    hgt.winter1[days,] <- hgt.winter[days,] - mean(hgt.winter[days,],na.rm = T)
  }
  
  save(hgt.summer1,hgt.winter1,hgt.spring1, hgt.fall1, 
       file=paste("CFSR.hgt.",i_region,".RData",sep=""))
  
  require(kohonen)
  som_grid <- somgrid(xdim = 4, ydim=4, topo="rectangular")  
  
  # Train the SOM model!
  system.time(som_model.hgt.winter <- som(hgt.winter1, 
                                          grid=som_grid, 
                                          rlen=100, 
                                          alpha=c(0.05,0.01), 
                                          n.hood = "circular",
                                          keep.data = TRUE ))
  
  system.time(som_model.hgt.summer <- som(hgt.summer1, 
                                          grid=som_grid, 
                                          rlen=100, 
                                          alpha=c(0.05,0.01), 
                                          n.hood = "circular",
                                          keep.data = TRUE ))
  
  system.time(som_model.hgt.spring <- som(hgt.spring1, 
                                          grid=som_grid, 
                                          rlen=100, 
                                          alpha=c(0.05,0.01), 
                                          n.hood = "circular",
                                          keep.data = TRUE ))
  
  system.time(som_model.hgt.fall <- som(hgt.fall1, 
                                        grid=som_grid, 
                                        rlen=100, 
                                        alpha=c(0.05,0.01), 
                                        n.hood = "circular",
                                        keep.data = TRUE ))
  
  som_grid <- somgrid(xdim = 3, ydim=3, topo="rectangular")  
  
  # Train the SOM model!
  system.time(som_model.hgt.winter_3_3 <- som(hgt.winter1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.summer_3_3 <- som(hgt.summer1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.spring_3_3 <- som(hgt.spring1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.fall_3_3 <- som(hgt.fall1, 
                                            grid=som_grid, 
                                            rlen=100, 
                                            alpha=c(0.05,0.01), 
                                            n.hood = "circular",
                                            keep.data = TRUE ))
  
  som_grid <- somgrid(xdim = 2, ydim=2, topo="rectangular")  
  
  # Train the SOM model!
  system.time(som_model.hgt.winter_2_2 <- som(hgt.winter1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.summer_2_2 <- som(hgt.summer1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.spring_2_2 <- som(hgt.spring1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.fall_2_2 <- som(hgt.fall1, 
                                            grid=som_grid, 
                                            rlen=100, 
                                            alpha=c(0.05,0.01), 
                                            n.hood = "circular",
                                            keep.data = TRUE ))
  
  #   class <- som_model.hgt.summer$unit.classif
  #   class <- som_model.hgt.winter$unit.classif
  #   class <- som_model.hgt.spring$unit.classif
  #   class <- som_model.hgt.fall$unit.classif
  save(som_model.hgt.summer,som_model.hgt.winter,som_model.hgt.spring,som_model.hgt.fall,
       som_model.hgt.summer_3_3,som_model.hgt.winter_3_3,som_model.hgt.spring_3_3,som_model.hgt.fall_3_3,
       som_model.hgt.summer_2_2,som_model.hgt.winter_2_2,som_model.hgt.spring_2_2,som_model.hgt.fall_2_2,
       file=paste("CFSR.hgt_SOM_",i_region,".RData",sep=""))
}

region <-  list(WNA=c(85,99,13,25),
                CNA=c(99,109,13,25),
                ENA=c(109,127,13,25),
                WEU=c(139,144,1,11,13,25),
                CEA=c(11,27,13,25),
                CA=c(27,43,13,25),
                EA=c(43,59,13,25),
                SA=c(25,49,37,49),
                NSA=c(117,133,29,41),
                SSA=c(117,133,13,29),
                NAF=c(137,144,1,21,37,49),
                SAF=c(137,144,1,21,25,37),
                OCE=c(21,44,21,37))


load(file="Lon_lat.RData")



for(i_region in 1:4){
  ######## region west ###############
  if(length(region[[i_region]])<5){
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[reg[1]:reg[2]]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[3]:reg[4]]
    grid <- expand.grid(y=lon,x=lat)
  }
  else{
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[c(reg[1]:reg[2],reg[3]:reg[4])]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[5]:reg[6]]
    grid <- expand.grid(y=lon,x=lat)
  }

  
  load(file=paste("CFSR.hgt_SOM_canada_",i_region,".RData",sep=""))
  pattern <- som_model.hgt.spring_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("test.som_2_2_spring_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.summer_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_summer_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.fall_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_fall_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.winter_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_winter_region_canada_",i_region),mcol = 2)
  
  pattern <- som_model.hgt.spring_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_spring_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.summer_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_summer_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.fall_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_fall_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.winter_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_winter_region_canada_",i_region),mcol = 3)
  
  pattern <- som_model.hgt.spring$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_spring_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.summer$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_summer_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.fall$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_fall_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.winter$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_winter_region_canada_",i_region),mcol = 4)
}

pattern.plot <- function(pattern,lon.west,lat.west,filename,mcol){ 
  
  data <- NULL
  for(j in 1:nrow(pattern)){
    data <- rbind(data,cbind(lon.west,lat.west,pattern[j,],j))
  }
  data <- as.data.frame(data)
  colnames(data) <- c("lon","lat","hgt","panel")
  
  library('data.table')
  library('maps')
  library('ggplot2')
  library('dplyr')
  library("PBSmapping")
  
  # plotting region
  xlim = range(lon.west)
  ylim = range(lat.west)
  
  npanels = nrow(pattern)
  
  # function to get map data for plotting
  get_clipped_map_data = function(mapname,xlim,ylim){
    require('data.table')
    require('PBSmapping')
    
    map = map_data(mapname)
    setnames(map, c("X","Y","PID","POS","region","subregion"))
    clipped_map = clipPolys(map, xlim=xlim,ylim=ylim, keepExtra=TRUE)
    setnames(clipped_map, c("X","Y","PID"), c("lon","lat","group"))
    return(as.data.table(as.data.frame(clipped_map)))
  }
  
  # statemap = get_clipped_map_data('state',xlim,ylim)
  worldmap = get_clipped_map_data('world',xlim,ylim)
  
  # generate random data 
  # grid_data = as.data.table(expand.grid(
  #   lon=seq(xlim[1],xlim[2],length.out=20),
  #   lat=seq(ylim[1],ylim[2],length.out=20),
  #   panel=1:4))
  # 
  # grid_data = mutate(group_by(grid_data,panel),value=rnorm(length(lat),sd=panel))
  
  p = ggplot()+
    # gridded data
    geom_raster(data=data,aes(x=lon,y=lat,fill=hgt), alpha = 0.8)+
    # contours 
    stat_contour(data = data, aes(x = lon, y = lat, z = hgt),col="black") +
    # states
    # geom_polygon(data=statemap,aes(lon,lat,group=group),fill=NA, color="grey50")+
    # panels
    facet_wrap(~panel,ncol=mcol)+
    # other countries
    geom_polygon(data=worldmap,aes(lon,lat,group=group),fill=NA, color="grey50")+
    # project
    # coord_map('lambert',xlim=xlim,ylim=ylim,lat0=ylim[1],lat1=ylim[2])+
    coord_quickmap() +
    # white background
    theme_bw()+
    # color range
    scale_fill_gradientn(colours=c("blue", "white","red" ), 
                           name="Geopotential\nHeight (m)")+
    # turn off plot elements
    theme(panel.border=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          line=element_blank(),legend.key.height=unit(7.5,"line"), 
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          strip.text = element_text(size=16))
  
  # print(p)
  ggsave(filename=paste0(filename,".pdf"),p,width=12,height=10)
}



######################### ERA data ##########################
setwd("C:/ExtremeAttribution/data/ERA_interim/hgt")
region <-  list(WNA=c(241,348,7,68),
                CNA=c(321,401,7,68),
                ENA=c(354,455,7,68))  ############### region for Canada study ##############

region <-  list(WNA=c(281,328,41,81),
                CNA=c(328,361,41,81),
                ENA=c(361,421,41,81),
                WEU=c(461,480,1,35,41,81),
                CEA=c(35,88,41,81),
                CA=c(88,141,41,81),
                EA=c(81,161,81,121),
                SA=c(86,172,85,129),
                NSA=c(338,441,107,148),
                SSA=c(338,441,148,201),
                NAF=c(455,480,1,68,81,121),
                SAF=c(455,480,1,68,121,161),
                OCE=c(135,215,121,175))


for(i_region in 1:length(region)){
  ######## region west ###############
  if(length(region[[i_region]])<5){
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[reg[1]:reg[2]]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[3]:reg[4]]
    grid <- expand.grid(y=lon,x=lat)
  }
  else{
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[c(reg[1]:reg[2],reg[3]:reg[4])]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[5]:reg[6]]
    grid <- expand.grid(y=lon,x=lat)
  }
  
  
  load(file=paste("ERA_Interim.hgt_SOM_canada_",i_region,"_1.RData",sep=""))
  pattern <- som_model.hgt.spring_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_spring_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.summer_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_summer_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.fall_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_fall_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.winter_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_winter_region_canada_",i_region),mcol = 2)
  
  pattern <- som_model.hgt.spring_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_spring_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.summer_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_summer_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.fall_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_fall_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.winter_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_winter_region_canada_",i_region),mcol = 3)
  
  pattern <- som_model.hgt.spring$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_spring_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.summer$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_summer_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.fall$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_fall_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.winter$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_winter_region_canada_",i_region),mcol = 4)
}


###################  JRA-55 ##################################
setwd("C:\\ExtremeAttribution\\data\\JARR\\hgt\\hgt.anomalies")
region <-  list(WNA=c(145,209,5,41),
                CNA=c(193,241,5,41),
                ENA=c(213,273,5,41))  ############### region for Canada study ##############

region <-  list(WNA=c(160,197,25,49),
                CNA=c(197,217,25,49),
                ENA=c(217,253,25,49),
                WEU=c(277,288,1,21,25,49),
                CEA=c(21,53,25,49),
                CA=c(53,85,25,49),
                EA=c(85,117,25,49),
                SA=c(49,97,49,73),
                NSA=c(233,265,65,89),
                SSA=c(233,265,89,121),
                NAF=c(277,288,1,41,49,73),
                SAF=c(273,288,1,41,73,97),
                OCE=c(81,129,73,115))
for(i_region in 1:length(region)){
  ######## region west ###############
  if(length(region[[i_region]])<5){
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[reg[1]:reg[2]]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[3]:reg[4]]
    grid <- expand.grid(y=lon,x=lat)
  }
  else{
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[c(reg[1]:reg[2],reg[3]:reg[4])]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[5]:reg[6]]
    grid <- expand.grid(y=lon,x=lat)
  }
  
  
  load(file=paste("JRA55.hgt_SOM_canada_",i_region,".RData",sep=""))
  pattern <- som_model.hgt.spring_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_spring_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.summer_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_summer_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.fall_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_fall_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.winter_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_winter_region_canada_",i_region),mcol = 2)
  
  pattern <- som_model.hgt.spring_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_spring_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.summer_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_summer_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.fall_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_fall_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.winter_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_winter_region_canada_",i_region),mcol = 3)
  
  pattern <- som_model.hgt.spring$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_spring_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.summer$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_summer_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.fall$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_fall_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.winter$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_winter_region_canada_",i_region),mcol = 4)
}

#################### MERRA #########################################
setwd("C:\\ExtremeAttribution\\data\\MERRA\\hgt\\hgt_anomalies")
region <-  list(WNA=c(1,65,104,141),
                CNA=c(48,97,104,141),
                ENA=c(68,129,104,141))  ############### region for Canada study ##############

region <-  list(WNA=c(24,52,96,121),
                CNA=c(52,73,96,121),
                ENA=c(73,109,96,121),
                WEU=c(132,165,96,121),
                CEA=c(165,229,96,121),
                CA=c(196,229,96,121),
                EA=c(229,261,96,121),
                SA=c(192,241,72,97),
                NSA=c(88,121,56,81),
                SSA=c(88,121,24,56),
                NAF=c(128,185,72,97),
                SAF=c(128,185,48,72),
                OCE=c(224,273,40,72))

for(i_region in 1:length(region)){
  ######## region west ###############
  if(length(region[[i_region]])<5){
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[reg[1]:reg[2]]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[3]:reg[4]]
    grid <- expand.grid(y=lon,x=lat)
  }
  else{
    load(file="Lon_lat.RData")
    reg <- region[[i_region]]
    lon <- lon[c(reg[1]:reg[2],reg[3]:reg[4])]
    lon[lon>=180] <- -(360-lon[lon>=180])
    lat <- lat[reg[5]:reg[6]]
    grid <- expand.grid(y=lon,x=lat)
  }
  
  
  load(file=paste("MERRA.hgt_SOM_canada_",i_region,".RData",sep=""))
  pattern <- som_model.hgt.spring_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_spring_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.summer_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_summer_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.fall_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_fall_region_canada_",i_region),mcol = 2)
  pattern <- som_model.hgt.winter_2_2$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_2_2_winter_region_canada_",i_region),mcol = 2)
  
  pattern <- som_model.hgt.spring_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_spring_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.summer_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_summer_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.fall_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_fall_region_canada_",i_region),mcol = 3)
  pattern <- som_model.hgt.winter_3_3$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_3_3_winter_region_canada_",i_region),mcol = 3)
  
  pattern <- som_model.hgt.spring$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_spring_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.summer$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_summer_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.fall$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_fall_region_canada_",i_region),mcol = 4)
  pattern <- som_model.hgt.winter$codes
  pattern.plot(pattern=pattern,lon.west=grid$y,lat.west=grid$x,filename=paste0("som_4_4_winter_region_canada_",i_region),mcol = 4)
}

#######################################################################
for(i_region in 1:length(region)){
  
  hgt.summer <- hgt.winter <- hgt.spring <- hgt.fall <- NULL
  
  for(i in 1:length(filename)){
    file <- nc_open(filename = filename[i])
    lat <- ncvar_get(file, varid="lat")
    lon <- ncvar_get(file, varid="lon")
    time <- ncvar_get(file, varid="time")
    if (length(region[[i_region]])>4){
      corner <- region[[i_region]]
      data <- ncvar_get(file, varid="HGT_L100",start=c(corner[1],corner[5],1),
                        count=c(corner[2]-corner[1]+1,corner[6]-corner[5]+1,-1))
      data1 <- ncvar_get(file, varid="HGT_L100",start=c(corner[3],corner[5],1),
                         count=c(corner[4]-corner[3]+1,corner[6]-corner[5]+1,-1))
      library(abind)
      data <- abind(data,data1,along = 1)
    }
    else{
      corner <- region[[i_region]]
      data <- ncvar_get(file, varid="HGT_L100",start=c(corner[1],corner[3],1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1))
    }
    nc_close(file)
    data1 <- array(0,dim=c(dim(data)[1],dim(data)[2],dim(data)[3]/4))
    for(iii in 1:(dim(data)[3]/4)){
      data1[,,iii] <- apply(X=data[,,((iii-1)*4+1):(iii*4)],MARGIN = c(1,2),FUN=mean,na.rm=T)
    }
    data.reshape <- matrix(NA,dim(data1)[3],dim(data1)[1]*dim(data1)[2])
    for (j in 1:dim(data1)[3]){
      data.reshape[j,] <- t(as.vector(data1[,,j]))
    }
    if(i%%72 %in% c(1:12,67:71,0)){
      hgt.winter <- rbind(hgt.winter,data.reshape)
    }
    else if(i%%72 %in% c(31:48)){
      hgt.summer <- rbind(hgt.summer,data.reshape)
    }
    else if(i%%72 %in% c(13:30)){
      hgt.spring <- rbind(hgt.spring,data.reshape) 
    }
    else{
      hgt.fall <- rbind(hgt.fall,data.reshape) 
    }
    
  }
  hgt.spring1 <- hgt.summer1 <- matrix(0,nrow=nrow(hgt.summer),ncol = ncol(hgt.summer))
  for(days in 1:92){
    for(jjj in seq(days,nrow(hgt.summer),length(year))){
      hgt.summer1[jjj,] <- hgt.summer[jjj,]- colMeans(hgt.summer[seq(days,nrow(hgt.summer),length(year)),],na.rm = T)
      hgt.spring1[jjj,] <- hgt.spring[jjj,]- colMeans(hgt.spring[seq(days,nrow(hgt.summer),length(year)),],na.rm = T)
    }
  }
  hgt.fall1 <- matrix(0,nrow=nrow(hgt.fall),ncol = ncol(hgt.fall))
  for(days in 1:91){
    for(jjj in seq(days,nrow(hgt.fall),length(year))){
      hgt.fall1[jjj,] <- hgt.fall[jjj,]- colMeans(hgt.fall[seq(days,nrow(hgt.fall),length(year)),],na.rm = T)
    }
  }
  
  days.winter <- NULL
  for(ii in 1:length(year) ){
    if(year[ii]%%4 !=0){
      days.winter <- c(days.winter,c(1:days_cum[2],(days_cum[11]+1):days_cum[12]))
    }
    else{
      days.winter <- c(days.winter,c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]))
    }
  }
  hgt.winter1 <- matrix(0,nrow=nrow(hgt.winter),ncol = ncol(hgt.winter))
  for(days in 1:nrow(hgt.winter)){
    hgt.winter1[days,] <- hgt.winter[days,] - colMeans(hgt.winter[days.winter==days.winter[days],],na.rm = T)
  }
  
  save(hgt.summer1,hgt.winter1,hgt.spring1, hgt.fall1, 
       file=paste("CFSR.hgt.canada.",i_region,".RData",sep=""))
  
  require(kohonen)
  som_grid <- somgrid(xdim = 4, ydim=4, topo="rectangular")  
  
  # Train the SOM model!
  system.time(som_model.hgt.winter <- som(hgt.winter1, 
                                          grid=som_grid, 
                                          rlen=100, 
                                          alpha=c(0.05,0.01), 
                                          n.hood = "circular",
                                          keep.data = TRUE ))
  
  system.time(som_model.hgt.summer <- som(hgt.summer1, 
                                          grid=som_grid, 
                                          rlen=100, 
                                          alpha=c(0.05,0.01), 
                                          n.hood = "circular",
                                          keep.data = TRUE ))
  
  system.time(som_model.hgt.spring <- som(hgt.spring1, 
                                          grid=som_grid, 
                                          rlen=100, 
                                          alpha=c(0.05,0.01), 
                                          n.hood = "circular",
                                          keep.data = TRUE ))
  
  system.time(som_model.hgt.fall <- som(hgt.fall1, 
                                        grid=som_grid, 
                                        rlen=100, 
                                        alpha=c(0.05,0.01), 
                                        n.hood = "circular",
                                        keep.data = TRUE ))
  
  som_grid <- somgrid(xdim = 3, ydim=3, topo="rectangular")  
  
  # Train the SOM model!
  system.time(som_model.hgt.winter_3_3 <- som(hgt.winter1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.summer_3_3 <- som(hgt.summer1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.spring_3_3 <- som(hgt.spring1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.fall_3_3 <- som(hgt.fall1, 
                                            grid=som_grid, 
                                            rlen=100, 
                                            alpha=c(0.05,0.01), 
                                            n.hood = "circular",
                                            keep.data = TRUE ))
  
  som_grid <- somgrid(xdim = 2, ydim=2, topo="rectangular")  
  
  # Train the SOM model!
  system.time(som_model.hgt.winter_2_2 <- som(hgt.winter1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.summer_2_2 <- som(hgt.summer1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.spring_2_2 <- som(hgt.spring1, 
                                              grid=som_grid, 
                                              rlen=100, 
                                              alpha=c(0.05,0.01), 
                                              n.hood = "circular",
                                              keep.data = TRUE ))
  
  system.time(som_model.hgt.fall_2_2 <- som(hgt.fall1, 
                                            grid=som_grid, 
                                            rlen=100, 
                                            alpha=c(0.05,0.01), 
                                            n.hood = "circular",
                                            keep.data = TRUE ))
  
  #   class <- som_model.hgt.summer$unit.classif
  #   class <- som_model.hgt.winter$unit.classif
  #   class <- som_model.hgt.spring$unit.classif
  #   class <- som_model.hgt.fall$unit.classif
  save(som_model.hgt.summer,som_model.hgt.winter,som_model.hgt.spring,som_model.hgt.fall,
       som_model.hgt.summer_3_3,som_model.hgt.winter_3_3,som_model.hgt.spring_3_3,som_model.hgt.fall_3_3,
       som_model.hgt.summer_2_2,som_model.hgt.winter_2_2,som_model.hgt.spring_2_2,som_model.hgt.fall_2_2,
       file=paste("CFSR.hgt_SOM_canada_",i_region,".RData",sep=""))
}
