

#################### MERRA2  data  #########################
###########################################################

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


library(ncdf4)
require(kohonen)
filename <- list.files()
year <- 1980:2016
month <- c("01","02","03","04","05","06","07","08","09","10","11","12")
days_noleap = c(31,28,31,30,31,30,31,31,30,31,30,31)
days_cum <- cumsum(days_noleap)
days_leap   = c(31,29,31,30,31,30,31,31,30,31,30,31)
days_cum_leap <- cumsum(days_leap)

####  Longitude 73:113  ####### Latitude  9:25
for(i_region in 1:length(region)){
  
  hgt.summer <- hgt.winter <-  hgt.spring <- hgt.fall <- data.reshape <- NULL
  for(i in 1:length(filename)){
    file <- nc_open(filename = filename[i])
    #   level <- ncvar_get(file, varid="Height")
    #   lat <- ncvar_get(file, varid="YDim")
    #   lon <- ncvar_get(file, varid="XDim")
    #   time <- ncvar_get(file, varid="TIME")
    if(length(region[[i_region]])>4){
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H")
      data <- data[c(corner[1]:corner[2],corner[3]:corner[4]),corner[5]:corner[6]]
    }
    else{
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H",start=c(corner[1],corner[3],1,1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1,-1))
    }
    data.reshape <- rbind(data.reshape, t(as.vector(data))) 
    rm(data)
    nc_close(file)
    
  }
  
  
  m <- 0
  for( i in 1:length(year)){
    if(year[i]%%4 !=0){
      data.winter <- data.reshape[m+c(1:days_cum[2],(days_cum[11]+1):days_cum[12]),]
      data.summer <- data.reshape[m+c((days_cum[2]+1):days_cum[5]),] 
      data.spring <- data.reshape[m+c((days_cum[5]+1):days_cum[8]),] 
      data.fall <- data.reshape[m+c((days_cum[8]+1):days_cum[11]),] 
      ## hgt[lon,lat,level,time] 
      hgt.summer <- rbind(hgt.summer,data.summer)  
      hgt.winter <- rbind(hgt.winter,data.winter)
      hgt.spring <- rbind(hgt.spring,data.spring) 
      hgt.fall <- rbind(hgt.fall,data.fall) 
      m <- m + 365
    }
    else{
      data.winter <- data.reshape[m+c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]),]
      data.summer <- data.reshape[m+c((days_cum_leap[5]+1):days_cum_leap[8]),]
      data.spring <- data.reshape[m+c((days_cum_leap[2]+1):days_cum_leap[5]),]
      data.fall <- data.reshape[m+c((days_cum_leap[8]+1):days_cum_leap[11]),]
      ## hgt[lon,lat,level,time] 
      hgt.summer <- rbind(hgt.summer,data.summer)
      hgt.winter <- rbind(hgt.winter,data.winter)
      hgt.spring <- rbind(hgt.spring,data.spring)
      hgt.fall <- rbind(hgt.fall,data.fall)
      m <- m + 366
    }
  }
  rm(data.reshape)
  
  hgt.spring1 <- hgt.summer1 <- matrix(0,nrow=nrow(hgt.summer),ncol = ncol(hgt.summer))
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
       file=paste("MERRA.hgt.",i_region,".RData",sep=""))
  rm(hgt.summer,hgt.winter,hgt.spring, hgt.fall)
  
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
       file=paste("MERRA.hgt_SOM_",i_region,".RData",sep=""))
  
}


region <-  list(WNA=c(1,121,261,354),
                CNA=c(91,191,261,354),
                ENA=c(129,241,261,354))  ############### region for Canada study ##############
filename <- list.files()

####  Longitude 73:113  ####### Latitude  9:25
for(i_region in 1){ ##:length(region)
  
  hgt.summer <- hgt.winter <-  hgt.spring <- hgt.fall <- data.reshape <- NULL
  
  for(i in 1:length(filename)){
    file <- nc_open(filename = filename[i])
      # level <- ncvar_get(file, varid="lev")
      # lat <- ncvar_get(file, varid="lat")
      # lon <- ncvar_get(file, varid="lon")
      # save(lon,lat,file="Lon_lat.RData")
      # time <- ncvar_get(file, varid="time")
    if(length(region[[i_region]])>4){
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H")
      data <- data[c(corner[1]:corner[2],corner[3]:corner[4]),corner[5]:corner[6]]
    }
    else{
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H",start=c(corner[1],corner[3],1,1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1,-1))
    }
    data.reshape <- rbind(data.reshape, t(as.vector(data))) 
    rm(data)
    nc_close(file)
    
  }
  
  load(file="Canada_all.RData")
  
  m <- 0
  for( i in 1:length(year)){
    if(year[i]%%4 !=0){
      data.winter <- data.reshape[m+c(1:days_cum[2],(days_cum[11]+1):days_cum[12]),]
      data.summer <- data.reshape[m+c((days_cum[2]+1):days_cum[5]),] 
      data.spring <- data.reshape[m+c((days_cum[5]+1):days_cum[8]),] 
      data.fall <- data.reshape[m+c((days_cum[8]+1):days_cum[11]),] 
      ## hgt[lon,lat,level,time] 
      hgt.summer <- rbind(hgt.summer,data.summer)  
      hgt.winter <- rbind(hgt.winter,data.winter)
      hgt.spring <- rbind(hgt.spring,data.spring) 
      hgt.fall <- rbind(hgt.fall,data.fall) 
      m <- m + 365
    }
    else{
      data.winter <- data.reshape[m+c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]),]
      data.summer <- data.reshape[m+c((days_cum_leap[5]+1):days_cum_leap[8]),]
      data.spring <- data.reshape[m+c((days_cum_leap[2]+1):days_cum_leap[5]),]
      data.fall <- data.reshape[m+c((days_cum_leap[8]+1):days_cum_leap[11]),]
      ## hgt[lon,lat,level,time] 
      hgt.summer <- rbind(hgt.summer,data.summer)
      hgt.winter <- rbind(hgt.winter,data.winter)
      hgt.spring <- rbind(hgt.spring,data.spring)
      hgt.fall <- rbind(hgt.fall,data.fall)
      m <- m + 366
    }
  }
  rm(data.reshape)
  
  hgt.spring1 <- hgt.summer1 <- matrix(0,nrow=nrow(hgt.summer),ncol = ncol(hgt.summer))
  
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
       file=paste("MERRA2.hgt.canada.",i_region,".RData",sep=""))
  
  # save(hgt.summer1,hgt.winter1,hgt.spring1, hgt.fall1, 
  #      file=paste("MERRA2.hgt.canada.","_All",".RData",sep=""))
  
  rm(hgt.summer,hgt.winter,hgt.spring, hgt.fall)
  
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
       file=paste("MERRA2.hgt_SOM_canada_",i_region,".RData",sep=""))
  
}


load(file=paste("MERRA2.hgt.canada.",i_region,".RData",sep=""))


for(i_region in 1){ ##:length(region)
  
  data.reshape <- NULL
  
  for(i in 1:5000){
    file <- nc_open(filename = filename[i])
    # level <- ncvar_get(file, varid="lev")
    # lat <- ncvar_get(file, varid="lat")
    # lon <- ncvar_get(file, varid="lon")
    # save(lon,lat,file="Lon_lat.RData")
    # time <- ncvar_get(file, varid="time")
    if(length(region[[i_region]])>4){
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H")
      data <- data[c(corner[1]:corner[2],corner[3]:corner[4]),corner[5]:corner[6]]
    }
    else{
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H",start=c(corner[1],corner[3],1,1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1,-1))
    }
    data.reshape <- rbind(data.reshape, t(as.vector(data))) 
    rm(data)
    nc_close(file)
    
  }
  
  save(data.reshape,file="Canada_1.RData")
}

for(i_region in 1){ ##:length(region)
  
  data.reshape <- NULL
  
  for(i in 5001:9999){
    file <- nc_open(filename = filename[i])
    # level <- ncvar_get(file, varid="lev")
    # lat <- ncvar_get(file, varid="lat")
    # lon <- ncvar_get(file, varid="lon")
    # save(lon,lat,file="Lon_lat.RData")
    # time <- ncvar_get(file, varid="time")
    if(length(region[[i_region]])>4){
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H")
      data <- data[c(corner[1]:corner[2],corner[3]:corner[4]),corner[5]:corner[6]]
    }
    else{
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H",start=c(corner[1],corner[3],1,1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1,-1))
    }
    data.reshape <- rbind(data.reshape, t(as.vector(data))) 
    rm(data)
    nc_close(file)
    
  }
  
  save(data.reshape,file="Canada_3.RData")
}



for(i_region in 1){ ##:length(region)
  
  data.reshape <- NULL
  
  for(i in 10001:length(filename)){
    file <- nc_open(filename = filename[i])
    # level <- ncvar_get(file, varid="lev")
    # lat <- ncvar_get(file, varid="lat")
    # lon <- ncvar_get(file, varid="lon")
    # save(lon,lat,file="Lon_lat.RData")
    # time <- ncvar_get(file, varid="time")
    if(length(region[[i_region]])>4){
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H")
      data <- data[c(corner[1]:corner[2],corner[3]:corner[4]),corner[5]:corner[6]]
    }
    else{
      corner <- region[[i_region]] 
      data <- ncvar_get(file,varid="H",start=c(corner[1],corner[3],1,1),
                        count=c(corner[2]-corner[1]+1,corner[4]-corner[3]+1,-1,-1))
    }
    data.reshape <- rbind(data.reshape, t(as.vector(data))) 
    rm(data)
    nc_close(file)
    
  }
  
  save(data.reshape,file="Canada_2.RData")
}


load(file=paste("MERRA2.hgt.canada.","_All",".RData",sep=""))

############## annual mean GPH anomalies ####################

hgt.spring <- hgt.summer <- hgt.fall <- hgt.winter <- matrix(NA,length(year),ncol(hgt.spring1))
m <- 0
for(i in 1:length(year)){
  hgt.spring[i,] <- colMeans(hgt.spring1[((i-1)*92+1):(i*92),])
  hgt.summer[i,] <- colMeans(hgt.summer1[((i-1)*92+1):(i*92),])
  hgt.fall[i,] <- colMeans(hgt.fall1[((i-1)*91+1):(i*91),])
  
  if(year[i] %%4 == 0){
    hgt.winter[i,] <- colMeans(hgt.winter1[m:(m+91),])
    m <- m + 91
  }
  else{
    hgt.winter[i,] <- colMeans(hgt.winter1[m:(m+90),])
    m <- m + 90
  }
}


lat <- lat[261:354]
lon <- lon[1:241]
grid <- expand.grid(y=lon,x=lat)
load(file="Lon_lat_Canada_all.RData")

library(raster)

############ spring ##########
sig <- trend <- matrix(NA, ncol(hgt.spring), 1 )
for(i in 1:ncol(hgt.spring)){
  trend[i] <- lm(hgt.spring[,i]~year)$coefficients[2]
  sig[i] <- mk.test(hgt.spring[,i])$p.value
}

dat <- data.frame(x=grid$y,y=grid$x,z=sig)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("sig_spring_",i,".nc") )

dat <- data.frame(x=grid$y,y=grid$x,z=trend)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("trend_spring_",i,".nc") )

########### summer ###############

sig <- trend <- matrix(NA, ncol(hgt.summer), 1 )
for(i in 1:ncol(hgt.summer)){
  trend[i] <- lm(hgt.summer[,i]~year)$coefficients[2]
  sig[i] <- mk.test(hgt.summer[,i])$p.value
}

dat <- data.frame(x=grid$y,y=grid$x,z=sig)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("sig_summer_",i,".nc") )

dat <- data.frame(x=grid$y,y=grid$x,z=trend)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("trend_summer_",i,".nc") )

#################### fall ###########################
sig <- trend <- matrix(NA, ncol(hgt.fall), 1 )
for(i in 1:ncol(hgt.fall)){
  trend[i] <- lm(hgt.fall[,i]~year)$coefficients[2]
  sig[i] <- mk.test(hgt.fall[,i])$p.value
}

dat <- data.frame(x=grid$y,y=grid$x,z=sig)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("sig_fall_",i,".nc") )

dat <- data.frame(x=grid$y,y=grid$x,z=trend)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("trend_fall_",i,".nc") )

################## winter ##########################
sig <- trend <- matrix(NA, ncol(hgt.winter), 1 )
for(i in 1:ncol(hgt.winter)){
  trend[i] <- lm(hgt.winter[,i]~year)$coefficients[2]
  sig[i] <- mk.test(hgt.winter[,i])$p.value
}

dat <- data.frame(x=grid$y,y=grid$x,z=sig)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("sig_winter_",i,".nc") )

dat <- data.frame(x=grid$y,y=grid$x,z=trend)
phgt <- rasterFromXYZ(xyz=dat)
crs(phgt) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# plot(phgt)
writeRaster(phgt,filename = paste0("trend_winter_",i,".nc") )


