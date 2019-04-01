

rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")

library(SDMTools);library(rgdal);library(raster);library(ggmap);library(rgeos);library(maptools)
library(dplyr);library(tidyr);library(png);library(grid);library(sp);library(kohonen)


## read precipitation and Geopotential Height data

###################### NARR Data Analysese ################################
load("C:\\ExtremeAttribution2\\NARR\\TP\\spring.precipitation.RData")
load(file="C:\\ExtremeAttribution\\data\\NARR\\hgt\\som_model.east.winter.RData")
pattern <- som_model.east.winter$unit.classif


load("C:\\ExtremeAttribution2\\NARR\\TP\\summer.precipitation.RData")

load("C:\\ExtremeAttribution2\\NARR\\TP\\fall.precipitation.RData")

load("C:\\ExtremeAttribution2\\NARR\\TP\\winter.precipitation.RData")

region <- c("west","central","east")
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL
for(j in 1:4){
  # filename1 <- paste0("C:\\ExtremeAttribution2\\NARR\\tp\\",season[j],".precipitation.RData")
  # file1 <- load(filename1)
  for(i in 1:3){
    
    pre.data <- get(paste0(region[i],".",season[j]))
    
    # filename <- paste0("C:\\ExtremeAttribution\\data\\NARR\\hgt\\som_model.",region[i],".",
    #                    season[j],".RData")
    # file <- load(filename)
    # data <- get(file)
    # results[[(j-1)*3+i]] <- winter.attribution(east.winter=get(paste0(region[i],".",season[j]))[1:length(data$unit.classif),], 
    #                                            pattern=data$unit.classif, 
    #                                            name=paste0(region[i],"_",season[j],"_"),
    #                                            season=season[j])
    # rm(som_model.west.winter,som_model.west.spring,som_model.west.summer,som_model.west.fall,
    #    som_model.central.winter,som_model.central.spring,som_model.central.summer,som_model.central.fall,
    #    som_model.east.winter,som_model.east.spring,som_model.east.summer,som_model.east.fall)
    
    for(k in 1:length(som)){
      if(k==3){
        load(paste("C:\\ExtremeAttribution\\data\\NARR\\hgt\\som_model.",region[i],".",season[j],".RData",sep=""))
        data <- get(paste0("som_model.",region[i],".",season[j]))
      }
      else{
        load(paste("C:\\ExtremeAttribution\\data\\NARR\\hgt\\som_model.",region[i],".",season[j],"_",som[k],".RData",sep=""))
        data <- get(paste0("som_model.",region[i],".",season[j]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data, 
                                                    pattern=data[1:nrow(pre.data)], 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1979,yeare=2014,n.pattern=(k+1)^2,
                                                    weight=F ))
      rm(som_model.west.spring,som_model.central.spring,som_model.east.spring,
         som_model.west.summer,som_model.central.summer,som_model.east.summer,
         som_model.west.fall,som_model.central.fall,som_model.east.fall,
         som_model.west.winter,som_model.central.winter,som_model.east.winter)
    } 
    
  }
}

# save(results.west.spring,results.central.spring,results.east.spring,
#      results.west.summer,results.central.summer,results.east.summer,
#      results.west.fall,results.central.fall,results.east.fall,
#      results.west.winter,results.central.winter,results.east.winter,
#      file="results.attribution.absolute_value.RData")

save(results,file="attribution.results.NARR.RData")


##########################################################################################################################
################################################# ERA-Interim ############################################################
##########################################################################################################################

rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL
load(file = paste0("C:\\ExtremeAttribution2\\ERA-Interim\\TP\\region.coordinations.ERA-Interim.RData"))
lat <- list(west.cor[,2],central.cor[,2],east.cor[,2])
for(j in 1:4){   ######## (j in 1:4) 
  load(file = paste0("C:\\ExtremeAttribution2\\ERA-Interim\\TP\\",season[j],".precipitation.RData"))
  for(i in 1:3){  ######### (i in 1:3)
    load(paste("C:\\ExtremeAttribution\\data\\ERA_interim\\hgt\\ERA_Interim.hgt_SOM_canada_",i,"_1.RData",sep=""))
    pre.data <- get(paste0(region[i],".",season[j]))
    for(ii in 1:nrow(pre.data)){
      for(jj in 1:ncol(pre.data)){
        if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
      }
    }
    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data[1:length(data),], 
                                                    pattern=data, 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1979,yeare=2014,
                                                    n.pattern=(k+1)^2,weight=T,wgt=cos(lat[[i]]/180*pi))) 
    } 
  }
}
save(results,file="ERA_interim.data.attribution.RData")

##########################################################################################################################
##########################################  NCEP2   data  ################################################################
##########################################################################################################################


rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL 
load(file = paste0("C:\\ExtremeAttribution2\\NCEP2\\TP\\region.coordinations.NCEP2.RData"))
lat <- list(west.cor[,2],central.cor[,2],east.cor[,2])
for(j in 1:4){   ######## (j in 1:4) 
  load(file = paste0("C:\\ExtremeAttribution2\\NCEP2\\TP\\",season[j],".precipitation.RData"))
  for(i in 1:3){  ######### (i in 1:3)
    load(paste("C:\\ExtremeAttribution\\data\\NCEP2\\hgt\\NCEP2.hgt_SOM_canada_",i,".RData",sep=""))
    pre.data <- get(paste0(region[i],".",season[j]))/10
    for(ii in 1:nrow(pre.data)){
      for(jj in 1:ncol(pre.data)){
        if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
      }
    }
    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data[1:length(data),], 
                                                    pattern=data, 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1979,yeare=2014,n.pattern=(k+1)^2,
                                                    weight=T,wgt=cos(lat[[i]]/180*pi) )) 
    } 
  }
}
save(results,file="NCEP2.data.attribution.RData")


##########################################################################################################################
##########################################  MERRA   data  ################################################################
##########################################################################################################################

rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL 
load(file = paste0("C:\\ExtremeAttribution2\\MERRA2\\TP\\region.coordinates.MERRA2.RData"))
lat <- list(west.cor[,2],central.cor[,2],east.cor[,2])
load(file = paste0("C:\\ExtremeAttribution2\\MERRA2\\TP\\Corrected_MERRA2.seasonal.regional.precipitation.data.RData"))

for(j in 1:4){   ######## (j in 1:4) 
  # load(file = paste0("C:\\ExtremeAttribution\\data\\MERRA\\tp\\",season[j],".precipitation.MERRA1.RData"))
  for(i in 1:3){  ######### (i in 1:3)
    load(paste("C:\\ExtremeAttribution2\\MERRA2\\GPH\\MERRA2.hgt_SOM_canada_",i,".RData",sep=""))
    pre.data <- get(paste0(region[i],".",season[j]))
    for(ii in 1:nrow(pre.data)){
      for(jj in 1:ncol(pre.data)){
        if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
      }
    }

    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data[1:length(data),], 
                                                    pattern=data, 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1979,yeare=2010,n.pattern=(k+1)^2,
                                                    weight=T,wgt=cos(lat[[i]]/180*pi) )) 
      rm(som_model.west.spring,som_model.central.spring,som_model.east.spring,
         som_model.west.summer,som_model.central.summer,som_model.east.summer,
         som_model.west.fall,som_model.central.fall,som_model.east.fall,
         som_model.west.winter,som_model.central.winter,som_model.east.winter)
    } 
  }
}
save(results,file="MERRA2.data.attribution.RData")

##########################################################################################################################
##########################################  JRA-55   data ################################################################
##########################################################################################################################


rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL
load(file = paste0("C:\\ExtremeAttribution2\\JARR\\TP\\region.coordinations.JARR.RData"))
lat <- list(west.cor[,2],central.cor[,2],east.cor[,2])
load("C:\\ExtremeAttribution2\\JARR\\TP\\JRA55.seasonal.regional.precipitation.data.RData")
for(j in 1:4){   ######## (j in 1:4) 
  # load(file = paste0("C:\\ExtremeAttribution\\data\\JARR\\tp\\",season[j],".precipitation.JARR.RData"))
  for(i in 1:3){  ######### (i in 1:3)
    load(paste("C:\\ExtremeAttribution\\data\\JARR\\hgt\\hgt_anomalies\\JRA55.hgt_SOM_canada_",i,".RData",sep=""))
    pre.data <- get(paste0(region[i],".",season[j]))/10
    for(ii in 1:nrow(pre.data)){
      for(jj in 1:ncol(pre.data)){
        if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
      }
    }
    pre.data <- pre.data
    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data[1:length(data),], 
                                                    pattern=data, 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1958,yeare=2013,n.pattern=(k+1)^2,
                                                    weight=T,wgt=cos(lat[[i]]/180*pi) )) 
      rm(som_model.west.spring,som_model.central.spring,som_model.east.spring,
         som_model.west.summer,som_model.central.summer,som_model.east.summer,
         som_model.west.fall,som_model.central.fall,som_model.east.fall,
         som_model.west.winter,som_model.central.winter,som_model.east.winter)
    } 
  }
}
save(results,file="JRA55.data.attribution.RData")

#########################################################################################
################################## ERA-Interim   data ###################################
################################ attribution McKenney data ##############################
#########################################################################################

rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL
# load(file = paste0("C:\\ExtremeAttribution\\data\\mck\\tp\\region.row&col.mck.RData"))
load(file = paste0("C:\\ExtremeAttribution\\data\\mck\\tp\\region.row&col.mck.RData"))
lat <- list(west.position[,2],central.position[,2],east.position[,2])
for(j in 1:4){   ######## (j in 1:4) 
    for(i in 1:3){  ######### (i in 1:3)
    # load(paste("C:\\ExtremeAttribution\\data\\ERA_interim\\hgt\\ERA_Interim.hgt_SOM_canada_",i,"_1.RData",sep=""))
    # load(file = paste0("C:\\ExtremeAttribution\\data\\mck\\",region[i],".",season[j],".mck.RData"))
    
    load(paste("ERA_Interim.hgt_SOM_canada_",i,"_1.RData",sep=""))
    load(file = paste0(region[i],".",season[j],".mck.RData"))
    
    pre.data <- get(paste0(region[i],".",season[j],"1"))
    # for(ii in 1:nrow(pre.data)){
    #   for(jj in 1:ncol(pre.data)){
    #     if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
    #   }
    # }
    rm(west.spring1,central.spring1,east.spring1,
       west.summer1,central.summer1,east.summer1,
       west.fall1,central.fall1,east.fall1,
       west.winter1,central.winter1,east.winter1)
    
    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data, 
                                                    pattern=data[1:nrow(pre.data)],  
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1979,yeare=2013,n.pattern=(k+1)^2,
                                                    weight=T,wgt=cos(lat[[i]]/180*pi) ))
      rm(som_model.west.spring,som_model.central.spring,som_model.east.spring,
         som_model.west.summer,som_model.central.summer,som_model.east.summer,
         som_model.west.fall,som_model.central.fall,som_model.east.fall,
         som_model.west.winter,som_model.central.winter,som_model.east.winter)
    } 
  }
}
save(results,file="ERA-Interim.mck.data.attribution.RData")

##########################################################################################################################
##########################################  CFSR   data ################################################################
##########################################################################################################################

load("C:\\ExtremeAttribution2\\results\\CFSR\\TP\\apcp.pre.1979-2010.RData")
load("C:\\ExtremeAttribution2\\results\\CFSR\\TP\\apcp.pre.2011-2016.RData")
apcp.west <- t(rbind(apcp.west1,apcp.west2*10))
apcp.central <- t(rbind(apcp.central1,apcp.central2*10))
apcp.east <- t(rbind(apcp.east1,apcp.east2*10))

year <- 1979:2016
month <- c("01","02","03","04","05","06","07","08","09","10","11","12")
days_noleap = c(31,28,31,30,31,30,31,31,30,31,30,31)
days_cum <- cumsum(days_noleap)
days_leap   = c(31,29,31,30,31,30,31,31,30,31,30,31)
days_cum_leap <- cumsum(days_leap)

west.winter <- central.winter <- east.winter <- west.spring <- central.spring <- east.spring <-
  west.summer <- central.summer <- east.summer <- west.fall <- central.fall <- east.fall <- NULL
m <- 0
for(i_year in 1:length(year)){  
  if(year[i_year]%%4 !=0){
    west.winter <- cbind(west.winter,apcp.west[,m+c(1:days_cum[2],(days_cum[11]+1):days_cum[12])])
    central.winter <- cbind(central.winter,apcp.central[,m+c(1:days_cum[2],(days_cum[11]+1):days_cum[12])])
    east.winter <- cbind(east.winter,apcp.east[,m+c(1:days_cum[2],(days_cum[11]+1):days_cum[12])])
    m <- m + 365
  }
  else {
    west.winter <- cbind(west.winter,apcp.west[,m+c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12])])
    central.winter <- cbind(central.winter,apcp.central[,m+c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12])])
    east.winter <- cbind(east.winter,apcp.east[,m+c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12])])
    m <- m + 366
  }
}

m <- 0
for(i_year in 1:length(year)){    
  if(year[i_year]%%4 !=0){
    west.spring <- cbind(west.spring,apcp.west[,m+c((days_cum[2]+1):days_cum[5])])
    central.spring <- cbind(central.spring,apcp.central[,m+c((days_cum[2]+1):days_cum[5])])
    east.spring <- cbind(east.spring,apcp.east[,m+c((days_cum[2]+1):days_cum[5])])
    m <- m + 365
  }
  else {
    west.spring <- cbind(west.spring,apcp.west[,m+c((days_cum_leap[2]+1):days_cum_leap[5])])
    central.spring <- cbind(central.spring,apcp.central[,m+c((days_cum_leap[2]+1):days_cum_leap[5])])
    east.spring <- cbind(east.spring,apcp.east[,m+c((days_cum_leap[2]+1):days_cum_leap[5])])
    m <- m + 366
  }
}

m <- 0
for(i_year in 1:length(year)){   
  if(year[i_year]%%4 !=0){
    west.summer <- cbind(west.summer,apcp.west[,m+c((days_cum[5]+1):days_cum[8])])
    central.summer <- cbind(central.summer,apcp.central[,m+c((days_cum[5]+1):days_cum[8])])
    east.summer <- cbind(east.summer,apcp.east[,m+c((days_cum[5]+1):days_cum[8])])
    m <- m + 365
  }
  else {
    west.summer <- cbind(west.summer,apcp.west[,m+c((days_cum_leap[5]+1):days_cum_leap[8])])
    central.summer <- cbind(central.summer,apcp.central[,m+c((days_cum_leap[5]+1):days_cum_leap[8])])
    east.summer <- cbind(east.summer,apcp.east[,m+c((days_cum_leap[5]+1):days_cum_leap[8])])
    m <- m + 366
  }
}

m <- 0
for(i_year in 1:length(year)){  
  if(year[i_year]%%4 !=0){
    west.fall <- cbind(west.fall,apcp.west[,m+c((days_cum[8]+1):days_cum[11])])
    central.fall <- cbind(central.fall,apcp.central[,m+c((days_cum[8]+1):days_cum[11])])
    east.fall <- cbind(east.fall,apcp.east[,m+c((days_cum[8]+1):days_cum[11])])
    m <- m + 365
  }
  else {
    west.fall <- cbind(west.fall,apcp.west[,m+c((days_cum_leap[8]+1):days_cum_leap[11])])
    central.fall <- cbind(central.fall,apcp.central[,m+c((days_cum_leap[8]+1):days_cum_leap[11])])
    east.fall <- cbind(east.fall,apcp.east[,m+c((days_cum_leap[8]+1):days_cum_leap[11])])
    m <- m + 366
  }
}
save(west.spring,west.summer,west.fall,west.winter,
     central.spring,central.summer,central.fall,central.winter,
     east.spring,east.summer,east.fall,east.winter,file = "seasonal.regional.precipitation.data.RData")

rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL 
load(file = paste0("C:\\ExtremeAttribution2\\CFSR\\TP\\region.coordinations.CFSR.RData"))
lat <- list(west.cor[,2],central.cor[,2],east.cor[,2])
load(file = paste0("C:\\ExtremeAttribution2\\CFSR\\TP\\seasonal.regional.precipitation.data.RData"))
for(j in 1:4){   ######## (j in 1:4)  
  
  for(i in 1:3){  ######### (i in 1:3)
    load(paste("C:\\ExtremeAttribution2\\CFSR\\GPH\\CFSR.hgt_SOM_canada",i,".RData",sep=""))
    pre.data <- t(get(paste0(region[i],".",season[j])))
    # for(ii in 1:nrow(pre.data)){
    #   for(jj in 1:ncol(pre.data)){
    #     if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
    #   }
    # }
    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data[1:length(data),], 
                                                    pattern=data, 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1979,yeare=2016,n.pattern=(k+1)^2,
                                                    weight=T,wgt=cos(lat[[i]]/180*pi) )) 
    } 
  } 
}
save(results,file="CFSR.data.attribution.RData")

##########################################################################################################################
##########################################  MERRA2   data ################################################################
##########################################################################################################################

rm(list=ls())
setwd("C:\\ExtremeAttribution2")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL 
load(file = paste0("C:\\ExtremeAttribution2\\results\\MERRA2\\TP\\region.coordinations.MERRA2.RData"))
lat <- list(west.cor[,2],central.cor[,2],east.cor[,2])
for(j in 1:4){   ######## (j in 1:4) 
  load(file = paste0("C:\\ExtremeAttribution2\\results\\MERRA2\\TP\\Corrected_MERRA2.seasonal.regional.precipitation.data.RData"))
  for(i in 1:3){  ######### (i in 1:3)
    # load(paste("C:\\ExtremeAttribution\\data\\MERRA\\hgt\\hgt_anomalies\\MERRA.hgt_SOM_canada_",i,".RData",sep=""))
    load(paste("C:\\ExtremeAttribution2\\results\\MERRA2\\GPH\\MERRA2.hgt_SOM_canada_",i,".RData",sep=""))
    pre.data <- get(paste0(region[i],".",season[j]))
    # for(ii in 1:nrow(pre.data)){
    #   for(jj in 1:ncol(pre.data)){
    #     if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
    #   }
    # }

    for(k in 1:length(som)){
      if(k==3){
        data <- get(paste0("som_model.hgt.",season[j]))
      }
      else{
        data <- get(paste0("som_model.hgt.",season[j],"_",som[k]))
      }
      data <- data$unit.classif
      results <- rbind (results, winter.attribution(east.winter= pre.data, 
                                                    pattern=data, 
                                                    name=paste0(region[i],"_",season[j],"_som_",k+1),
                                                    season=season[j],years=1980,yeare=2016,n.pattern=(k+1)^2,
                                                    weight=T,wgt=cos(lat[[i]]/180*pi) ))
      rm(som_model.west.spring,som_model.central.spring,som_model.east.spring,
         som_model.west.summer,som_model.central.summer,som_model.east.summer,
         som_model.west.fall,som_model.central.fall,som_model.east.fall,
         som_model.west.winter,som_model.central.winter,som_model.east.winter)
    } 
  }
}
save(results,file="MERRA2.data.attribution.RData")

