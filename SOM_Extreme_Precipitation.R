
rm(list=ls())
setwd("C:\\ExtremeAttribution")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")

library(SDMTools)
library(rgdal)
library(raster)
library(ggmap)
library(rgeos)
library(maptools)
# library(dplyr)
# library(tidyr)
library(png)
library(grid)
library(sp)
## read precipitation data

library(ncdf4)
canada <- readOGR(dsn="F:/data/Canada data/Map/province1","provincell")

year <- 1979:2015
month <- c("01","02","03","04","05","06","07","08","09","10","11","12")
days_noleap = c(31,28,31,30,31,30,31,31,30,31,30,31)
days_cum <- cumsum(days_noleap)
days_leap   = c(31,29,31,30,31,30,31,31,30,31,30,31)
days_cum_leap <- cumsum(days_leap)


west.winter <- central.winter <- east.winter <- west.spring <- central.spring <- east.spring <-
  west.summer <- central.summer <- east.summer <- west.fall <- central.fall <- east.fall <- NULL

for(i_year in 1:length(year)){
  con <- paste("C:\\Users\\Tan\\Documents\\NARR\\apcp.",year[i_year],".nc",sep="")
  file <- nc_open(filename=con)
  #     level <- ncvar_get(file,varid="level")
  #     lat <- ncvar_get(file,varid="lat")
  #     lon <- ncvar_get(file,varid="lon")
  time <- ncvar_get(file,varid="time")
  apcp <- ncvar_get(file,varid="apcp",start=c(1,1,1),count=c(-1,-1,-1))  ######### level = 17 for 500hPa 
  apcp.west <- apcp[1:173,109:dim(apcp)[2],]
  apcp.central <- apcp[163:245,109:dim(apcp)[2],]
  apcp.east <- apcp[182:349,109:dim(apcp)[2],]
  
  data.west <- matrix(NA,length(time),dim(apcp.west)[1]*dim(apcp.west)[2])
  data.central <- matrix(NA,length(time),dim(apcp.central)[1]*dim(apcp.central)[2])
  data.east <- matrix(NA,length(time),dim(apcp.east)[1]*dim(apcp.east)[2])
  for (j in 1:length(time)){
    data.west[j,] <- t(as.vector(apcp.west[,,j]))
    data.central[j,] <- t(as.vector(apcp.central[,,j]))
    data.east[j,] <- t(as.vector(apcp.east[,,j]))
  }
  
  if(year[i_year]%%4 !=0){
      west.winter <- rbind(west.winter,data.west[c(1:days_cum[2],(days_cum[11]+1):days_cum[12]),])
      central.winter <- rbind(central.winter,data.central[c(1:days_cum[2],(days_cum[11]+1):days_cum[12]),])
      east.winter <- rbind(east.winter,data.east[c(1:days_cum[2],(days_cum[11]+1):days_cum[12]),])
  }
  else {
    west.winter <- rbind(west.winter,data.west[c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]),])
    central.winter <- rbind(central.winter,data.central[c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]),])
    east.winter <- rbind(east.winter,data.east[c(1:days_cum_leap[2],(days_cum_leap[11]+1):days_cum_leap[12]),])
  }
 
}
save(west.winter,central.winter,east.winter,file="winter.precipitation.RData")
rm(west.winter,central.winter,east.winter)

for(i_year in 1:length(year)){
  con <- paste("C:\\Users\\Tan\\Documents\\NARR\\apcp.",year[i_year],".nc",sep="")
  file <- nc_open(filename=con)
  #     level <- ncvar_get(file,varid="level")
  #     lat <- ncvar_get(file,varid="lat")
  #     lon <- ncvar_get(file,varid="lon")
  time <- ncvar_get(file,varid="time")
  apcp <- ncvar_get(file,varid="apcp",start=c(1,1,1),count=c(-1,-1,-1))  ######### level = 17 for 500hPa 
  apcp.west <- apcp[1:173,109:dim(apcp)[2],]
  apcp.central <- apcp[163:245,109:dim(apcp)[2],]
  apcp.east <- apcp[182:349,109:dim(apcp)[2],]
  
  data.west <- matrix(NA,length(time),dim(apcp.west)[1]*dim(apcp.west)[2])
  data.central <- matrix(NA,length(time),dim(apcp.central)[1]*dim(apcp.central)[2])
  data.east <- matrix(NA,length(time),dim(apcp.east)[1]*dim(apcp.east)[2])
  for (j in 1:length(time)){
    data.west[j,] <- t(as.vector(apcp.west[,,j]))
    data.central[j,] <- t(as.vector(apcp.central[,,j]))
    data.east[j,] <- t(as.vector(apcp.east[,,j]))
  }
  if(year[i_year]%%4 !=0){
    west.spring <- rbind(west.spring,data.west[c((days_cum[2]+1):days_cum[5]),])
    central.spring <- rbind(central.spring,data.central[c((days_cum[2]+1):days_cum[5]),])
    east.spring <- rbind(east.spring,data.east[c((days_cum[2]+1):days_cum[5]),])
  }
  else {
    west.spring <- rbind(west.spring,data.west[c((days_cum_leap[2]+1):days_cum_leap[5]),])
    central.spring <- rbind(central.spring,data.central[c((days_cum_leap[2]+1):days_cum_leap[5]),])
    east.spring <- rbind(east.spring,data.east[c((days_cum_leap[2]+1):days_cum_leap[5]),])
  }
  
}
save(west.spring,central.spring,east.spring,file="spring.precipitation.RData")
rm(west.spring,central.spring,east.spring)

for(i_year in 1:length(year)){
  con <- paste("C:\\Users\\Tan\\Documents\\NARR\\apcp.",year[i_year],".nc",sep="")
  file <- nc_open(filename=con)
  #     level <- ncvar_get(file,varid="level")
  #     lat <- ncvar_get(file,varid="lat")
  #     lon <- ncvar_get(file,varid="lon")
  
  
  time <- ncvar_get(file,varid="time")
  apcp <- ncvar_get(file,varid="apcp",start=c(1,1,1),count=c(-1,-1,-1))  ######### level = 17 for 500hPa 
  apcp.west <- apcp[1:173,109:dim(apcp)[2],]
  apcp.central <- apcp[163:245,109:dim(apcp)[2],]
  apcp.east <- apcp[182:349,109:dim(apcp)[2],]
  
  data.west <- matrix(NA,length(time),dim(apcp.west)[1]*dim(apcp.west)[2])
  data.central <- matrix(NA,length(time),dim(apcp.central)[1]*dim(apcp.central)[2])
  data.east <- matrix(NA,length(time),dim(apcp.east)[1]*dim(apcp.east)[2])
  for (j in 1:length(time)){
    data.west[j,] <- t(as.vector(apcp.west[,,j]))
    data.central[j,] <- t(as.vector(apcp.central[,,j]))
    data.east[j,] <- t(as.vector(apcp.east[,,j]))
  }
  if(year[i_year]%%4 !=0){
    west.summer <- rbind(west.summer,data.west[c((days_cum[5]+1):days_cum[8]),])
    central.summer <- rbind(central.summer,data.central[c((days_cum[5]+1):days_cum[8]),])
    east.summer <- rbind(east.summer,data.east[c((days_cum[5]+1):days_cum[8]),])
  }
  else {
    west.summer <- rbind(west.summer,data.west[c((days_cum_leap[5]+1):days_cum_leap[8]),])
    central.summer <- rbind(central.summer,data.central[c((days_cum_leap[5]+1):days_cum_leap[8]),])
    east.summer <- rbind(east.summer,data.east[c((days_cum_leap[5]+1):days_cum_leap[8]),])
  }
  
}
save(west.summer,central.summer,east.summer,file="summer.precipitation.RData")
rm(west.summer,central.summer,east.summer)

for(i_year in 1:length(year)){
  con <- paste("C:\\Users\\Tan\\Documents\\NARR\\apcp.",year[i_year],".nc",sep="")
  file <- nc_open(filename=con)
  #     level <- ncvar_get(file,varid="level")
  #     lat <- ncvar_get(file,varid="lat")
  #     lon <- ncvar_get(file,varid="lon")
  time <- ncvar_get(file,varid="time")
  apcp <- ncvar_get(file,varid="apcp",start=c(1,1,1),count=c(-1,-1,-1))  ######### level = 17 for 500hPa 
  apcp.west <- apcp[1:173,109:dim(apcp)[2],]
  apcp.central <- apcp[163:245,109:dim(apcp)[2],]
  apcp.east <- apcp[182:349,109:dim(apcp)[2],]
  
  data.west <- matrix(NA,length(time),dim(apcp.west)[1]*dim(apcp.west)[2])
  data.central <- matrix(NA,length(time),dim(apcp.central)[1]*dim(apcp.central)[2])
  data.east <- matrix(NA,length(time),dim(apcp.east)[1]*dim(apcp.east)[2])
  for (j in 1:length(time)){
    data.west[j,] <- t(as.vector(apcp.west[,,j]))
    data.central[j,] <- t(as.vector(apcp.central[,,j]))
    data.east[j,] <- t(as.vector(apcp.east[,,j]))
  }
  if(year[i_year]%%4 !=0){
    west.fall <- rbind(west.fall,data.west[c((days_cum[8]+1):days_cum[11]),])
    central.fall <- rbind(central.fall,data.central[c((days_cum[8]+1):days_cum[11]),])
    east.fall <- rbind(east.fall,data.east[c((days_cum[8]+1):days_cum[11]),])
  }
  else {
    west.fall <- rbind(west.fall,data.west[c((days_cum_leap[8]+1):days_cum_leap[11]),])
    central.fall <- rbind(central.fall,data.central[c((days_cum_leap[8]+1):days_cum_leap[11]),])
    east.fall <- rbind(east.fall,data.east[c((days_cum_leap[8]+1):days_cum_leap[11]),])
  }
  
}
save(west.fall,central.fall,east.fall,file="fall.precipitation.RData")
rm(west.fall,central.fall,east.fall)

rm(data.central,data.west,data.east)

########## west precipitation ###########

############# SOM linked with extreme precipitation ##########
library(kohonen)
load(file="som_model.east.winter.RData")
pattern <- som_model.east.winter$unit.classif
load(file="winter.precipitation.RData")
rm(central.winter,west.winter)

year.day <- matrix(NA,35,1)
year.day[1,] <- 90  ########## winter just has data from 1979-2013
for(i in 2:35){
  if(year[i+1]%%4 == 0){
    year.day[i,] <- year.day[i-1,] + 91
  }
  else{
    year.day[i,] <- year.day[i-1,] + 90
  }
}

pat.occ <- matrix(NA,16,35)
for(i in 1:16){
  pat.occ[i,1] <- length(pattern[61:year.day[1]][pattern[61:year.day[1]]==i])
  for(j in 2:35){
    pat.occ[i,j] <- length(pattern[year.day[j]:year.day[j-1]][pattern[year.day[j]:year.day[j-1]]==i])
  }
} 
avg.pat.occ <- rowMeans(pat.occ)               ########  fi
time <- 1:35
pat.occ.trend <- matrix(NA, 16,2)
for(i in 1:16){ 
  a <- lm (pat.occ[i,]~time) 
  pat.occ.trend[i,1] <- a$coefficients[2]      ######## delta fi
  pat.occ.trend[i,2] <- summary(a)$coefficients[2,4]
} 


############# mean precipitation ############
for(i in nrow(east.winter)){
  for(j in 1:ncol(east.winter)){
    if (east.winter[i,j]< 0.001 & !is.na(east.winter[i,j])){
      east.winter[i,j]<- 0
    }
  }
}

east.winter.mean <- colMeans(east.winter,na.rm=T)       ########## mean precipitation averaged on all patterns
link <- link.anom <- link.delta <- NULL
precip.pattern <- matrix(NA,16,35)

for(i in 1:16){
  link[[i]] <- east.winter[pattern==i,]                                    ######    pi
  link.delta[[i]] <- (colMeans(link[[i]],na.rm=T) - east.winter.mean)      ######    
  link.anom[[i]] <- (colMeans(link[[i]],na.rm=T) - east.winter.mean)/east.winter.mean*100
  
  precip.pattern [i,1] <- mean(east.winter[pattern[61:year.day[1]]==i,],na.rm=T)
  for(j in 2:35){
    precip.pattern [i,j] <- mean(east.winter[pattern[year.day[j]:year.day[j-1]]==i,],na.rm=T)
  }
}

for(i in 1:16){
  for(j in 1:35){
    if(is.na(precip.pattern [i,j])){
      precip.pattern [i,j] <- 0
    }
  }
}
precip.pattern.mean <- rowMeans(precip.pattern, na.rm=T)

precip.pattern.trend <- matrix(NA, 16,2)
for(i in 1:16){ 
  a <- lm (precip.pattern[i,]~time) 
  precip.pattern.trend[i,1] <- a$coefficients[2]      ######## delta pi
  precip.pattern.trend[i,2] <- summary(a)$coefficients[2,4]
} 

############ contribution of dynamics and thermodynamics to precipitiaton
thermo <- dynamic <- combined <- present <- 0
for(i in 1:16){
  thermo <- thermo + avg.pat.occ[i]*precip.pattern.trend[i,1]     #### fi * delta pi
  dynamic <- dynamic + pat.occ.trend[i,1]*precip.pattern.mean[i]   #### delta fi * pi
  combined <- combined + pat.occ.trend[i,1]*precip.pattern.trend[i,1]  #######  delta pi * delta fi
  present <- present + avg.pat.occ[i]*precip.pattern.mean[i]                #######  fi * pi
}  


############# occurrences of extreme precipitation ###############
# threshold <- apply (X=east.winter,MARGIN = 2, FUN=quantile,probs=c(0.95,0.99),na.rm=T)

threshold <- matrix(NA,2,ncol(east.winter))
for(i in 1:ncol(east.winter)){
  threshold[,i] <- quantile(east.winter[,i][!is.na(east.winter[,i]) & east.winter[,i] > 0.01],probs=c(0.95,0.99),na.rm=T)
}


heavy <- matrix(NA,ncol(east.winter),2)
for(i in 1:ncol(east.winter)){
  heavy[i,1] <- length(east.winter[,i][east.winter[,i]>threshold[1,i]  &  !is.na(east.winter[,i])])
  heavy[i,2] <- length(east.winter[,i][east.winter[,i]>threshold[2,i]  &  !is.na(east.winter[,i])])
}

heavy.occ <- extreme.occ <- matrix(0,16,35)
for(i in 1:16){
  heavy.occ[i,1] <- mean(heavy[,1][61:year.day[1]][pattern[61:year.day[1]]==i],na.rm=T)
  extreme.occ[i,1] <- mean(heavy[,2][61:year.day[1]][pattern[61:year.day[1]]==i],na.rm=T)
  for(j in 2:35){
    heavy.occ[i,j] <- mean(heavy[,1][year.day[j]:year.day[j-1]][pattern[year.day[j]:year.day[j-1]]==i],na.rm=T)
    extreme.occ[i,j] <- mean(heavy[,2][year.day[j]:year.day[j-1]][pattern[year.day[j]:year.day[j-1]]==i],na.rm=T)
  }
}

for(i in 1:16){
  for(j in 1:35){
    if(is.na(heavy.occ[i,j])){
      heavy.occ[i,j] <- 0
    }
    if(is.na(extreme.occ[i,j])){
      extreme.occ[i,j] <- 0
    }
  }
}

heavy.occ.mean <- rowMeans(heavy.occ)
extreme.occ.mean <- rowMeans(extreme.occ)

##########
heavy.occ.trend <- extreme.occ.trend <- matrix(NA,16,2)
for(i in 1:16){
  a <- lm (heavy.occ[i,]~time)
  b <- lm (extreme.occ[i,]~time)
  heavy.occ.trend[i,1] <- a$coefficients[2]      ######## delta fi
  heavy.occ.trend[i,2] <- summary(a)$coefficients[2,4]
  extreme.occ.trend[i,1] <- b$coefficients[2]      ######## delta fi
  extreme.occ.trend[i,2] <- summary(b)$coefficients[2,4]
}

############ contribution of dynamics and thermodynamics to extreme precipitiaton
thermo.heavy <- dynamic.heavy <- combined.heavy <- 0
thermo.extreme <- dynamic.extreme <- combined.extreme <- 0
for(i in 1:16){
  thermo.heavy <- thermo.heavy + avg.pat.occ[i]/91*heavy.occ.trend[i,1]/ncol(east.winter)             #### fi * delta pi
  dynamic.heavy <- dynamic.heavy + pat.occ.trend[i,1]/91*heavy.occ.mean[i]/ncol(east.winter)          #### delta fi * pi
  combined.heavy <- combined.heavy + pat.occ.trend[i,1]/91*heavy.occ.trend[i,1]/ncol(east.winter)     #######  delta pi * delta fi
  thermo.extreme <- thermo.extreme + avg.pat.occ[i]/91*extreme.occ.trend[i,1]/ncol(east.winter)       #### fi * delta pi
  dynamic.extreme <- dynamic.extreme + pat.occ.trend[i,1]/91*extreme.occ.mean[i]/ncol(east.winter)       #### delta fi * pi
  combined.extreme <- combined.extreme + pat.occ.trend[i,1]/91*extreme.occ.trend[i,1]/ncol(east.winter)  #######  delta pi * delta fi
}

############ NARR precipitation #########################
library(kohonen)
load(file="som_model.east.winter.RData")
pattern <- som_model.east.winter$unit.classif
load(file="winter.precipitation.RData")
rm(central.winter,west.winter)


rm(list=ls())
setwd("C:\\ExtremeAttribution")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east")
season <- c("spring","summer","fall","winter")
results <- NULL
for(j in 1:4){
  filename1 <- paste0("C:\\ExtremeAttribution\\data\\NARR\\tp\\",season[j],".precipitation.RData")
  file1 <- load(filename1)
  for(i in 1:3){
    filename <- paste0("C:\\ExtremeAttribution\\data\\NARR\\hgt\\som_model.",region[i],".",
                       season[j],".RData")
    file <- load(filename)
    data <- get(file)
    results[[(j-1)*3+i]] <- winter.attribution(east.winter=t(get(paste0(region[i],".",season[j]))), 
                                               pattern=data$unit.classif, 
                                               name=paste0(region[i],"_",season[j],"_"),
                                               season=season[j])
    rm(som_model.west.winter,som_model.west.spring,som_model.west.summer,som_model.west.fall,
       som_model.central.winter,som_model.central.spring,som_model.central.summer,som_model.central.fall,
       som_model.east.winter,som_model.east.spring,som_model.east.summer,som_model.east.fall)
  }
}

# save(results.west.spring,results.central.spring,results.east.spring,
#      results.west.summer,results.central.summer,results.east.summer,
#      results.west.fall,results.central.fall,results.east.fall,
#      results.west.winter,results.central.winter,results.east.winter,
#      file="results.attribution.absolute_value.RData")

save(results,file="attribution.results.NARR.RData")

##################### NARR precipitation with NCEP reanalysis #2 geopotential heights  #############
region <- c("west","central","east")
season <- c("spring","summer","fall","winter")
results <- NULL
for(j in 1:4){
  filename1 <- paste0(season[j],".precipitation.RData")
  file1 <- load(filename1)
  for(i in 1:3){
    load(paste("NCEP2.hgt_SOM_",i,".RData",sep=""))
    data <- get(paste0("som_model.hgt.",season[j]))
    data <- data$unit.classif
    results[[(j-1)*3+i]] <- winter.attribution(east.winter=t(get(paste0(region[i],".",season[j]))), 
                                               pattern=data, 
                                               name=paste0(region[i],"_",season[j],"_ncep2_narr"),
                                               season=season[j])
  }
  rm(west.winter,central.winter,east.winter,
     west.spring,central.spring,east.spring,
     west.summer,central.summer,east.summer,
     west.fall,central.fall,east.fall)
}

save(results,file="attribution.results.NCEP2_NARR.RData")

########################## observed data ########################


region <- c("west","central","east")
season <- c("spring","summer","fall","winter")
results <- NULL
for(j in 1:4){   ######## (j in 1:4)
  for(i in 1:3){  ######### (i in 1:3)
    filename <- paste0("som_model.",region[i],".",
                       season[j],".RData")
    file <- load(filename)
    data <- get(paste0("som_model.hgt.",season[j]))
    load(file="C:\\ExtremeAttribution\\precipitation.data.RData")
    east.winter <- get(paste0(region[i],".",season[j]))
    results[[(j-1)*3+i]] <- winter.attribution(east.winter= east.winter,
                                              pattern=data$unit.classif[1:nrow(east.winter)], 
                                              name=paste0(region[i],"_",season[j],"_observed"),
                                              season=season[j],yeare=2005)
    rm(som_model.west.winter,som_model.west.spring,som_model.west.summer,som_model.west.fall,
       som_model.central.winter,som_model.central.spring,som_model.central.summer,som_model.central.fall,
       som_model.east.winter,som_model.east.spring,som_model.east.summer,som_model.east.fall)
    
  }
}
save(results,file="observed.data.attribution.RData")


######################## Mckenney Precipitation dataset ####################
rm(list=ls())
year <- 1979:2013
month <- c("01","02","03","04","05","06","07","08","09","10","11","12")
days_noleap = c(31,28,31,30,31,30,31,31,30,31,30,31)
days_cum <- cumsum(days_noleap)
days_leap   = c(31,29,31,30,31,30,31,31,30,31,30,31)
days_cum_leap <- cumsum(days_leap)

west.spring <- west.summer <- west.fall <- west.winter <-
  central.spring <- central.summer <- central.fall <- central.winter <- 
  east.spring <- east.summer <- east.fall <- east.winter <- NULL
for(i_year in 1:length(year)){
  load(file = paste0("C:\\ExtremeAttribution\\Mcken\\precipitation.",year[i_year],".RData"))
  if(year[i_year]%% 4 == 0){
    west.spring <- data.west[(days_cum_leap[2]+1):(days_cum_leap[5]),]
    west.summer <- data.west[(days_cum_leap[5]+1):(days_cum_leap[8]),]
    west.fall <- data.west[(days_cum_leap[8]+1):(days_cum_leap[11]),] 
    west.winter <- data.west[c(1:days_cum_leap[2],(days_cum_leap[11]+1):(days_cum_leap[12])),]
    central.spring <- data.central[(days_cum_leap[2]+1):(days_cum_leap[5]),]
    central.summer <- data.central[(days_cum_leap[5]+1):(days_cum_leap[8]),]
    central.fall <- data.central[(days_cum_leap[8]+1):(days_cum_leap[11]),] 
    central.winter <- data.central[c(1:days_cum_leap[2],(days_cum_leap[11]+1):(days_cum_leap[12])),]
    east.spring <- data.east[(days_cum_leap[2]+1):(days_cum_leap[5]),]
    east.summer <- data.east[(days_cum_leap[5]+1):(days_cum_leap[8]),]
    east.fall <- data.east[(days_cum_leap[8]+1):(days_cum_leap[11]),] 
    east.winter <- data.east[c(1:days_cum_leap[2],(days_cum_leap[11]+1):(days_cum_leap[12])),]
  }
  else{
    east.spring <- data.east[(days_cum[2]+1):(days_cum[5]),]
    east.summer <- data.east[(days_cum[5]+1):(days_cum[8]),]
    east.fall <- data.east[(days_cum[8]+1):(days_cum[11]),]
    east.winter <- data.east[c(1:days_cum[2],(days_cum[11]+1):(days_cum[12])),]
    west.spring <- data.west[(days_cum[2]+1):(days_cum[5]),]
    west.summer <- data.west[(days_cum[5]+1):(days_cum[8]),]
    west.fall <- data.west[(days_cum[8]+1):(days_cum[11]),]
    west.winter <- data.west[c(1:days_cum[2],(days_cum[11]+1):(days_cum[12])),]
    central.spring <- data.central[(days_cum[2]+1):(days_cum[5]),]
    central.summer <- data.central[(days_cum[5]+1):(days_cum[8]),]
    central.fall <- data.central[(days_cum[8]+1):(days_cum[11]),]
    central.winter <- data.central[c(1:days_cum[2],(days_cum[11]+1):(days_cum[12])),]
  }
  save(west.spring,central.spring,east.spring,file=paste0("spring.precipitation.mck.",year[i_year],".RData"))
  save(west.summer,central.summer,east.summer,file=paste0("summer.precipitation.mck.",year[i_year],".RData"))
  save(west.fall,central.fall,east.fall,file=paste0("fall.precipitation.mck.",year[i_year],".RData"))
  save(west.winter,central.winter,east.winter,file=paste0("winter.precipitation.mck.",year[i_year],".RData"))
}

west.spring1 <- west.summer1 <- west.fall1 <- west.winter1 <-
  central.spring1 <- central.summer1 <- central.fall1 <- central.winter1 <- 
  east.spring1 <- east.summer1 <- east.fall1 <- east.winter1 <- NULL

for(i_year in 1:length(year)){
  load(file=paste0("spring.precipitation.mck.",year[i_year],".RData"))
  # west.spring1 <- rbind(west.spring1,west.spring)
  central.spring1 <- rbind(central.spring1,central.spring)
   # east.spring1 <- rbind(east.spring1,east.spring)
}
save(west.spring1,file="west.spring.mck.RData")
rm(west.spring1)
save(central.spring1,file="central.spring.mck.RData")
rm(central.spring1)
save(east.spring1,file="east.spring.mck.RData")
rm(east.spring1)

for(i_year in 1:length(year)){
  load(file=paste0("winter.precipitation.mck.",year[i_year],".RData"))
  west.winter1 <- rbind(west.winter1,west.winter)
  # central.winter1 <- rbind(central.winter1,central.winter)
#   east.winter1 <- rbind(east.winter1,east.winter)
}
save(west.winter1,file="west.winter.mck.RData")
rm(west.winter1)
save(central.winter1,file="central.winter.mck.RData")
rm(central.winter1)
save(east.winter1,file="east.winter.mck.RData")
rm(east.winter1)

for(i_year in 1:length(year)){
  load(file=paste0("summer.precipitation.mck.",year[i_year],".RData"))
  # west.summer1 <- rbind(west.summer1,west.summer)
  # central.summer1 <- rbind(central.summer1,central.summer)
  east.summer1 <- rbind(east.summer1,east.summer)
}
save(west.summer1,file="west.summer.mck.RData")
rm(west.summer1)
save(central.summer1,file="central.summer.mck.RData")
rm(central.summer1)
save(east.summer1,file="east.summer.mck.RData")
rm(east.summer1)

for(i_year in 1:length(year)){
  load(file=paste0("fall.precipitation.mck.",year[i_year],".RData"))
  # west.fall1 <- rbind(west.fall1,west.fall)
  central.fall1 <- rbind(central.fall1,central.fall)
#   east.fall1 <- rbind(east.fall1,east.fall)
}
save(west.fall1,file="west.fall.mck.RData")
rm(west.fall1)
save(central.fall1,file="central.fall.mck.RData")
rm(central.fall1)
save(east.fall1,file="east.fall.mck.RData")
rm(east.fall1)

################################ attribution McKenney data ##############################
#########################################################################################
rm(list=ls())

rm(list=ls())
setwd("C:\\ExtremeAttribution")
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
region <- c("west","central","east") 
season <- c("spring","summer","fall","winter")
som <- c("2_2","3_3","")
results <- NULL
# load(file = paste0("C:\\ExtremeAttribution\\data\\NARR\\tp\\region.row&col.narr.RData"))
# lat <- list(west.position[,2],central.position[,2],east.position[,2])
for(j in 1:4){   ######## (j in 1:4) 
  load(file = paste0("C:\\ExtremeAttribution\\data\\NARR\\tp\\",season[j],".precipitation.RData"))
  for(i in 1:3){  ######### (i in 1:3)
    pre.data <- get(paste0(region[i],".",season[j]))
    for(ii in 1:nrow(pre.data)){
      for(jj in 1:ncol(pre.data)){
        if(pre.data[ii,jj]<0 | is.na(pre.data[ii,jj])) pre.data[ii,jj] <- 0
      }
    }
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
                                                    season=season[j],years=1979,yeare=2013,n.pattern=(k+1)^2,
                                                    weight=F ))
      rm(som_model.west.spring,som_model.central.spring,som_model.east.spring,
         som_model.west.summer,som_model.central.summer,som_model.east.summer,
         som_model.west.fall,som_model.central.fall,som_model.east.fall,
         som_model.west.winter,som_model.central.winter,som_model.east.winter)
    } 
  }
  rm(west.spring,central.spring,east.spring,
     west.summer,central.summer,east.summer,
     west.fall,central.fall,east.fall,
     west.winter,central.winter,east.winter)
}
save(results,file="NARR.data.attribution.RData")


source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\west.spring.mck.RData")
load(file= paste0("som_model.west.spring.RData"))
west.spirng <- winter.attribution(east.winter=west.spring1,
                       pattern=som_model.west.spring$unit.classif[1:3220], 
                       name="west.spring.mck",
                       season="spring",yeare=2013)
save(west.spring,file="attribution.west.spring.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\central.spring.mck.RData")
load(file= paste0("som_model.central.spring.RData"))
central.spirng <- winter.attribution(east.winter=central.spring1,
                                  pattern=som_model.central.spring$unit.classif[1:3220], 
                                  name="central.spring.mck",
                                  season="spring",yeare=2013)
save(central.spring,file="attribution.central.spring.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\east.spring.mck.RData")
load(file= paste0("som_model.east.spring.RData"))
east.spirng <- winter.attribution(east.winter=east.spring1,
                                  pattern=som_model.east.spring$unit.classif[1:3220], 
                                  name="east.spring.mck",
                                  season="spring",yeare=2013)
save(east.spring,file="attribution.east.spring.mck.RData")

source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\west.summer.mck.RData")
load(file= paste0("som_model.west.summer.RData"))
west.spirng <- winter.attribution(east.winter=west.summer1,
                                  pattern=som_model.west.summer$unit.classif[1:3220], 
                                  name="west.summer.mck",
                                  season="summer",yeare=2013)
save(west.summer,file="attribution.west.summer.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\central.summer.mck.RData")
load(file= paste0("som_model.central.summer.RData"))
central.spirng <- winter.attribution(east.winter=central.summer1,
                                     pattern=som_model.central.summer$unit.classif[1:3220], 
                                     name="central.summer.mck",
                                     season="summer",yeare=2013)
save(central.summer,file="attribution.central.summer.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\east.summer.mck.RData")
load(file= paste0("som_model.east.summer.RData"))
east.spirng <- winter.attribution(east.winter=east.summer1,
                                  pattern=som_model.east.summer$unit.classif[1:3220], 
                                  name="east.summer.mck",
                                  season="summer",yeare=2013)
save(east.summer,file="attribution.east.summer.mck.RData")

source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\west.fall.mck.RData")
load(file= paste0("som_model.west.fall.RData"))
west.spirng <- winter.attribution(east.winter=west.fall1,
                                  pattern=som_model.west.fall$unit.classif[1:3220], 
                                  name="west.fall.mck",
                                  season="fall",yeare=2013)
save(west.fall,file="attribution.west.fall.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\central.fall.mck.RData")
load(file= paste0("som_model.central.fall.RData"))
central.spirng <- winter.attribution(east.winter=central.fall1,
                                     pattern=som_model.central.fall$unit.classif[1:3220], 
                                     name="central.fall.mck",
                                     season="fall",yeare=2013)
save(central.fall,file="attribution.central.fall.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\east.fall.mck.RData")
load(file= paste0("som_model.east.fall.RData"))
east.spirng <- winter.attribution(east.winter=east.fall1,
                                  pattern=som_model.east.fall$unit.classif[1:3220], 
                                  name="east.fall.mck",
                                  season="fall",yeare=2013)
save(east.fall,file="attribution.east.fall.mck.RData")

source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\west.winter.mck.RData")
load(file= paste0("som_model.west.winter.RData"))
west.spirng <- winter.attribution(east.winter=west.winter1,
                                  pattern=som_model.west.winter$unit.classif[1:3220], 
                                  name="west.winter.mck",
                                  season="winter",yeare=2013)
save(west.winter,file="attribution.west.winter.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\central.winter.mck.RData")
load(file= paste0("som_model.central.winter.RData"))
central.spirng <- winter.attribution(east.winter=central.winter1,
                                     pattern=som_model.central.winter$unit.classif[1:3220], 
                                     name="central.winter.mck",
                                     season="winter",yeare=2013)
save(central.winter,file="attribution.central.winter.mck.RData")

rm(list=ls())
source("C:\\Users\\Tan\\Dropbox\\Thermodynamic and dynamic precipitation extreme change\\Code\\DynamicAttribution.R")
load("C:\\ExtremeAttribution\\data\\NARR\\tp\\east.winter.mck.RData")
load(file= paste0("som_model.east.winter.RData"))
east.spirng <- winter.attribution(east.winter=east.winter1,
                                  pattern=som_model.east.winter$unit.classif[1:3220], 
                                  name="east.winter.mck",
                                  season="winter",yeare=2013)
save(east.winter,file="attribution.east.winter.mck.RData")
