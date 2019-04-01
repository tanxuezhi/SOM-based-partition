

winter.attribution <- function(east.winter,pattern, name, season="winter",
                               years=1979,yeare=2014,n.pattern=16,
                               wgt=0,weight=F){
  
  year <- years:yeare
  
  if(season=="winter"){
    year.day <- matrix(NA,length(year),1)
    if(years%%4==0){
      year.day[1,] <- 91  ########## winter just has data from 1979-2013
    }
    else{
      year.day[1,] <- 90  ########## winter just has data from 1979-2013
    }
    for(i in 2:length(year)){
      if(year[i]%%4 == 0){
        year.day[i,] <- year.day[i-1,] + 91
      }
      else{
        year.day[i,] <- year.day[i-1,] + 90
      }
    }
  }
  
  if(season=="summer" | season=="spring")  year.day <- seq(92,92*length(year),92)
  if(season=="fall") year.day <- seq(91,91*length(year),91)

  
  pat.occ <- matrix(NA,n.pattern,length(year))
  for(i in 1:n.pattern){
    pat.occ[i,1] <- length(pattern[1:year.day[1]][pattern[1:year.day[1]]==i])
    for(j in 2:length(year)){
      pat.occ[i,j] <- length(pattern[(year.day[j-1]+1):year.day[j]][pattern[(year.day[j-1]+1):year.day[j]]==i])
    }
  } 
  avg.pat.occ <- rowMeans(pat.occ)               ########  fi
  time <- 1:length(year)
  pat.occ.trend <- matrix(NA, n.pattern,2)
  for(i in 1:n.pattern){ 
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
  
  if(weight){
    east.winter.mean <- colMeans(east.winter,na.rm=T)       ########## mean precipitation averaged on all patterns
    link <- link.anom <- NULL
    link.delta <- NULL
    precip.pattern <- matrix(NA,n.pattern,length(year))
    
    for(i in 1:n.pattern){
      link[[i]] <- east.winter[pattern==i,]                                    ######    pi
      # link.delta[[i]] <- (colMeans(link[[i]],na.rm=T) - east.winter.mean)      ######    
      link.anom[[i]] <- (colMeans(link[[i]],na.rm=T) - east.winter.mean)/east.winter.mean*100
      
      dddd <- colMeans(east.winter[pattern[1:year.day[1]]==i,],na.rm=T)  
      precip.pattern [i,1] <- weighted.mean(x=dddd, w=wgt)          ####### regional mean value
      for(j in 2:length(year)){
        dddd <- colMeans(east.winter[pattern[(year.day[j-1]+1):year.day[j]]==i,],na.rm=T)
        precip.pattern [i,j] <- weighted.mean(x=dddd, w=wgt)        ####### regional mean value  
      }
    }
    rm(link)
    
    
    for(i in 1:n.pattern){
      for(j in 1:length(year)){
        if(is.na(precip.pattern [i,j])){
          precip.pattern [i,j] <- 0
        }
      }
    }
    precip.pattern.mean <- rowMeans(precip.pattern, na.rm=T)
    
    precip.pattern.trend <- matrix(NA, n.pattern,2)
    for(i in 1:n.pattern){ 
      a <- lm (precip.pattern[i,]~time) 
      precip.pattern.trend[i,1] <- a$coefficients[2]      ######## delta pi
      precip.pattern.trend[i,2] <- summary(a)$coefficients[2,4]
    } 
    
    save(avg.pat.occ, 
         pat.occ, pat.occ.trend, 
         link.anom, 
         precip.pattern, precip.pattern.trend,
         file=paste0("link.anom.",name,".RData"))
    
    
    ############ contribution of dynamics and thermodynamics to precipitiaton
    thermo <- dynamic <- combined <- present <- 0
    for(i in 1:n.pattern){
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
    
    heavy.sp <- extreme.sp <- matrix(0,ncol(east.winter),n.pattern)
    for(i in 1:ncol(east.winter)){
      for(j in 1:n.pattern){
        heavy.sp[i,j] <- length(east.winter[,i][east.winter[,i]>threshold[1,i]  &  !is.na(east.winter[,i]) & pattern==j ])
        extreme.sp[i,j] <- length(east.winter[,i][east.winter[,i]>threshold[2,i]  &  !is.na(east.winter[,i]) & pattern==j])
      }
    }
    
    
    heavy.occ <- extreme.occ <- matrix(0,n.pattern,length(year))
    occ.heavy <- occ.extreme <- array(0,dim=c(n.pattern,ncol(east.winter),length(year)))
    for(i in 1:n.pattern){
      
      if(pat.occ[i,1]>1){
        aaa <- east.winter[1:year.day[1],] 
        aaa <- aaa[pattern[1:year.day[1]]==i,]
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[1,ii]  &  !is.na(aaa[,ii])])
        }
        occ.heavy[i,,1] <- heavy.year
        heavy.occ[i,1] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[2,ii]  &  !is.na(aaa[,ii])])
        }
        occ.extreme[i,,1] <- heavy.year
        extreme.occ[i,1] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
      }
      else if (pat.occ[i,1]>0){
        aaa <- east.winter[1:year.day[1],] 
        aaa <- aaa[pattern[1:year.day[1]]==i,]
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[1,ii]  &  !is.na(aaa[ii])])
        }
        occ.heavy[i,,1] <- heavy.year
        heavy.occ[i,1] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[2,ii]  &  !is.na(aaa[ii])])
        }
        occ.extreme[i,,1] <- heavy.year
        extreme.occ[i,1] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
      }
      
      
      # heavy.occ[i,1] <- mean(heavy[,1][1:year.day[1]][pattern[1:year.day[1]]==i],na.rm=T)
      # extreme.occ[i,1] <- mean(heavy[,2][1:year.day[1]][pattern[1:year.day[1]]==i],na.rm=T)
      
      for(j in 2:length(year)){
        
        if(pat.occ[i,j]>1){
          if(year.day[j] > nrow(east.winter)) year.day[j] <- nrow(east.winter)  ########## used for data that do not have whole winter data in MERRA
          aaa <- east.winter[(year.day[j-1]+1):year.day[j],] 
          aaa <- aaa[pattern[(year.day[j-1]+1):year.day[j]]==i,]         
          
          bbb <- (year.day[j-1]+1):year.day[j]
          if(length(bbb[pattern[(year.day[j-1]+1):year.day[j]]==i])<2){        ############ used for data have diffent length between GPH and Precipitation
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[1,ii]  &  !is.na(aaa[ii])])
            }
            occ.heavy[i,,j] <- heavy.year
            heavy.occ[i,j] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T) 
            
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[2,ii]  &  !is.na(aaa[ii])])
            }
            occ.extreme[i,,j] <- heavy.year
            extreme.occ[i,j] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
          }
          else{
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[1,ii]  &  !is.na(aaa[,ii])])
            }
            occ.heavy[i,,j] <- heavy.year
            heavy.occ[i,j] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T) 
            
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[2,ii]  &  !is.na(aaa[,ii])])
            }
            occ.extreme[i,,j] <- heavy.year
            extreme.occ[i,j] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
          }
          
        }
        else if (pat.occ[i,j]>0){
          if(year.day[j] > nrow(east.winter)) year.day[j] <- nrow(east.winter)
          aaa <- east.winter[(year.day[j-1]+1):year.day[j],] 
          aaa <- aaa[pattern[(year.day[j-1]+1):year.day[j]]==i,]         
          
          heavy.year <- matrix(0,ncol(east.winter),1)
          for(ii in 1:ncol(east.winter)){
            heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[1,ii]  &  !is.na(aaa[ii])])
          }
          occ.heavy[i,,j] <- heavy.year
          heavy.occ[i,j] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T) 
          
          heavy.year <- matrix(0,ncol(east.winter),1)
          for(ii in 1:ncol(east.winter)){
            heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[2,ii]  &  !is.na(aaa[ii])])
          }
          occ.extreme[i,,j] <- heavy.year
          extreme.occ[i,j] <- weighted.mean(x=heavy.year,w=wgt,na.rm=T)
        }
        
       }
    }
    
    rm(heavy.year)
    
    for(i in 1:n.pattern){
      for(j in 1:length(year)){
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
    heavy.occ.trend <- extreme.occ.trend <- matrix(NA,n.pattern,2)
    for(i in 1:n.pattern){
      a <- lm (heavy.occ[i,]~time)
      b <- lm (extreme.occ[i,]~time)
      heavy.occ.trend[i,1] <- a$coefficients[2]      ######## delta fi
      heavy.occ.trend[i,2] <- summary(a)$coefficients[2,4]
      extreme.occ.trend[i,1] <- b$coefficients[2]      ######## delta fi
      extreme.occ.trend[i,2] <- summary(b)$coefficients[2,4]
    }
    
    save(heavy, heavy.sp, extreme.sp,occ.heavy,occ.extreme,
         heavy.occ, heavy.occ.trend,
         extreme.occ, extreme.occ.trend,
         file=paste0("heavy.precipitation.occurrence.",name,".RData") )
    
    ############ contribution of dynamics and thermodynamics to extreme precipitiaton
    thermo.heavy <- dynamic.heavy <- combined.heavy <- 0
    thermo.extreme <- dynamic.extreme <- combined.extreme <- 0
    for(i in 1:n.pattern){
      thermo.heavy <- thermo.heavy + avg.pat.occ[i]/91*heavy.occ.trend[i,1]/ncol(east.winter)             #### fi * delta pi
      dynamic.heavy <- dynamic.heavy + pat.occ.trend[i,1]/91*heavy.occ.mean[i]/ncol(east.winter)          #### delta fi * pi
      combined.heavy <- combined.heavy + pat.occ.trend[i,1]/91*heavy.occ.trend[i,1]/ncol(east.winter)     #######  delta pi * delta fi
      thermo.extreme <- thermo.extreme + avg.pat.occ[i]/91*extreme.occ.trend[i,1]/ncol(east.winter)       #### fi * delta pi
      dynamic.extreme <- dynamic.extreme + pat.occ.trend[i,1]/91*extreme.occ.mean[i]/ncol(east.winter)       #### delta fi * pi
      combined.extreme <- combined.extreme + pat.occ.trend[i,1]/91*extreme.occ.trend[i,1]/ncol(east.winter)  #######  delta pi * delta fi
    }
    results <- c( thermo, dynamic,combined, present, thermo.heavy, dynamic.heavy, combined.heavy,
                  thermo.extreme,dynamic.extreme,combined.extreme)
  }
  
  else{
    east.winter.mean <- colMeans(east.winter,na.rm=T)       ########## mean precipitation averaged on all patterns
    link <- link.anom <- NULL
    link.delta <- NULL
    precip.pattern <- matrix(NA,n.pattern,length(year))
    
    for(i in 1:n.pattern){
      link[[i]] <- east.winter[pattern==i,]                                    ######    pi
      # link.delta[[i]] <- (colMeans(link[[i]],na.rm=T) - east.winter.mean)      ######    
      link.anom[[i]] <- (colMeans(link[[i]],na.rm=T) - east.winter.mean)/east.winter.mean*100
      
      precip.pattern [i,1] <- mean(east.winter[pattern[1:year.day[1]]==i,],na.rm=T)
      for(j in 2:length(year)){
        precip.pattern [i,j] <- mean(east.winter[pattern[(year.day[j-1]+1):year.day[j]]==i,],na.rm=T)
      }
    }
    rm(link)
    
    
    for(i in 1:n.pattern){
      for(j in 1:length(year)){
        if(is.na(precip.pattern [i,j])){
          precip.pattern [i,j] <- 0
        }
      }
    }
    precip.pattern.mean <- rowMeans(precip.pattern, na.rm=T)
    
    precip.pattern.trend <- matrix(NA, n.pattern,2)
    for(i in 1:n.pattern){ 
      a <- lm (precip.pattern[i,]~time) 
      precip.pattern.trend[i,1] <- a$coefficients[2]      ######## delta pi
      precip.pattern.trend[i,2] <- summary(a)$coefficients[2,4]
    } 
    
    save(avg.pat.occ, 
         pat.occ, pat.occ.trend, 
         link.anom, 
         precip.pattern, precip.pattern.trend,
         file=paste0("link.anom.",name,".RData"))
    
    
    ############ contribution of dynamics and thermodynamics to precipitiaton
    thermo <- dynamic <- combined <- present <- 0
    for(i in 1:n.pattern){
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
    
    heavy.sp <- extreme.sp <- matrix(0,ncol(east.winter),n.pattern)
    for(i in 1:ncol(east.winter)){
      for(j in 1:n.pattern){
        heavy.sp[i,j] <- length(east.winter[,i][east.winter[,i]>threshold[1,i]  &  !is.na(east.winter[,i]) & pattern==j ])
        extreme.sp[i,j] <- length(east.winter[,i][east.winter[,i]>threshold[2,i]  &  !is.na(east.winter[,i]) & pattern==j])
      }
    }
    
    
    heavy.occ <- extreme.occ <- matrix(0,n.pattern,length(year))
    for(i in 1:n.pattern){
      
      if(pat.occ[i,1]>1){
        aaa <- east.winter[1:year.day[1],] 
        aaa <- aaa[pattern[1:year.day[1]]==i,]
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[1,ii]  &  !is.na(aaa[,ii])])
        }
        heavy.occ[i,1] <- mean(heavy.year,na.rm=T)
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[2,ii]  &  !is.na(aaa[,ii])])
        }
        extreme.occ[i,1] <- mean(heavy.year,na.rm=T)
      }
      else if (pat.occ[i,1]>0){
        aaa <- east.winter[1:year.day[1],] 
        aaa <- aaa[pattern[1:year.day[1]]==i,]
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[1,ii]  &  !is.na(aaa[ii])])
        }
        heavy.occ[i,1] <- mean(heavy.year,na.rm=T)
        
        heavy.year <- matrix(0,ncol(east.winter),1)
        for(ii in 1:ncol(east.winter)){
          heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[2,ii]  &  !is.na(aaa[ii])])
        }
        extreme.occ[i,1] <- mean(heavy.year,na.rm=T)
      }
      
      
      # heavy.occ[i,1] <- mean(heavy[,1][1:year.day[1]][pattern[1:year.day[1]]==i],na.rm=T)
      # extreme.occ[i,1] <- mean(heavy[,2][1:year.day[1]][pattern[1:year.day[1]]==i],na.rm=T)
      
      for(j in 2:length(year)){
        
        if(pat.occ[i,j]>1){
          if(year.day[j] > nrow(east.winter)) year.day[j] <- nrow(east.winter)  ########## used for data that do not have whole winter data in MERRA
          aaa <- east.winter[(year.day[j-1]+1):year.day[j],] 
          aaa <- aaa[pattern[(year.day[j-1]+1):year.day[j]]==i,]         
          
          bbb <- (year.day[j-1]+1):year.day[j]
          if(length(bbb[pattern[(year.day[j-1]+1):year.day[j]]==i])<2){        ############ used for data have diffent length between GPH and Precipitation
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[1,ii]  &  !is.na(aaa[ii])])
            }
            heavy.occ[i,j] <- mean(heavy.year,na.rm=T) 
            
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[2,ii]  &  !is.na(aaa[ii])])
            }
            extreme.occ[i,j] <- mean(heavy.year,na.rm=T)
          }
          else{
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[1,ii]  &  !is.na(aaa[,ii])])
            }
            heavy.occ[i,j] <- mean(heavy.year,na.rm=T) 
            
            heavy.year <- matrix(0,ncol(east.winter),1)
            for(ii in 1:ncol(east.winter)){
              heavy.year[ii] <- length(aaa[,ii][aaa[,ii]>threshold[2,ii]  &  !is.na(aaa[,ii])])
            }
            extreme.occ[i,j] <- mean(heavy.year,na.rm=T)
          }
          
        }
        else if (pat.occ[i,j]>0){
          if(year.day[j] > nrow(east.winter)) year.day[j] <- nrow(east.winter)
          aaa <- east.winter[(year.day[j-1]+1):year.day[j],] 
          aaa <- aaa[pattern[(year.day[j-1]+1):year.day[j]]==i,]         
          
          heavy.year <- matrix(0,ncol(east.winter),1)
          for(ii in 1:ncol(east.winter)){
            heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[1,ii]  &  !is.na(aaa[ii])])
          }
          heavy.occ[i,j] <- mean(heavy.year,na.rm=T) 
          
          heavy.year <- matrix(0,ncol(east.winter),1)
          for(ii in 1:ncol(east.winter)){
            heavy.year[ii] <- length(aaa[ii][aaa[ii]>threshold[2,ii]  &  !is.na(aaa[ii])])
          }
          extreme.occ[i,j] <- mean(heavy.year,na.rm=T)
        }
        
        # heavy.occ[i,j]   <- mean(heavy[,1][(year.day[j-1]+1):year.day[j]][pattern[(year.day[j-1]+1):year.day[j]]==i],na.rm=T)
        # extreme.occ[i,j] <- mean(heavy[,2][(year.day[j-1]+1):year.day[j]][pattern[(year.day[j-1]+1):year.day[j]]==i],na.rm=T)
      }
    }
    
    rm(heavy.year)
    
    for(i in 1:n.pattern){
      for(j in 1:length(year)){
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
    heavy.occ.trend <- extreme.occ.trend <- matrix(NA,n.pattern,2)
    for(i in 1:n.pattern){
      a <- lm (heavy.occ[i,]~time)
      b <- lm (extreme.occ[i,]~time)
      heavy.occ.trend[i,1] <- a$coefficients[2]      ######## delta fi
      heavy.occ.trend[i,2] <- summary(a)$coefficients[2,4]
      extreme.occ.trend[i,1] <- b$coefficients[2]      ######## delta fi
      extreme.occ.trend[i,2] <- summary(b)$coefficients[2,4]
    }
    
    save(heavy, heavy.sp, extreme.sp,
         heavy.occ, heavy.occ.trend,
         extreme.occ, extreme.occ.trend,
         file=paste0("heavy.precipitation.occurrence.",name,".RData") )
    
    ############ contribution of dynamics and thermodynamics to extreme precipitiaton
    thermo.heavy <- dynamic.heavy <- combined.heavy <- 0
    thermo.extreme <- dynamic.extreme <- combined.extreme <- 0
    for(i in 1:n.pattern){
      thermo.heavy <- thermo.heavy + avg.pat.occ[i]/91*heavy.occ.trend[i,1]/ncol(east.winter)             #### fi * delta pi
      dynamic.heavy <- dynamic.heavy + pat.occ.trend[i,1]/91*heavy.occ.mean[i]/ncol(east.winter)          #### delta fi * pi
      combined.heavy <- combined.heavy + pat.occ.trend[i,1]/91*heavy.occ.trend[i,1]/ncol(east.winter)     #######  delta pi * delta fi
      thermo.extreme <- thermo.extreme + avg.pat.occ[i]/91*extreme.occ.trend[i,1]/ncol(east.winter)       #### fi * delta pi
      dynamic.extreme <- dynamic.extreme + pat.occ.trend[i,1]/91*extreme.occ.mean[i]/ncol(east.winter)       #### delta fi * pi
      combined.extreme <- combined.extreme + pat.occ.trend[i,1]/91*extreme.occ.trend[i,1]/ncol(east.winter)  #######  delta pi * delta fi
    }
    
    results <- c( thermo, dynamic,combined, present, thermo.heavy, dynamic.heavy, combined.heavy,
                  thermo.extreme,dynamic.extreme,combined.extreme)
  }
  
  results
}


