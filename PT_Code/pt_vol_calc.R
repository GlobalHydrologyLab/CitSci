#library
#install.packages("tidyverse")
#install.packages("dplyr")
library("tidyverse")
library("dplyr")

#set working directory
setwd("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_volume_calcs/")

#make a directory
hgt_dir = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_data/"
area_dir = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/areas/"

#list all of the files
hgt_files = list.files(hgt_dir, pattern='.csv')
area_files = list.files(area_dir, pattern='.csv')

#where to put new files
outpath = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_volume_calcs/"


for(ind in 20:length(area_files)){
  #read in area file
  toRead = paste("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_data/", hgt_files[ind], sep='')
  hgt_vals = read.csv(toRead, header = TRUE, sep = ",", quote =  "")
  
  # these line will find associated height file to current area file
  pattern = hgt_files[ind]
  hgt_ind = which(area_files == pattern)
  area_vals = read.csv(paste(area_dir, area_files[hgt_ind], sep=''))
  
  
  hgt_m = hgt_vals$LEVEL
  
  #dates
  if (hgt_files[hgt_ind] == "BWW2.csv") {
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else if (hgt_files[hgt_ind] == "CAW2.csv"){
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else if (hgt_files[hgt_ind] == "BEW2.csv"){
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else if (hgt_files[hgt_ind] == "LAW2.csv"){
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else if (hgt_files[hgt_ind] == "RSW2.csv"){
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else if (hgt_files[hgt_ind] == "THW2.csv"){
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else if (hgt_files[hgt_ind] == "SAW2.csv"){
    hgt_date = as.Date(substring(hgt_vals$X.Date, 2), "%m/%d/%Y")
  }else{
    hgt_date = as.Date(hgt_vals$Date, "%m/%d/%y")
  }
  
  area_km2 = rev(na.omit(area_vals[,4]))
  area_m2 = na.omit(area_km2 * 1000000)
  area_dates = rev(na.omit(area_vals[,1]))
  
  
  ## find closest height date to area date
  ad = as.vector(area_dates)
  adh = na.omit(as.Date(ad, format = "%b %d,%Y"))

  hd = hgt_date
  hda = na.omit(hd)
  
  #remove outliers in hgts
  both = cbind(hgt_m, hda)
  
  outlier2 = any(hgt_m < 0.05)
  out2 = which (hgt_m < 0.05)
  if (outlier2 == TRUE){
    both = (both[-c(out2), ])
  } else {
    both = both
  }
  
  hgt_m = both[,1]
  hgt_date = both[,2]
  hda = as.Date(hgt_date, "1970-01-01")
  
  #outlier
  bothA = cbind(area_m2, adh)
  
  outlierA = any(area_m2 < 0.1)
  outA = which (area_m2 < 0.1)
  if (outlierA == TRUE){
    bothA = (bothA[-c(outA), ])
  } else {
    bothA = bothA
  }
  
  outlierB = any(area_m2 < median(area_m2)-(median(area_m2)*0.25))
  outB = which (area_m2 < median(area_m2)-(median(area_m2)*0.25))
  if (outlierB == TRUE){
    bothA = (bothA[-c(outB), ])
  } else {
    bothA = bothA
  }
  
  outlierC = any(area_m2 > median(area_m2)+(median(area_m2)*0.25))
  outC = which (area_m2 > median(area_m2)+(median(area_m2)*0.25))
  if (outlierC == TRUE){
    bothA = (bothA[-c(outC), ])
  } else {
    bothA = bothA
  }
  
  
  area_m2 = bothA[,1]
  area_date = bothA[,2]
  adh = as.Date(area_date, "1970-01-01")
  
  
  #PLUS OF MINUS A DAY
  
  # make data frame of hgts
  hgts_J = cbind(hgt_m, hda)
  
  # add date plus one to hgts
  hda_p1 = hda + 1
  hgts_J = cbind(hgt_m, hda, hda_p1)
  
  # add dates minus one to hgts
  hda_m1 = hda - 1
  hgts_J = as.data.frame(cbind(hgt_m, hda, hda_p1, hda_m1))
  
  # make data frame of areas
  area_J = cbind(area_m2, adh)
  
  # add date plus one to hgts
  adh_p1 = adh + 1
  area_J = cbind(area_m2, adh, adh_p1)
  
  # add dates minus one to hgts
  adh_m1 = adh - 1
  area_J = as.data.frame(cbind(area_m2, adh, adh_p1, adh_m1))
  
  #add another area date so that it stays during the join
  adhCopy= adh
  area_J = as.data.frame(cbind(area_m2, adh, adhCopy, adh_p1, adh_m1))
  
  # inner join for date to date
  IJ_ha = inner_join(hgts_J, area_J, by = c("hda" = "adhCopy"))
  
  # inner join for date plus one A to date H
  IJ_hp1a = inner_join(hgts_J, area_J, by = c("hda_p1" = "adhCopy"))
  
  # inner join for date minus one A to date H
  IJ_hm1a = inner_join(hgts_J, area_J, by = c("hda_m1" = "adhCopy"))
  
  # inner join for date plus one H to date A - DONT NEED
  IJ_hap1 = inner_join(area_J, hgts_J, by = c("adh_p1" = "hda"))
  
  # inner join for date minus one H to date A - DONT NEED
  IJ_ham1 = inner_join(area_J, hgts_J, by = c("adh_m1" = "hda"))
  
  # combine all the data together
  mergeAH = rbind(IJ_ha, IJ_hp1a, IJ_hm1a)
  
  if (length(mergeAH$hda) != 0 ){
    
    #average duplicates
    uh = unique(mergeAH$hda)
    dhi = unique(mergeAH$hda[duplicated(mergeAH$hda)])
    dh = dhi[which(!is.na(dhi))]
    hgt_m_repl = vector()
    for(j in 1:length(uh)) {
      if(any(uh[j]==dh)) {
        wi = which((uh[j]==dh) == TRUE)
        d = dh[wi]
        dind= which(mergeAH$hda==d)
        avghj = mean(mergeAH$hgt_m[dind])
        hgt_m_repl = replace(x=mergeAH$hgt_m, list=dind, values=avghj)
        #print(dind) ### to see the indices being averaged each loop around
      }else{
        hgt_m_repl = as.numeric(mergeAH$hgt_m)
      }
    }
    hgt_m_repl = as.numeric(hgt_m_repl)
    
    mergeAH[, "hgt_m"] <- hgt_m_repl
    
    #remove duplicates
    
    w1 = rep(0, length=length(mergeAH$hda))
    if (length(mergeAH$hda) >= 2){
      for(i in 2:length(mergeAH$hda)) {
        if((mergeAH$hda[i]==mergeAH$hda[i-1]) == TRUE) {
          w1[i] = i
        }else{
          w1[i] = 0
        }
      }
    }
    woDup = data.frame()
    for(j in 1:length(uh)) {
      if(any(uh[j]==dh)) {
        woDup = mergeAH[-c(w1), ]
      } else{
        woDup = mergeAH
      }
    }
    
    if(length(woDup) != 0){
      
      
      #pull out measurements for calculations
      woDup = woDup[sort.list(woDup[,2]), ]
      hV = woDup$hgt_m
      aV = woDup$area_m2
      dateInd = woDup$hda
      dateInd = as.Date(dateInd, "1970-01-01")
      
      
      
      #save to list
      save = cbind(dateInd,hV,aV)
      savef=save
      
      ## take difference
      #area
      area_diff_saves = vector()
      for(i in 1:length(aV)){
        area_diff = aV[i] + aV[1]
        #print(area_diff)
        area_diff_saves[i]=area_diff
      }
      
      #hgt
      hgt_diff_saves = as.numeric()
      for(i in 1:length(hV)){
        hgt_diff = hV[i] - hV[1]
        #print(hgt_diff)
        hgt_diff_saves[i]=hgt_diff
      }
      
      ## calculate volume change
      vol_change_saves = vector()
      for(i in 1:length(dateInd)){
        vol_change = hgt_diff_saves[i]*(area_diff_saves[i]/2)
        #print(vol_change)
        vol_change_saves[i]=vol_change
      }
      
      to_save = cbind(save, area_diff_saves, hgt_diff_saves, vol_change_saves)
      
      #create table
      table_data = data.frame(Date=dateInd,Heightm=hV, Aream2=aV,ChangeInHeightm=hgt_diff_saves,ChangeInAream2=area_diff_saves, ChangeInVolumem3=vol_change_saves)
      #create file
     
      #write in file             
      write.table(table_data, file = paste0(outpath,pattern), sep=',', append = FALSE, row.names=FALSE) 
      #delete original file
      
      
      ### plot data
      R = substr(pattern, 1, 4)
      #datehgt = cbind(as.Date(dateInd, "1970-01-01"), vol_change_saves)
      date_plot = dateInd
      jpeg(paste0(outpath,R,'.jpg'))
      plot(date_plot, vol_change_saves, pch=7, xlab = '', ylab = expression(paste("Change in Volume (m"^"3", ")")), col='blue') # for posters cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      lines(date_plot, vol_change_saves, lty=1, lwd = 4, col = 'blue')  #cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
      title(main = paste("Change in", R, "Lake Water Volume"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      dev.off()
      
      print(ind)
    } else {
      
    }
  }else{
    print("no results") 
  }
}

