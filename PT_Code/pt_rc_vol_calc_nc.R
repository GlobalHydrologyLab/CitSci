##rating curve

library("tidyverse")
library("dplyr")
library("ggplot2")
library ("gridExtra")

#set working directory
setwd("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/vols/")

#where to put new files
outpath = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/"
outpath2 = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/vols/"

#make a vol directory
vol_dir = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_volume_calcs/"

#list all of the vol files
#NC
vol_files = list.files(vol_dir, pattern='N2.csv')

#make a level directory
hgt_dir = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_data/logger_data_clean/"

#list all of the level files
#NC
hgt_files = list.files(hgt_dir, pattern='N2.csv')

#par(mfrow = c(4, 3))

#doesnt always run through on its own and needs some prompting for no reason

for (jnd in 1:length(vol_files)) {
  print(jnd)
  toRead1 = paste("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_volume_calcs/",vol_files[jnd], sep='')
  lk11 = read.csv(toRead1)
  
  forcurve = cbind.data.frame(volume = lk11$ChangeInVolumem3,height = lk11$ChangeInHeightm)
  
  #plot
  pattern1 = substr(vol_files[jnd], 1, 4)
  R = paste(pattern1)
  g = ggplot(forcurve, aes(x=height, y=volume)) + geom_point() #+ geom_smooth(method="lm")
  #plot (g)
  # Equation of the line : 
  reg<-lm(volume ~ height, data = forcurve)
  coeff=coefficients(reg)
  eq = paste0("y = ", round(coeff[2],1), "*x + ", round(coeff[1],1))
  figure = g + labs(title="Rating Curve", subtitle=R, y="Change in Volume", x="Change in Height", caption = eq) +
    theme(axis.title.x = element_text(size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          plot.title = element_text(size = 20))
  #print(figure)
  #print(coef(lm(volume ~ height, data = forcurve)))
  fig = figure + stat_smooth(method="lm", se=FALSE)
  print(fig)
  
  #export plot
  ggsave(fig, filename = paste0(outpath2,R,'.jpg'))
  
  #############extrapolate################################
  
  #read in heights
  
  toRead2 = paste("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/logger_data/logger_data_clean/",hgt_files[jnd], sep='')
  lk22 = read.csv(toRead2)
  
  hV = na.omit(lk22$LEVEL)
  hd = as.vector(lk22$Date)
  hdaa = as.Date(lk22$Date, format = "%m/%d/%y") 
  hdaaa = na.omit(hdaa)
  hhda = as.Date(hdaaa, "1970-01-01")
  
  #remove outliers in hgts
  both = cbind(hV, hhda)
  
  
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
    mergeAH = cbind.data.frame(hgt_m, hda)
  
  #average height duplicates
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
  
  #remove the hgt duplicates
  
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
  if(any(w1>0)){
    woDup = mergeAH[-c(w1), ]
  }else{
    woDup = mergeAH
  }
  
  hgt_m2 = woDup[,1] 
  hda2 = woDup[,2]
  
  hgt_diff_saves = as.numeric()
  for(i in 1:length(hgt_m2)){
    hgt_diff = hgt_m2[i] - hgt_m2[1]
    #print(hgt_diff)
    hgt_diff_saves[i]=hgt_diff
  }
  
  hdates_p2 = hda2 + 1
  hdates_m2 = hda2 - 1
  lk2 = cbind.data.frame(Height = hgt_diff_saves, Date1 = hda2, Date_copy = hda2, Date_plus1 = hdates_p2, Date_minus1 = hdates_m2)
  
  
  #pull out of the equation
  #plug in all heights to equation
  
  volRC = list()
  slope = round(coeff[2],1)
  int = round(coeff[1],1)
  for(i in 1:length(lk2$Height)){
    x = slope*lk2$Height
    y = x + int
    volRC = y
  }
  
  
  #put those volumes back with the dates
  lkNew = cbind.data.frame(Height = hgt_diff_saves, VolumeRC = volRC, Date1 = hda2, Date_plus1 = hdates_p2, Date_minus1 = hdates_m2)
  
  #export chart
  headers = data.frame("^Height", "^Volume(RC)", "Date", "Date_plus1", "Date_minus1")
  write.table(headers, file = paste0(R,'_Rating_Curve_Outputs'), sep=',', append = TRUE, row.names=FALSE, col.names=FALSE) 
  write.table(lkNew, file = paste0(R,'_Rating_Curve_Outputs'), sep=',', append = TRUE, row.names=FALSE, col.names = FALSE) 
  
  #plot
  pattern1 = substr(hgt_files[jnd], 1, 4)
  R = paste(pattern1)
  gh = ggplot(lkNew, aes(x=Date1, y=Height)) + geom_point() + labs(title="Change in Height", subtitle=R, y="Height (m)", x="Dates")+
    theme(axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          plot.title = element_text(size = 20))
  #plot (gh)
  gv = ggplot(lkNew, aes(x=Date1, y=VolumeRC)) + geom_point() + labs(title="Change in Volume", subtitle=R, y="Volume (m3)", x="Dates", caption = eq)+
    theme(axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          plot.title = element_text(size = 20))
  #plot (gv)
  rchv = grid.arrange(gh, gv, ncol=2) 
  print(rchv)
  print(coef(lm(volume ~ height, data = forcurve)))
  #gh + stat_smooth(method="lm", se=FALSE)
  
  #export plot
  ggsave(rchv, filename = paste0(outpath2,R,'_RCHgtVol.jpg'))
  
  
}     

