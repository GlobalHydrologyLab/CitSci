#library
#install.packages("tidyverse")
#install.packages("dplyr")
library("tidyverse")
library("dplyr")

#set working directory
setwd("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/corrs_spearman/all/")

#make a directory
vol_dir = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/vols/"


#list all of the files (REMINDER: CHANGE CORR FILE NAME FOR EACH STATE - lines 29, 32,86)
#for NC
#vol_files = list.files(vol_dir, pattern='N2_Rating_Curve_Outputs')
#for WA
#vol_files = list.files(vol_dir, pattern='W2_Rating_Curve_Outputs')
#for Il
vol_files = list.files(vol_dir, pattern='L2_Rating_Curve_Outputs')

#where to put new files
outpath = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/corrs_spearman/all/"

#create table to store pearsons correlation coefficients
table_data = data.frame()
#create a file to store the table
#file.create("IL_volume_correlation_coefficients")
#make headers
headers = data.frame("Lake 1", "Lake 2", "Spearmans Correlation Coefficient", "P-Value", "n")
write.table(headers, file = "IL_PT_RC_volume_correlation_coefficients", sep=',', append = TRUE, row.names=FALSE, col.names=FALSE) 

#create a file to store the table for sig vals
#file.create("IL_significant_correlations")
write.table(headers, file = "IL_PT_RC_volume_significant_correlations_05", sep=',', append = TRUE, row.names=FALSE, col.names=FALSE) 

write.table(headers, file = "IL_PT_RC_volume_significant_correlations_01", sep=',', append = TRUE, row.names=FALSE, col.names=FALSE) 


ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}



for (jnd in 1:length(vol_files)) {
  print(jnd)
  toRead1 = paste("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/vols/",vol_files[jnd], sep='')
  lk11 = read.csv(toRead1)
  for (ind in 1:length(vol_files)) {
    print(ind)
    #read in hgt file
    if(ind == jnd){
      next 
    }
    
    #print(c(jnd, ind))
    toRead2 = paste("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_rating_curve/vols/",vol_files[ind], sep='')
    lk22 = read.csv(toRead2)
    
    
    hdates1 = as.Date(lk11$Date, format = "%Y-%m-%d", origin = "1904-01-01")
    hdates1_p1 = as.Date(lk11$Date_plus1, format = "%Y-%m-%d", origin = "1904-01-01")
    hdates1_m1 = as.Date(lk11$Date_minus1, format = "%Y-%m-%d", origin = "1904-01-01")
    lk1 = cbind.data.frame(Volume = lk11$X.Volume.RC., Date = hdates1, Date_plus1 = hdates1_p1, Date_minus1 = hdates1_m1)
    hdates2 = as.Date(lk22$Date, format = "%Y-%m-%d", origin = "1904-01-01")
    hdates2_p2 = as.Date(lk22$Date_plus1, format = "%Y-%m-%d", origin = "1904-01-01")
    hdates2_m2 = as.Date(lk22$Date_minus1, format = "%Y-%m-%d", origin = "1904-01-01")
    lk2 = cbind.data.frame(Volume = lk22$X.Volume.RC., Date = hdates2, Date_plus1 = hdates2_p2, Date_minus1 = hdates2_m2)
    
    uh1 = lk1 %>% group_by(Date) %>% summarise(Volume = mean(Volume)) %>% ungroup()
    uh1p1 = lk1 %>% group_by(Date_plus1) %>% summarise(Volume = mean(Volume)) %>% ungroup()
    uh1m1 = lk1 %>% group_by(Date_minus1) %>% summarise(Volume = mean(Volume)) %>% ungroup()
    uh2 = lk2 %>% group_by(Date) %>% summarise(Volume = mean(Volume)) %>% ungroup()
    uh2p1 = lk2 %>% group_by(Date_plus1) %>% summarise(Volume = mean(Volume)) %>% ungroup()
    uh2m1 = lk2 %>% group_by(Date_minus1) %>% summarise(Volume = mean(Volume)) %>% ungroup()
    
    uh11 = cbind.data.frame(Date1 = uh1$Date, Volume1 = uh1$Volume, DatePlus11 = uh1p1$Date_plus1, DateMinus11 = uh1m1$Date_minus1)
    uh22 = cbind.data.frame(Date2 = uh2$Date, Volume2 = uh2$Volume, DatePlus12 = uh2p1$Date_plus1, DateMinus12 = uh2m1$Date_minus1)
    uh22 = cbind.data.frame(uh22, Date2Copy = uh22[, 1])
    
    #do a join date to date
    lk1lk2 =  inner_join(uh11, uh22, by = c("Date1" = "Date2Copy"))
    
    # inner join for lk1 plus one to lk2
    lk1_p1 = inner_join(uh11, uh22, by = c("DatePlus11" = "Date2Copy"))
    
    # inner join for lk1 minus one to lk2
    lk1_m1 = inner_join(uh11, uh22, by = c("DateMinus11" = "Date2Copy"))
    
    # inner join for lk2 plus one to lk1 - DONT NEED
    ##lk2_p1 = inner_join(uh11, uh22, by = c("Date" = "Date_plus1"))
    
    # inner join for lk1 minus one to lk2 - DONT NEED
    ##lk2_m1 = inner_join(uh11, uh22, by = c("Date" = "Date_minus1"))
    
    h_levels = rbind(lk1lk2, lk1_p1, lk1_m1)
    
    if (length(h_levels$Date1) > 9 ){
      #plot
      pattern1 = substr(vol_files[jnd], 1, 4)
      pattern2 = substr(vol_files[ind], 1, 4)
      R = paste(pattern1,'vs', pattern2, "regression (RC)")
      print(R)
     # figure = ggplotRegression(lm(Volume2 ~ Volume1, data = h_levels)) + labs(x = pattern1, y = pattern2)
     # ggsave(figure, filename = paste0(outpath,R,'.jpg'))
      
      #jpeg(paste(outpath,R,'.jpg'))
      #plot(h_levels$Volume1, h_levels$Volume2, pch=7, xlab = '', ylab = expression(paste("Change in Volume (m"^"3", ")")), col='blue') # for posters cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      #abline(lm(Volume2 ~ Volume1, data = h_levels), lty=1, lwd = 4, col = 'blue')  #cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
      #title(main = paste("Change in", R, "Lake Water Volume"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #  dev.off()
      
      
      #get and export pearsons
      h1 = h_levels$Volume1
      h2 = h_levels$Volume2
      corr = cor(h1, h2, method = "spearman")
      corr = round(corr, 2)
      test = cor.test(h1, h2, method = "spearman")
      pval = test$p.value
      pval = round(pval, 3)
      n = length(h_levels$Volume1)
      #create table
      table_data = data.frame(Lake1 = pattern1, Lake2 = pattern2, cor = corr, pvalue = pval, n = n)
      #write the correlation coeffiecint in file             
      write.table(table_data, file = "IL_PT_RC_volume_correlation_coefficients", sep=',', append = TRUE, row.names=FALSE, col.names = FALSE) 
      
      if (!is.na(corr)){
        sig = any(pval <= 0.05)
        if (sig == TRUE){
          table_data2 = table_data
          write.table(table_data, file = "IL_PT_RC_volume_significant_correlations_05", sep=',', append = TRUE, row.names=FALSE, col.names = FALSE)
        } else {
          
        }
      }
      if (!is.na(corr)){
        sig = any(pval <= 0.01)
        if (sig == TRUE){
          table_data2 = table_data
          write.table(table_data, file = "IL_PT_RC_volume_significant_correlations_01", sep=',', append = TRUE, row.names=FALSE, col.names = FALSE)
        } else {
          
        }
      }
      
    }else{
      next()
    }
  }
}
