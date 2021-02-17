#library
#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("chron")
library("tidyverse")
library("dplyr")
library("chron")


#set working directory
setwd("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/")

#where to put new files
outpath = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/"

#make a directory
logger_dir = "/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/"

#list all of the files 
logger_files = list.files(logger_dir, pattern='.csv')

#read them in 
lk1 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/BEW2.csv", header = FALSE)
lk2 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/BTN2.csv", header = FALSE)
lk3 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/BWW2.csv", header = FALSE)
lk4 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/CAW2.csv", header = FALSE)
lk5 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/CFN2.csv", header = FALSE)
lk6 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/DMN2.csv", header = FALSE)
lk7 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/FDN2.csv", header = FALSE)
lk8 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/GRN2.csv", header = FALSE)
lk9 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/HPN2.csv", header = FALSE)
lk10 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/JNN2.csv", header = FALSE)
lk11 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/LAW2.csv", header = FALSE)
lk12 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/MTN2.csv", header = FALSE)
lk13 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/NC1009.csv", header = FALSE)
lk14 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/RSW2.csv", header = FALSE)
lk15 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/SAN2.csv", header = FALSE)
lk16 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/SAW2.csv", header = FALSE)
lk17 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/STN2.csv", header = FALSE)
lk18 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/THW2.csv", header = FALSE)
lk19 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/VCN2.csv", header = FALSE)
lk20 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/WHN2.csv", header = FALSE)
lk21 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/DAL2.csv", header = FALSE)
lk22 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/DFL2.csv", header = FALSE)
lk23 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/GAL2.csv", header = FALSE)
lk24 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/HAL2.csv", header = FALSE)
lk25 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/KLL2.csv", header = FALSE)
lk26 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/LNL2.csv", header = FALSE)
lk27 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/LYL2.csv", header = FALSE)
lk28 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/MCL2.csv", header = FALSE)
lk29 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/MPL2.csv", header = FALSE)
lk30 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/RDL2.csv", header = FALSE)
lk31 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/TML2.csv", header = FALSE)
lk32 = read.csv("/Users/rinab14/Documents/gradschool/lab/pressure_transducer_data/pt_accuracy_outputs/YSL2.csv", header = FALSE)




#measurements
vs = rbind(lk1, lk2)
vs = rbind(vs, lk3)
vs = rbind(vs, lk4)
vs = rbind(vs, lk5)
vs = rbind(vs, lk6)
vs = rbind(vs, lk7)
vs = rbind(vs, lk8)
vs = rbind(vs, lk9)
vs = rbind(vs, lk10)
vs = rbind(vs, lk11)
vs = rbind(vs, lk12)
vs = rbind(vs, lk13)
vs = rbind(vs, lk14)
vs = rbind(vs, lk15)
vs = rbind(vs, lk16)
vs = rbind(vs, lk17)
vs = rbind(vs, lk18)
vs = rbind(vs, lk19)
vs = rbind(vs, lk20)
vs = rbind(vs, lk21)
vs = rbind(vs, lk22)
vs = rbind(vs, lk23)
vs = rbind(vs, lk24)
vs = rbind(vs, lk25)
vs = rbind(vs, lk26)
vs = rbind(vs, lk27)
vs = rbind(vs, lk28)
vs = rbind(vs, lk29)
vs = rbind(vs, lk30)
vs = rbind(vs, lk31)
vs = rbind(vs, lk32)

#length
len = length(vs$V1)

#R squared
r = cor(vs$V1, vs$V4, method = "pearson")
r2 = round(r*r, 3)
# error is actual minus predicted
actual = vs$V1
predicted = vs$V4
error = actual - predicted
# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}
theRMSE = round(rmse(error), 3)
# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}
theMAE = round(mae(error), 3)


# Equation of the line : 
reg<-lm(V4 ~ V1, data = vs)
coeff=coefficients(reg)
eq = paste0("y = ", round(coeff[2],1), "*x + ", round(coeff[1],1))
#plot
lb1 <- paste("r^2 == ", r2)
gh = ggplot(vs, aes(x=V1, y=V4)) + expand_limits(x=0, y=0) + geom_point() + geom_smooth(method="lm") + 
  labs( y="levels by citizen scientists (m)", x="levels from loggers (m)", caption = eq)+ #title="Accuracy of Citizen Scientists", subtitle="Washington, North Carolina, and Illinois Lakes",
  theme_bw()+
  theme(legend.position = "right",legend.title=element_blank(),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.text.x=element_text(size=15, angle = 0, hjust = .5),
        axis.text.y=element_text(size=15, angle = 45, hjust = .5),
        axis.title.x = element_text(size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(size = 17, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        plot.title = element_text(size = 23, hjust=0.5, face="bold", color="black"),
        plot.subtitle = element_text(size = 15, hjust=0.5, color="black"))
g = gh + annotate('text', x = -Inf, vjust = 1.5, hjust = -.1, y = Inf, size=6, label = (paste("     MAE =", theMAE, "m \n   RMSE = ", theRMSE, "m \n         n = ", len) )) +   #+ annotate('text', x = -Inf, vjust = 2.5, hjust = -.4, y = Inf, label = (paste("MAE = ", theMAE) )) + annotate('text', x = -Inf, vjust = 4, hjust = -.325, y = Inf, label = (paste("RMSE = ", theRMSE) ))+
  annotate("text", x = 0.18, y = 1.04, label=lb1, parse=TRUE, color = 'black', size = 6)

plot (g)

#export plot
ggsave(g, filename = paste0(outpath,'All_citsci_accuracy.jpg'))



