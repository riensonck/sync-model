#######################################
## Importing Data and Summarizing it ##
#######################################

setwd("~/repos/github_sync-model/Data/1s-Response-Time-Data") # Setting the working directory to the folder containing the data
#setwd("~/repos/github_sync-model/Data/Lapsing-Bounds-Data/") # Setting the working directory to the folder containing the data
files <-list.files(pattern="*.csv")   # collecting all the .csv files in the current working directory
data <- do.call(rbind, lapply(files, function(x){
  df <- read.csv(x)
  str <- gsub("[^0-9]", "",  x) # extracting the numbers from the file name
  # some files have a longer str than others, because the theta frequency values range from 1-20
  # TODO: change the Python code, to save the files names as ex. 01 instead of 1 for theta. 
  if (nchar(str) == 7){ # extracting theta values 1 to 9
    df$Hz <- substr(str, start = 2, stop = 2) # extracting the theta frequency value
    df$thres <- substr(str, start = 5, stop = 5) #  extracting the threshold value 
    df$Gd <- substr(str, start = 6, stop = 6) # extracting the drift rate of the neurons (inverse of the drift of the drift diffusion model)
  } else if (nchar(str == 8)){ # extracting theta values 10 to 20
    df$Hz <- substr(str, start = 2, stop = 3) # extracting the theta frequency value
    df$thres <- substr(str, start = 6, stop = 6) #  extracting the threshold value 
    df$Gd <- substr(str, start = 7, stop = 7) # extracting the drift rate of the neurons (inverse of the drift of the drift diffusion model)
  } else {
    print("file: " + str(x) + " not included")
  }
  return(df)
}))

# need Hz as factor for these plots
data$Hz <- factor(data$Hz, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))

sub_cor <- subset(data, accuracy == 1)
sub_err <- subset(data, accuracy == 0)

rt_err <- data[which(data$accuracy == 0),]
rt_err$rt <- rt_err$rt * -1
rt_cor <- data[which(data$accuracy == 1),]

library(ggplot2)
p2.1 <- ggplot(rt_cor, aes(x = rt, color = Hz)) + 
  facet_grid(cols = vars(thres), rows = vars(Gd), labeller = label_both) +
  geom_density() +
  geom_density(data = rt_err) +
  labs(title = 'RT Distribution ~ frequency(Hz), gamma distance (Gd), threshold', x='reaction time', color = "Hz") 
p2.1
ggsave(paste0('~/repos/github_sync-model/Images/Collapsing-Bounds/rt-distribution_Hz-isi-thres','.jpg'),plot=p2.1, units = 'in', width=14, height =20)
