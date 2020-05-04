library(ggplot2) # library for plotting
library(ggpubr)
library(dplyr) # library to join dataframes
library(cowplot) # has the function background_grid
library(tidybayes)
library(RColorBrewer) # for color gradient
library(reshape2)
theme_set(theme_tidybayes() + background_grid()) # changing the theme of the plots

################
## Functions  ##
################

##   Source Code: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
##   Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr) # loading plyr locally will mask dplyr globally, and not give rise to conflicting issues
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  # loading the library (plyr) and (dplyr) together will give rise to conflicting functions, and this will result in an Error: All arguments must be named
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## EZ Diffusion Model
## Code adapted from Eric-Jan Wagenmakers, Han L.J. van der Maas , and Raoul P. P. P. Grasman (2007)
## Gives the parameters of the EZ diffusion model
## theta: theta frequencies to be included as a column in the output dataset of this function
## thres: threshold to be included as a column in the output dataset of this function
## drift: neuron drift to be included as a column in the output dataset of this function
## isi: inter-stimulus interval (isi) to be included as a column in the output dataset of this function
## Pc: the proportion of correct decisions
## VRT: the variance of response times for correct decisions
## MRT: mean response time for correct decisions
## s: scaling parameter, the default vaue equals 0.1
## use.isi: boolean, set TRUE if you want to include the isi as a column in the output dataset of this function
get.vaTer = function(theta, drift, thres,  isi = 0, Pc, VRT, MRT, s = 0.1, use.isi = FALSE){
  # The default value for the scaling parameter s equals 0.1
  s2 <- s^2
  # If Pc equals 0, .5, or 1, the method will not work, and
  # an edge correction is required.
  if (Pc == 0){
    cat("Oops, Pc == 0!\n")
  }
  if (Pc == 0.5){
    cat("Oops, Pc == 0.5!\n")
    Pc = 0.5001
  }
  if (Pc == 1){
    cat("Oops, Pc == 1!\n")
  }
  # The function “qlogis” calculates the logit.
  L <- qlogis(Pc)
  x <- L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
  v <- sign(Pc - 0.5) * s * x^(1/4)
  # This gives drift rate.
  a = s2 * qlogis(Pc)/v
  # This gives boundary separation.
  y = -v * a/s2
  MDT = (a/(2 * v)) * (1 - exp(y))/(1 + exp(y))
  Ter = MRT - MDT
  # This gives nondecision time.
  if(use.isi){ # returns the isi as a column in the dataset
    return(c(theta, drift, thres, isi, round(v,4), round(a, 4), round(Ter,4)))
  } else {
    return(c(theta, drift, thres,round(v,4), round(a, 4), round(Ter,4)))
  }
}

# function that results in the sample points where the density plot changes
time_change <- function(data, thetas, drifts, threshs, accuracy){
  df <- vector()
  for(th in thetas){
    for(thr in threshs){
      for(dr in drifts){
        for (acc in accuracy){
          sub <- subset(data, theta == th & drift == dr & thres == thr & accuracy == acc) 
          d <- density(sub$rt)  # reaction time density plot
          d1 <- sign(diff(d$y)/diff(d$x)) # derivative of the density plot
          change <- which(diff(d1) != 0) # the timepoints where the density derivative changes sign
          # using relative timepoints (difference between timepoints) instead of absolute
          T1 <-  d$x[change[2]] - d$x[change[1]]
          if (!is.na(T1)){
            if (T1 > 200){
              T1 <- NaN
              T2 <- NaN
              T3 <- NaN
            } else {
              T2 <-  d$x[change[3]] -  d$x[change[2]] 
              T3 <-  d$x[change[3]] -  d$x[change[1]]
            }
          } else {
            T2 <- NaN
            T3 <- NaN
          }
          V1 <- c(T1, T2) # vector 1
          V2 <- c(T2, T3) # vector 2
          mag <- function(x) sqrt(sum(x^2)) # magnitude function
          rad <- acos(V1 %*% V2 / (mag(V1) * mag(V2))) # calculate the radians of the angle 
          rad2deg <- function(rad) {(rad * 180) / (pi)} # radians to degrees function
          deg <- rad2deg(rad) 
          df <- rbind(df, c(th, thr, dr, acc, "T2 - T1", T1, deg))
          df <- rbind(df, c(th, thr, dr, acc, "T3 - T2", T2, deg))
          df <- rbind(df, c(th, thr, dr, acc, "T3 - T1", T3, deg))
        }
      }
    }
  }
  return(df)
}
time_change <- function(data, thetas, drifts, threshs, accuracy){
  df <- vector()
  for(th in thetas){
    for(thr in threshs){
      for(dr in drifts){
        for (acc in accuracy){
          sub <- subset(data, theta == th & drift == dr & thres == thr & accuracy == acc) 
          d <- density(sub$rt)  # reaction time density plot
          d1 <- sign(diff(d$y)/diff(d$x)) # derivative of the density plot
          change <- which(diff(d1) != 0) # the timepoints where the density derivative changes sign
          # using relative timepoints (difference between timepoints) instead of absolute
          T1 <-  d$x[change[2]] - d$x[change[1]]
          if (!is.na(T1)){
            if (T1 > 200){
              T1 <- NaN
              T2 <- NaN
              T3 <- NaN
            } else {
              T2 <-  d$x[change[3]] -  d$x[change[2]] 
              T3 <-  d$x[change[3]] -  d$x[change[1]]
            }
          } else {
            T2 <- NaN
            T3 <- NaN
          }
          # vector (T1, T2), the x-values are randomly chosen, since they don't matter
          point1 <- c(50., T1) 
          point2 <- c(100., T2) 
          diff <- (point2[2] - point1[2]) / (point2[1] - point1[1])
          rad = atan(diff)
          rad2deg <- function(rad) {(rad * 180) / (pi)} # radians to degrees function
          deg <- rad2deg(rad) 
          df <- rbind(df, c(th, thr, dr, acc, deg))
        }
      }
    }
  }
  return(df)
}

#####################
## Importing Data  ##
#####################

setwd("~/repos/github_sync-model/Data/Collapsing-Bounds-Data/") # Setting the working directory to the folder containing the data
#setwd("~/repos/github_sync-model/Data/Lapsing-Bounds-Data/") # Setting the working directory to the folder containing the data
files <-list.files(pattern="*.csv")   # collecting all the .csv files in the current working directory
data <- do.call(rbind, lapply(files, function(x){
  df <- read.csv(x)
  str <- gsub("[^0-9]", "",  x) # extracting the numbers from the file name
  # some files have a longer str than others, because the theta frequency values range from 1-20
  # TODO: change the Python code, to save the files names as ex. 01 instead of 1 for theta. 
  if (nchar(str) == 7){ # extracting theta values 1 to 9
    df$theta <- substr(str, start = 2, stop = 2) # extracting the theta frequency value
    df$thres <- substr(str, start = 5, stop = 5) #  extracting the threshold value 
    df$drift <- substr(str, start = 6, stop = 6) # extracting the drift rate of the neurons (inverse of the drift of the drift diffusion model)
  } else if (nchar(str == 8)){ # extracting theta values 10 to 20
    df$theta <- substr(str, start = 2, stop = 3) # extracting the theta frequency value
    df$thres <- substr(str, start = 6, stop = 6) #  extracting the threshold value 
    df$drift <- substr(str, start = 7, stop = 7) # extracting the drift rate of the neurons (inverse of the drift of the drift diffusion model)
  } else {
    print("file: " + str(x) + " not included")
  }
  return(df)
  }))

####################################
## Summarizing data - without isi ##
####################################
# used for the plots without isi
# including the isi column in the dataset would result in extra points being plotted for each isi in many plots, which might make some plot interpretations harder. 

# Summarizing the RT data
RT_summary <- summarySE(data, measurevar = "rt", groupvars = c("theta", "drift", "thres"))
names(RT_summary) <- c("theta", "drift", "thres", "N", "RT", "RT_sd", "RT_se", "RT_ci")

# Summarizing the accuracy data
ACC_summary <- summarySE(data, measurevar = "accuracy", groupvars = c("theta", "drift", "thres"))
names(ACC_summary) <- c("theta", "drift", "thres", "N", "ACC", "ACC_sd", "ACC_se", "ACC_ci")

# Summarizing and transforming to drift v, boundary a, and nondecision-time Ter
EZ_data <- matrix(ncol=6, nrow = 112) # nrows needs to be equal to the amount of rows in RT_summary
MRT <- ddply(data, .(theta, drift, thres) , summarize, mean = mean(rt)) # mean response time for correct and uncorrect decisions for each theta
Pc <- ddply(data, .(theta, drift, thres) , summarize, mean = mean(accuracy)) # Pc, proportion of correct decisions
VRT <- ddply(data, .(theta, drift, thres) , summarize, sd = sd(accuracy)) # VRT, variance of response times for correct decisions
for (i in 1:nrow(Pc)){
  vaTer <- get.vaTer(theta = Pc$theta[i], drift = Pc$drift[i], thres = Pc$thres[i], Pc = Pc$mean[i], VRT = VRT$sd[i], MRT = MRT$mean[i])
  EZ_data[i,] <- vaTer
}
colnames(EZ_data) <- c("theta", "drift", "thres", "v", "a", "Ter") # v - Mean drift rate, a - Boundary separation, Ter - Mean of the nondecision component of processing
EZ_data <- data.frame(EZ_data) 
# all columns are of class factor, have to change them to the numeric class. Changing from factor to numeric is meaningless in R,
# and R will just give each factor a numeric postiion value, so the real value is lost. Using a different way: 
EZ_data$theta <- as.numeric(levels(EZ_data$theta))[EZ_data$theta]
EZ_data$v <- as.numeric(levels(EZ_data$v))[EZ_data$v]
EZ_data$a <- as.numeric(levels(EZ_data$a))[EZ_data$a]
EZ_data$Ter <- as.numeric(levels(EZ_data$Ter))[EZ_data$Ter]
# merging together
summary <- left_join(RT_summary, ACC_summary, by = c("theta", "drift", "thres"))
summary$theta <- as.numeric(summary$theta)
summary <- left_join(summary, EZ_data, by = c("theta", "drift", "thres"))


##################################
## Summarizing data - with isi ##
#################################
# used for the plots with isi

# Summarizing the RT data
RT_summary <- summarySE(data, measurevar = "rt", groupvars = c("theta", "drift", "thres", "isi"))
names(RT_summary) <- c("theta", "drift", "thres", "isi", "N", "RT", "RT_sd", "RT_se", "RT_ci")

# Summarizing the accuracy data
ACC_summary <- summarySE(data, measurevar = "accuracy", groupvars = c("theta", "drift", "thres", "isi"))
names(ACC_summary) <- c("theta", "drift", "thres", "isi", "N", "ACC", "ACC_sd", "ACC_se", "ACC_ci")

# Summarizing and transforming to drift v, boundary a, and nondecision-time Ter
EZ_data <- matrix(ncol=7, nrow = 7040)
MRT <- ddply(data, .(theta, drift, thres, isi) , summarize, mean = mean(rt)) # mean response time for correct and uncorrect decisions for each theta
Pc <- ddply(data, .(theta, drift, thres, isi) , summarize, mean = mean(accuracy)) # Pc, proportion of correct decisions
VRT <- ddply(data, .(theta, drift, thres, isi) , summarize, sd = sd(accuracy)) # VRT, variance of response times for correct decisions
for (i in 1:nrow(Pc)){
  vaTer <- get.vaTer(Pc$theta[i], Pc$drift[i], Pc$thres[i], Pc$isi[i], Pc$mean[i], VRT$sd[i], MRT$mean[i], use.isi = TRUE)
  EZ_data[i,] <- vaTer
}
colnames(EZ_data) <- c("theta", "drift", "thres", "isi", "v", "a", "Ter") # v - Mean drift rate, a - Boundary separation, Ter - Mean of the nondecision component of processing
EZ_data <- data.frame(EZ_data) 
# all columns are of class factor, have to change them to the numeric class. Changing from factor to numeric is meaningless in R,
# and R will just give each factor a numeric postiion value, so the real value is lost. Using a different way: 
EZ_data$theta <- as.numeric(levels(EZ_data$theta))[EZ_data$theta]
EZ_data$v <- as.numeric(levels(EZ_data$v))[EZ_data$v]
EZ_data$isi <- as.numeric(levels(EZ_data$isi))[EZ_data$isi]
EZ_data$a <- as.numeric(levels(EZ_data$a))[EZ_data$a]
EZ_data$Ter <- as.numeric(levels(EZ_data$Ter))[EZ_data$Ter]

summary_isi <- left_join(RT_summary, ACC_summary, by = c("theta", "drift", "thres", "isi"))
summary_isi$theta <- as.numeric(summary_isi$theta)
summary_isi <- left_join(summary_isi, EZ_data, by = c("theta", "drift", "thres", "isi"))

###############
## Plotting  ##
###############

############################################
## Parameter Search - Reaction Time (RT)  ##
############################################

# Plot1: Mean RT ~ theta, drift, threshold
p1 <- ggplot(summary, aes(x=theta, y = RT, color=drift, group = drift)) + 
  facet_wrap( ~ thres, labeller =  label_both) +
  geom_point() + 
  labs(title = 'Mean RT ~ theta, drift, threshold', y='reaction time') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=RT-RT_ci, ymax=RT+RT_ci), width=.1)
p1

sub1 <- subset(summary_isi, theta == 3 & drift == 1 & thres == 3)
p2 <- ggplot(sub1, aes(x = isi, y = v, group = factor(theta), color= factor(theta))) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs(title = '', x='isi', y ="drift v", color = "Hz") +
  theme_linedraw()
p2 

n_samples <-length(sub1$v) 
fft <- fft(sub1$v)
amp <- abs(fft)
mask_good_freq <- seq(1, n_samples/2 + 1)
amp[mask_good_freq]

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Plot 2
thetas <- sort(unique(as.numeric(data$theta)))
drifts <- sort(unique(as.numeric(data$drift)))
threshs <- sort(unique(as.numeric(data$thres)))
accuracy <- sort(unique(as.numeric(data$accuracy)))
summary_change <- time_change(data, thetas, drifts, threshs, accuracy)
summary_change <- as.data.frame(summary_change)

colnames(summary_change) <- c("theta", "thres", "drift", "accuracy", "time", "rel_timepoint", "degrees")
summary_change$theta <- as.numeric(levels(summary_change$theta))[summary_change$theta]
summary_change$rel_timepoint <- as.numeric(levels(summary_change$rel_timepoint))[summary_change$rel_timepoint]
summary_change$degrees <- as.numeric(levels(summary_change$degrees))[summary_change$degrees]
summary_change$time <- factor(summary_change$time, levels = c("T2 - T1","T3 - T2", "T3 - T1"))

sub2 <- subset(data, theta == 3 & accuracy == 1)
d <- density(sub2$rt)
plot(d)

sub2.1 <- subset(summary_change, theta < 4 & thres == 3)
sub2.2 <- subset(summary_change, theta > 3 & theta < 9 & thres == 3)
sub2.3 <- subset(summary_change, theta > 8 & theta < 13 & thres == 3)
sub2.4 <- subset(summary_change, theta > 12 & theta < 15 & thres == 3)

sub2 <- subset(summary_change, (theta == 15 | theta == 8 | theta == 2 | theta == 3) & thres == 3 )
sub2 <- subset(summary_change, thres == 3)
p2 <- ggplot(sub2, aes(x = time, y = rel_timepoint, group = factor(theta), color= factor(theta))) + 
  facet_grid(cols = vars(drift), rows = vars(accuracy), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs(title = '', x='time intervals', y ="time difference", color = "Hz") +
  theme_linedraw()
p2 

aggregate(summary_change$degrees, by=list(summary_change$theta, summary_change$accuracy), FUN = mean, na.rm=TRUE)
a <- aggregate(summary_change$degrees, by=list(summary_change$theta, summary_change$accuracy, summary_change$drift), FUN = mean, na.rm=TRUE)
colnames(a) <- c("theta", "accuracy", "drift", "degrees")
p2 <- ggplot(a, aes(x = factor(theta), y = degrees, group = factor(drift), color= factor(theta))) + 
  facet_grid(cols = vars(drift), rows = vars(accuracy), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs( x='Hz', y ="angle degrees", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p2 

a <- aggregate(summary_change$degrees, by=list(summary_change$theta, summary_change$accuracy, summary_change$drift, summary_change$thres), FUN = mean, na.rm=TRUE)
colnames(a) <- c("theta", "accuracy", "drift", "thres", "degrees")
p2 <- ggplot(a, aes(x = factor(theta), y = degrees, group = accuracy, color= accuracy)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs( x='Hz', y ="angle degrees", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=15)) +
p2 

colnames(summary_change) <- c("theta", "thres", "drift", "accuracy", "degrees")
summary_change$theta <- as.numeric(levels(summary_change$theta))[summary_change$theta]
summary_change$degrees <- as.numeric(levels(summary_change$degrees))[summary_change$degrees]
sub2 <- subset(summary_change, thres == 3)
p2 <- ggplot(sub2, aes(x = factor(theta), y = degrees, group = factor(drift), color= factor(theta))) + 
  facet_grid(cols = vars(drift), rows = vars(accuracy), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs( x='Hz', y ="angle degrees", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = -60, size = 4, color = "black") +
  theme(text = element_text(size=15)) +
  theme_linedraw() +
  geom_hline(aes(yintercept = 0)) 
p2 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Plot 2: RT Distribution ~ theta, drift, threshold
data$theta <- factor(data$theta, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
sub2 <- subset(data, theta == 5 | theta == 9)
p2 <- ggplot(sub2, aes(x = rt, color= theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_density() +
  labs(title = 'RT Distribution ~ theta, drift, threshold', x='reaction time', color = "theta") 
p2 
ggsave(paste0('~/repos/github_sync-model/Images/1s-Response-Time/p2','.jpg'),plot=p2, units = 'in', width=14, height =20)
# ggsave(paste0('~/repos/github_sync-model/Images/Lapsing-Bounds/p2','.jpg'),plot=p2, units = 'in', width=14, height =20)

# Plot 2.1: RT Distribution ~ theta, drift, threshold
data$theta <- factor(data$theta, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
s2 <- subset(data, thres == 3)
p2 <- ggplot(data, aes(x = rt, color= theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_density() +
  labs(title = 'RT Distribution ~ theta, drift, threshold 3', x='reaction time', color = "theta") 
p2 
ggsave(paste0('~/repos/github_sync-model/Images/1s-Response-Time/p2','.jpg'),plot=p2, units = 'in', width=14, height =20)

s2 <- subset(data, drift == 1)
p2 <- ggplot(data, aes(x = rt, color= theta)) + 
  facet_grid(cols = vars(thres), rows = vars(isi), labeller = label_both) +
  geom_density() +
  labs(title = 'RT Distribution ~ theta, drift 1, threshold', x='reaction time', color = "theta") 
p2 
ggsave(paste0('~/repos/github_sync-model/Images/1s-Response-Time/p2-drift1','.jpg'),plot=p2, units = 'in', width=14, height =20)

## Plot 3: RT Correct vs Error Distribution ~ Theta
data$distribution <- as.factor(data$accuracy)
data$distribution <- mapvalues(data$distribution, from = c(0, 1), to = c("Error", "Correct"))
mu <- ddply(data, c("distribution", "theta"), summarise, grp.mean = mean(rt))
p3 <- ggplot(data, aes(x = rt, color = distribution) ) + geom_density() +
  geom_vline(data = mu, aes(xintercept=grp.mean, color=distribution),linetype="dashed", size = 1) + # Add mean lines
  facet_wrap(~ theta, labeller = label_both) +
  labs(title = ' RT Correct vs Error Distribution ~ theta', x='reaction time') 
p3

## Plot 4: RT Correct vs Error Distribution ~ drift, threshold
data$distribution <- as.factor(data$accuracy)
data$distribution <- mapvalues(data$distribution, from = c(0, 1), to = c("Error", "Correct"))
mu1 <- ddply(data, c("distribution", "drift", "thres"), summarise, grp.mean = mean(rt))
p4 <- ggplot(data, aes(x = rt, color = distribution) ) + geom_density() +
  geom_vline(data = mu1, aes(xintercept=grp.mean, color=distribution),linetype="dashed", size = 1) + # Add mean lines
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  labs(title = ' RT Correct vs Error Distribution ~ drift, threshold', x='reaction time') 
p4

## Plot 5: RT Correct vs Error Distribution ~ theta and different responses
data$response <- mapvalues(data$response, from = c(0, 1, 2, 3), to = c("LH-LF", "LH-RF", "RH-LF", "RH-RF"))
# left hand: left finger vs right finger, comparing the mean RTs of correct and error responses, no bias, Check!
data$theta <- as.numeric(data$theta)
data_sub <- subset(data, response != -1 & theta < 11)
data_sub$theta <- factor(data_sub$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
mu2 <- ddply(data_sub, .(response,distribution, theta), summarize, mean = mean(rt))
p5 <- ggplot(data_sub, aes(x = rt, color = distribution) ) + geom_density() +
  geom_vline(data = mu3, aes(xintercept=mean, color=distribution),linetype="dashed") + 
  facet_grid(cols = vars(response), rows = vars(theta), labeller = label_both) 
p5

########################################
## Parameter Search - Accuracy (ACC)  ##
########################################

# Plot 6: Mean Accuracy ~ theta, drift, threshold
p6 <- ggplot(summary, aes(x=theta, y = ACC, color=drift, group = drift)) + 
  facet_wrap( ~ thres, labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right     geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1) + 
  labs(title = 'Mean Accuracy ~ theta, drift, threshold', y ="accuracy") 
p6

# Plot 7: Accuracy Distribution ~ drift, threshold
myPalette <- colorRampPalette(brewer.pal(9, "PuBu"))
sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0.6, 1))
p7 <- ggplot(data, aes(x = accuracy, fill = ..prop..)) + 
  sc +
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) + 
  geom_bar(aes(y = ..prop..), colour="black") +
  labs(title = 'Accuracy Distribution ~ drift, threshold', x='accuracy') +
  scale_x_continuous(breaks=c(0,1)) +
  ylim(0,1) +
  #geom_line(aes(y= ..prop..),stat= "count", position = position_dodge(0.1), color = "magenta4", size = 1) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.1)
p7 

# Plot 8: Accuracy Distribution ~ theta
sub8 <- data
sub8$theta <- factor(sub8$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
myPalette <- colorRampPalette(brewer.pal(9, "PuBu"))
sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0.5, 1))
p8 <- ggplot(sub8, aes(x = accuracy, fill = ..prop..)) + 
  sc +
  facet_wrap(~ theta, labeller = label_both) +
  geom_bar(aes(y = ..prop..), colour="black") +
  labs(title = 'Accuracy Distribution ~ theta', x='accuracy') +
  scale_x_continuous(breaks=c(0,1)) +
  ylim(0,1) +
  #geom_line(aes(y= ..prop..),stat= "count", position = position_dodge(0.1), color = "magenta4", size = 1) +
  geom_text(aes( label = scales::percent(round(..prop..,2)),
                 y= ..prop.. ), stat= "count", vjust = -.1)
p8 

# Plot 9: accuracy - ISI (Delta frequencies)
isi <- subset(summary_isi, theta < 4)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
p9 <- ggplot(isi, aes(x=isi, y = ACC, color=theta, group = theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Delta)', y='accuracy', color ="delta") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p9

# Plot 10: accuracy - ISI (Theta frequencies)
isi <- subset(summary_isi, theta < 8 & theta > 3)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
p10 <- ggplot(isi, aes(x=isi, y = ACC, color=theta, group = theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Theta)', y='accuracy', color ="theta") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p10

# Plot 11: accuracy - ISI (Alpha frequencies)
isi <- subset(summary_isi, theta < 13 & theta > 7)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
p11 <- ggplot(isi, aes(x=isi, y = ACC, color=theta, group = theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Alpha)', y='accuracy', color ="alpha") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p11

###############################
## Parameter Search - drift  ##
###############################
# Plot 12: Drift v ~ theta, drift, threshold
p12 <- ggplot(summary, aes(x=theta, y = v, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Drift v ~ theta, drift, threshold', y='drift v') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p12

# Plot 13: Boundary a ~ theta, drift, threshold
p13 <- ggplot(summary, aes(x=theta, y = a, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Boundary a ~ theta, drift, threshold', y='boundary a') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.215, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p13

# Plot 14: Nondecision time Ter ~ theta, drift, threshold3
sub14 <- subset(summary, thres == 3)
p14 <- ggplot(sub14, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 230, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p14

# Plot 15: Nondecision time Ter ~ theta, drift, threshold4
sub15 <- subset(summary, thres == 4)
p15 <- ggplot(sub15, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 300, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p15

# Plot 16: Nondecision time Ter ~ theta, drift, threshold5
sub16 <- subset(summary, thres == 5)
p16 <- ggplot(sub16, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y =550, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p16

# Plot 17: Nondecision time Ter ~ theta, drift, threshold6
sub17 <- subset(summary, thres == 6)
p17 <- ggplot(sub17, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 650, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p17

##################################
## Parameter Correlation Plots  ##
##################################

# Plot 18: correlation Rt - accuracy ~ drift, thres
p18 <- ggplot(summary, aes(x=RT, y = ACC, color=factor(theta), group = thres)) + 
  geom_point(shape = 1,colour = "black", size = 2, position = "jitter") +
  geom_point(size = 1, position = "jitter") + 
  labs(title = 'Correlation RT - Accuracy', y='accuracy', x = "RT", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 200, label.y = 0.4, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p18

#Plot 19: correlation RT-accuracy ~ theta 15
sub19 <- subset(summary, theta == 15)
p19 <- ggplot(sub19, aes(x=RT, y = ACC, color=factor(thres))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation RT - Accuracy ~ theta 18', y='accuracy', x = "RT", color = "threshold") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 300, label.y = 0.8, size = 4, color = "black") +
  theme(text = element_text(size=15)) +
  facet_wrap(~ drift, labeller = label_both )
 p19
 
#Plot 20: correlation RT-accuracy ~ theta 6
sub20 <- subset(summary, theta == 6)
p20 <- ggplot(sub20, aes(x=RT, y = ACC, color=factor(thres))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation RT - Accuracy ~ theta 6', y='accuracy', x = "RT", color = "threshold") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 300, label.y = 0.8, size = 4, color = "black") +
  theme(text = element_text(size=15)) +
  facet_wrap(~ drift, labeller = label_both )
p20

# Plot 21: Correlation drift v - boundary a
p21 <- ggplot(summary, aes(x=v, y = a, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - boundary a', y='boundary a', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.05, label.y = 0.21, size = 5, color = "black") +
  theme(text = element_text(size=15)) 
p21

# Plot 22: Correlation drift v - Ter
p22 <- ggplot(summary, aes(x=v, y = Ter, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Ter', y='Ter', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.15, label.y = 500, size = 5, color = "black") +
  theme(text = element_text(size=15)) 
p22

# Plot 23: Correlation drift v - boundary a ~ threshold
p23 <- ggplot(summary, aes(x=v, y = Ter, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Ter ~ threshold', y='Ter', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.05, label.y = 300, size = 5, color = "black") +
  theme(text = element_text(size=15))  + 
  facet_wrap(~ thres, labeller = label_both)
p23

# Plot 24: Correlation drift v - accuracy
p24 <- ggplot(summary, aes(x=v, y = ACC, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Accuracy', y='accuracy', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.05, label.y = 1, size = 5, color = "black") +
  theme(text = element_text(size=15))  
p24