  library(ggplot2) # library for plotting
  library(ggpubr)
  library(dplyr) # library to join dataframes
  library(cowplot) # has the function background_grid
  library(tidybayes)
  library(RColorBrewer) # for color gradient
  library(reshape2)
  library(ggstance)
  require(gridExtra)
  library(grid)
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
  get.vaTer = function(Hz, Gd, thres,  isi = 0, Pc, VRT, MRT, s = 0.1, use.isi = FALSE){
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
      return(c(Hz, Gd, thres, isi, round(v,4), round(a, 4), round(Ter,4)))
    } else {
      return(c(Hz, Gd, thres,round(v,4), round(a, 4), round(Ter,4)))
    }
  }
  
  # function that results in the sample points where the density plot changes
  time_change <- function(data, frequencies, Gds, threshs){
    df <- vector()
    for(freq in frequencies){
      for(thr in threshs){
        for(gd in Gds){
            sub <- subset(data, Hz == freq & Gd == gd & thres == thr) 
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
            df <- rbind(df, c(freq, thr, gd, "T2 - T1", T1, deg))
            df <- rbind(df, c(freq, thr, gd, "T3 - T2", T2, deg))
            df <- rbind(df, c(freq, thr, gd, "T3 - T1", T3, deg))
          
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
  
  # This function summarizes all the data into two dataframes, one without and one with isi (separating them makes it eassier to plot)
  summarize_data <- function(data, use.isi = FALSE){
    if(use.isi){
      # Summarizing the RT data
      RT_summary.isi <- summarySE(data, measurevar = "rt", groupvars = c("Hz", "Gd", "thres", "isi"))
      names(RT_summary.isi) <- c("Hz", "Gd", "threshold", "isi", "N", "RT", "RT_sd", "RT_se", "RT_ci")
      # Summarizing the accuracy data
      ACC_summary.isi <- summarySE(data, measurevar = "accuracy", groupvars = c("Hz", "Gd", "thres", "isi"))
      names(ACC_summary.isi) <- c("Hz", "Gd", "threshold", "isi", "N", "ACC", "ACC_sd", "ACC_se", "ACC_ci")
      
      # Summarizing and transforming to drift v, boundary a, and nondecision-time Ter
      EZ_data.isi <- matrix(ncol=7, nrow = nrow(RT_summary.isi)) # nrows needs to be equal to the amount of rows in RT_summary
      MRT.isi <- ddply(data, .(Hz, Gd, thres, isi) , summarize, mean = mean(rt)) # mean response time for correct and uncorrect decisions for each theta
      Pc.isi <- ddply(data, .(Hz, Gd, thres, isi) , summarize, mean = mean(accuracy)) # Pc, proportion of correct decisions
      VRT.isi <- ddply(data, .(Hz, Gd, thres, isi) , summarize, sd = sd(accuracy)) # VRT, variance of response times for correct decision
      for (i in 1:nrow(Pc.isi)){
        vaTer <- get.vaTer(Pc.isi$Hz[i], Pc.isi$Gd[i], Pc.isi$thres[i], Pc.isi$isi[i], Pc.isi$mean[i], VRT.isi$sd[i], MRT.isi$mean[i], use.isi = TRUE)
        EZ_data.isi[i,] <- vaTer
      }
      colnames(EZ_data.isi) <- c("Hz", "Gd", "threshold", "isi", "v", "a", "Ter") # v - Mean drift rate, a - Boundary separation, Ter - Mean of the nondecision component of processing
      EZ_data.isi <- data.frame(EZ_data.isi) 
      
      # all columns are of class factor, have to change them to the numeric class. Changing from factor to numeric is meaningless in R,
      # and R will just give each factor a numeric postiion value, so the real value is lost. Using a different way: 
      EZ_data.isi$Hz <- as.numeric(levels(EZ_data.isi$Hz))[EZ_data.isi$Hz]
      EZ_data.isi$v <- as.numeric(levels(EZ_data.isi$v))[EZ_data.isi$v]
      EZ_data.isi$isi <- as.numeric(levels(EZ_data.isi$isi))[EZ_data.isi$isi]
      EZ_data.isi$a <- as.numeric(levels(EZ_data.isi$a))[EZ_data.isi$a]
      EZ_data.isi$Ter <- as.numeric(levels(EZ_data.isi$Ter))[EZ_data.isi$Ter]
      
      # merging together
      summary.isi <- left_join(RT_summary.isi, ACC_summary.isi, by = c("Hz", "Gd", "threshold", "isi"))
      summary.isi$Hz <- as.numeric(summary.isi$Hz)
      summary.isi <- left_join(summary.isi, EZ_data.isi, by = c("Hz", "Gd", "threshold", "isi"))
      return(summary.isi)
      
    } else {
      # Summarizing the RT datas
      RT_summary <- summarySE(data, measurevar = "rt", groupvars = c("Hz", "Gd", "thres"))
      names(RT_summary) <- c("Hz", "Gd", "threshold", "N", "RT", "RT_sd", "RT_se", "RT_ci")
      # Summarizing the accuracy data
      ACC_summary <- summarySE(data, measurevar = "accuracy", groupvars = c("Hz", "Gd", "thres"))
      names(ACC_summary) <- c("Hz", "Gd", "threshold", "N", "ACC", "ACC_sd", "ACC_se", "ACC_ci")
      
      # Summarizing and transforming to drift v, boundary a, and nondecision-time Ter
      EZ_data <- matrix(ncol=6, nrow = nrow(RT_summary))
      MRT <- ddply(data, .(Hz, Gd, thres) , summarize, mean = mean(rt)) # mean response time for correct and uncorrect decisions for each theta
      Pc <- ddply(data, .(Hz, Gd, thres) , summarize, mean = mean(accuracy)) # Pc, proportion of correct decisions
      VRT <- ddply(data, .(Hz, Gd, thres) , summarize, sd = sd(accuracy)) # VRT, variance of response times for correct decisions
      for (i in 1:nrow(Pc)){
        vaTer <- get.vaTer(Hz = Pc$Hz[i], Gd = Pc$Gd[i], thres = Pc$thres[i], Pc = Pc$mean[i], VRT = VRT$sd[i], MRT = MRT$mean[i])
        EZ_data[i,] <- vaTer
      }
      colnames(EZ_data) <- c("Hz", "Gd", "threshold", "v", "a", "Ter") # v - Mean drift rate, a - Boundary separation, Ter - Mean of the nondecision component of processing
      EZ_data <- data.frame(EZ_data)
      
      # all columns are of class factor, have to change them to the numeric class. Changing from factor to numeric is meaningless in R,
      # and R will just give each factor a numeric postiion value, so the real value is lost. Using a different way: 
      EZ_data$Hz <- as.numeric(levels(EZ_data$Hz))[EZ_data$Hz]
      EZ_data$v <- as.numeric(levels(EZ_data$v))[EZ_data$v]
      EZ_data$a <- as.numeric(levels(EZ_data$a))[EZ_data$a]
      EZ_data$Ter <- as.numeric(levels(EZ_data$Ter))[EZ_data$Ter]
      
      # merging together
      summary <- left_join(RT_summary, ACC_summary, by = c("Hz", "Gd", "threshold"))
      summary$Hz <- as.numeric(summary$Hz)
      summary <- left_join(summary, EZ_data, by = c("Hz", "Gd", "threshold"))
      return(summary)
    }
  }
  
  #######################################
  ## Importing Data and Summarizing it ##
  #######################################
  
  setwd("~/repos/github_sync-model/Data/model_data/collapsing_bounds_data/") # Setting the working directory to the folder containing the data
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
      #df$Hz <- substr(str, start = 2, stop = 3) # extracting the theta frequency value
      df$Hz <- substr(str, start = 2, stop = 2) # extracting the gd value
      # df$thres <- substr(str, start = 6, stop = 6) #  extracting the threshold value 
      df$thres <- substr(str, start = 5, stop = 5) #  extracting the gd value 
      #df$Gd <- substr(str, start =7 , stop = 7) # extracting the drift rate of the neurons (inverse of the drift of the drift diffusion model)
      df$Gd <- substr(str, start  = 6 , stop = 7) # extracting the drift rate of the neurons (inverse of the drift of the drift diffusion model)
    } else {
      print("file: " + str(x) + " not included")
    }
    return(df)
    }))
  # renaming columns to less confusing names. 
  # Frequency range is from 2-15Hz, so this isn't only the theta range. 
  # renaming from "theta" to "Hz" 
  names(data)[11] <- c("Hz")
  # neuronal drift in the model,renaming "drift" to "Gd" (Gamma distance)
  names(data)[13] <- c("Gd")
  
  # Summarizing the data
  summary <- summarize_data(data)
  summary.isi <- summarize_data(data, use.isi = TRUE)

#########################################################
## Plotting - figures for parameter space exploration  ##
#########################################################

############################################
## Parameter Search - Reaction Time (RT)  ##
############################################

##############
## Mean RT  ##
##############

# Plot 1.1: Mean RT ~ Hz, Gd, threshold
sub1.1 <- subset(summary, threshold == 4)
p1.1 <- ggplot(sub1.1, aes(x=Hz, y = RT, color=Gd, group = Gd)) + 
  facet_wrap( ~ threshold, labeller =  label_both) +
  geom_point() + 
  labs(title = 'Mean RT', y='reaction time', x = 'Hz', color = 'Gd') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=RT-RT_ci, ymax=RT+RT_ci), width=.1) +
  theme(text = element_text(size=15)) 
p1.1
ggsave(paste0('~/repos/github_sync-model/Results/Collapsing-Bounds/rt-mean-gamma-distance-paper','.jpg'),plot=p1.1, units = 'in', width=10, height =8)


# Plot 1.2: Mean RT ~ Hz, Gd, threshold = 4 (single example plot)
sub1.2 <- subset(summary, threshold == 4,)
p1.2 <- ggplot(sub1.2, aes(x=Hz, y = RT, color=Gd, group = Gd)) + 
  facet_wrap( ~ threshold, labeller =  label_both) +
  geom_point() + 
  labs(title = 'Mean RT ~ frequency (Hz), gamma distance (Gd), threshold', y='reaction time', x = 'Hz', color = 'Gd') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=RT-RT_ci, ymax=RT+RT_ci), width=.1)
p1.2

######################
## RT Distribution  ##
######################

# need Hz as factor for these plots
data$Hz <- factor(data$Hz, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))

# Plot 2.1: RT Distribution ~ Hz, isi, threshold
p2.1 <- ggplot(data, aes(x = rt, color= as.factor(Hz))) + 
  facet_grid(cols = vars(thres), rows = vars(isi), labeller = label_both) +
  geom_density() +
  labs(title = 'RT Distribution ~ frequency(Hz), gamma distance (Gd), threshold', x='reaction time', color = "Hz") 
p2.1
ggsave(paste0('~/repos/github_sync-model/Images/Collapsing-Bounds/rt-distribution_Hz-isi-thres','.jpg'),plot=p2.1, units = 'in', width=14, height =20)

# Plot 2.2: RT Distribution ~ Hz, Gd, threshold, (example plot using only couple frequencies) 
sub2.2 <- subset(data, (Hz == 4 | Hz == 8 | Hz == 15) & thres != 1)
p2.2 <- ggplot(data, aes(x = rt, color= Hz)) + 
  facet_grid(cols = vars(thres), rows = vars(Gd), labeller = label_both) +
  geom_density() +
  labs(title = 'RT Distribution ~ frequency (Hz), gamma distance (Gd), threshold', x='reaction time', color = "Hz") 
p2.2
ggsave(paste0('~/repos/github_sync-model/Images/1s-Response-Time/p2','.jpg'),plot=p2, units = 'in', width=14, height =20)
# ggsave(paste0('~/repos/github_sync-model/Images/Lapsing-Bounds/p2','.jpg'),plot=p2, units = 'in', width=14, height =20)

######################
## RT Peak Timing   ##
######################
frequencies <- sort(unique(as.numeric(levels(data$Hz))[data$Hz]))
Gds <- sort(unique(as.numeric(data$Gd)))
threshs <- sort(unique(as.numeric(data$thres)))
#accuracy <- sort(unique(as.numeric(data$accuracy)))
summary_change <- time_change(data, frequencies, Gds, threshs)
summary_change <- as.data.frame(summary_change)

colnames(summary_change) <- c("Hz", "thres", "Gd", "time", "rel_timepoint", "degrees")
summary_change$Hz <- as.numeric(levels(summary_change$Hz))[summary_change$Hz]
summary_change$rel_timepoint <- as.numeric(levels(summary_change$rel_timepoint))[summary_change$rel_timepoint]
summary_change$degrees <- as.numeric(levels(summary_change$degrees))[summary_change$degrees]
summary_change$time <- factor(summary_change$time, levels = c("T2 - T1","T3 - T2", "T3 - T1"))

# Plot 3.1: 
sub2.1 <- subset(summary_change, theta < 4 & thres == 3)
sub2.2 <- subset(summary_change, theta > 3 & theta < 9 & thres == 3)
sub2.3 <- subset(summary_change, theta > 8 & theta < 13 & thres == 3)
sub2.4 <- subset(summary_change, theta > 12 & theta < 15 & thres == 3)

sub2 <- subset(summary_change, (Hz == 15 | Hz == 8 | Hz == 4) & thres != 1 & thres != 2)
sub2$Hz <- factor(sub2$Hz, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
p2 <- ggplot(sub2, aes(x = time, y = rel_timepoint, group = factor(Hz), color= factor(Hz))) + 
  facet_grid(cols = vars(Gd), rows = vars(thres), labeller = label_both) +
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
  theme(text = element_text(size=15)) 
p2 

colnames(summary_change) <- c("theta", "thres", "drift", "accuracy", "degrees")
#summary_change$theta <- as.numeric(levels(summary_change$theta))[summary_change$theta]
summary_change$degrees <- as.numeric(levels(summary_change$degrees))[summary_change$degrees]
sub2 <- subset(summary_change, thres == 3)
p2 <- ggplot(sub2, aes(x = factor(theta), y = degrees, group = factor(drift), color= factor(theta))) + 
  facet_grid(cols = vars(drift), rows = vars(accuracy), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs( x='Hz', y ="angle degrees", color = "Hz") +
  theme_linedraw() +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = -60, size = 4, color = "black") +
  #theme(text = element_text(size=15)) +
  theme_linedraw() +
  geom_hline(aes(yintercept = 0)) 
p2 

##############################
## RT Correct vs uncorrect  ##
##############################
# Creating a new variable for the next plots
data$distribution <- as.factor(data$accuracy)
data$distribution <- mapvalues(data$distribution, from = c(0, 1), to = c("Error", "Correct"))

## Plot 4.1: RT Correct vs Error Distribution ~ Hz
mu4.1 <- ddply(data, c("distribution", "Hz"), summarise, grp.mean = mean(rt))
p4.1 <- ggplot(data, aes(x = rt, color = distribution) ) + geom_density() +
  geom_vline(data = mu4.1, aes(xintercept=grp.mean, color=distribution),linetype="dashed", size = 1) + # Add mean lines
  facet_wrap(~ Hz, labeller = label_both) +
  labs(title = ' RT Correct vs Error Distribution ~ frequency (Hz)', x='reaction time') 
p4.1

## Plot 4.2: RT Correct vs Error Distribution ~ Gd, threshold
mu4.2 <- ddply(data, c("distribution", "Gd", "Hz"), summarise, grp.mean = mean(rt))
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
#####################
## Mean Accuracy  ##
####################
# new colors for some of the next plots
myPalette <- colorRampPalette(brewer.pal(9, "PuBu"))
sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0.6, 1))

# Plot 5.1: Mean Accuracy ~ Hz, Gd, threshold
sub5.1 <- subset(summary, threshold == 4)
p5.1 <- ggplot(sub5.1, aes(x=Hz, y = ACC, color=Gd, group = Gd)) + 
  facet_wrap( ~ threshold, labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right     geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1) + 
  labs(title = 'Mean Accuracy', y ="accuracy") +
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)+ 
  theme(text = element_text(size=15)) 
p5.1

############################
## Accuracy Distribution  ##
############################

# Plot 6.1: Accuracy Distribution ~ Gd, threshold
p6.1 <- ggplot(data, aes(x = accuracy, fill = ..prop..)) + 
  sc +
  facet_grid(cols = vars(Gd), rows = vars(thres), labeller = label_both) + 
  geom_bar(aes(y = ..prop..), colour="black") +
  labs(title = 'Accuracy Distribution ~ gamma distance (Gd), threshold', x='accuracy') +
  scale_x_continuous(breaks=c(0,1)) +
  ylim(0,1) +
  #geom_line(aes(y= ..prop..),stat= "count", position = position_dodge(0.1), color = "magenta4", size = 1) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.1)
p6.1

# Plot 6.2: Accuracy Distribution ~ Hz
p6.2 <- ggplot(data, aes(x = accuracy, fill = ..prop..)) + 
  sc +
  facet_wrap(~ Hz, labeller = label_both) +
  geom_bar(aes(y = ..prop..), colour="black") +
  labs(title = 'Accuracy Distribution ~ frequency (Hz)', x='accuracy') +
  scale_x_continuous(breaks=c(0,1)) +
  ylim(0,1) +
  #geom_line(aes(y= ..prop..),stat= "count", position = position_dodge(0.1), color = "magenta4", size = 1) +
  geom_text(aes( label = scales::percent(round(..prop..,2)),
                 y= ..prop.. ), stat= "count", vjust = -.1)
p6.2

###################
## Accuracy ISI  ##
###################

# Plot 7.1: accuracy - ISI (Delta frequencies)
isi7.1 <- subset(summary.isi, Hz < 4)
isi7.1$Hz <- factor(isi7.1$Hz, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
p7.1 <- ggplot(isi7.1, aes(x=isi, y = ACC, color=Hz, group = Hz)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Delta)', y='accuracy', color ="delta") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p7.1

# Plot 7.2: accuracy - ISI (Theta frequencies)
isi7.2 <- subset(summary.isi, Hz < 8 & Hz > 3)
isi7.2$Hz <- factor(isi7.2$Hz, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
p7.2 <- ggplot(isi7.2, aes(x=isi, y = ACC, color=Hz, group = Hz)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Theta)', y='accuracy', color ="theta") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p7.2

# Plot 7.3: accuracy - ISI (Alpha frequencies)
isi7.3 <- subset(summary.isi, Hz < 13 & Hz > 7)
isi7.3$Hz <- factor(isi7.3$Hz, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
p7.3 <- ggplot(isi7.3, aes(x=isi, y = ACC, color=Hz, group = Hz)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Theta)', y='accuracy', color ="alpha") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p7.3

###############################
## Parameter Search - drift  ##
###############################
# Plot 8.1: Drift v ~ Hz, Gd, threshold
sub8.1 <- subset(summary, threshold == 4)
p8.1 <- ggplot(sub8.1, aes(x=Hz, y = v, color= Gd, group = Gd)) + 
  facet_grid(cols =vars(threshold), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Drift v', y='drift v', x='Hz') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p8.1
ggsave(paste0('~/repos/github_sync-model/Results/Collapsing-Bounds/drift-performance-bump','.jpg'),plot=p8.1, units = 'in', width=14, height =6)



# Plot 8.2: Drift v ~  Hz, Gd, threshold == 3 (example plot)
sub8.2 <- subset(summary, threshold == 3)
p8.2 <- ggplot (sub8.2, aes(x=Hz, y = v, color= Gd, group = Gd)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Drift v ~ frequency (Hz), gamma distance (Gd), threshold', y='drift v') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p8.2

##################################
## Parameter Search - boundary  ##
##################################
# Plot 9.1: Boundary a ~ Hz, Gd, threshold
sub9.1 <- subset(summary, threshold == 4)
  p9.1 <- ggplot(sub9.1, aes(x= Hz, y = a, color= Gd, group = Gd)) + 
    facet_grid(cols = vars(threshold), labeller = label_both) +
    geom_point() + 
    geom_line(position = position_dodge(0.1)) + 
    labs(title = 'Boundary a', y='boundary a') +
    #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
    #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.215, size = 4, color = "black") +
    theme(text = element_text(size=15)) 
  p9.1
  ggsave(paste0('~/repos/github_sync-model/Results/Collapsing-Bounds/boundary-performance-bump','.jpg'),plot=p9.1, units = 'in', width=14, height =6)

# Plot 9.1: Boundary a ~ Hz, Gd, threshold
p9.2 <- ggplot(summary, aes(x= Hz, y = a, color= Gd, group = Gd)) + 
  facet_grid( rows = vars(threshold), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Boundary a ~ frequency (Hz), gamma distance (Gd), threshold', y='boundary a') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.215, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p9.2
  ggsave(paste0('~/repos/github_sync-model/Results/Collapsing-Bounds/boundary-threshold-map','.jpg'),plot=p9.2, units = 'in', width=12, height =7)
    

###########################################
## Parameter Search - Nondecision Time   ##
###########################################
# Plot 10.1: Nondecision time Ter ~ Hz, Gd, threshold
sub10.1 <- subset(summary, threshold == 4)
p10.1 <- ggplot(sub10.1, aes(x=Hz, y = Ter, color=Gd, group = Gd)) + 
  facet_grid(cols = vars(threshold), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter', y='Ter') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 230, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p10.1

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

# Plot 11.1: correlation Rt - accuracy ~ Gd, thres
# shows a three-way correlation
p11.1 <- ggplot(summary, aes(x=RT, y = ACC, color=factor(Hz), group = threshold)) + 
  geom_point(shape = 1,colour = "black", size = 2, position = "jitter") +
  geom_point(size = 1, position = "jitter") + 
  labs(title = 'Correlation RT - Accuracy', y='accuracy', x = "RT", color = "Hz") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 180, label.y = 0.4, size = 4, color = "black", position = position_dodge((width =300))) +
  theme(text = element_text(size=15)) 
p11.1


# Plot 11.2: Correlation drift v - boundary a
p11.2 <- ggplot(summary, aes(x=v, y = a, color=factor(Hz))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - boundary a', y='boundary a', x = "drift v", color = "Hz") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.05, label.y = 0.20, size = 5, color = "black") +
  theme(text = element_text(size=15)) 
p11.2

# Plot 11.3: Correlation drift v - Ter
p11.3 <- ggplot(summary, aes(x=v, y = Ter, group = factor(threshold), color=factor(Hz))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Ter', y='Ter', x = "drift v", color = "Hz") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.12, label.y = 100, size = 5, color = "black", position = position_dodgev(height=250)) +
  theme(text = element_text(size=15)) 
p11.3

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


sub1.2 <- subset(summary.isi, Hz == 4 & Gd == 1 & threshold == 3)
p2 <- ggplot(sub1.2, aes(x = isi, y = v, group = factor(Hz), color= factor(Hz))) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point(position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs(title = '', x='isi', y ="drift v", color = "Hz") +
  theme_linedraw()
p2 

##########
## FFT  ##
##########
frequencies <- seq(2,10) # 2 to 10 Hz
threshs <- 4 #seq(1,4)
Gds <-  5 #seq(1,5)
df = vector()
# ISI are separated by 50ms (1/0.050) = 20Hz 
fs <- 20
for(freq in frequencies){
  for(thr in threshs){
    for(gd in Gds){
      sub <- subset(summary.isi, Hz == freq & Gd == gd & threshold == thr)
      sub.select <-select(sub, RT, ACC, v, a, Ter)
      f.results <- vector()  
      for (i in names(sub.select)){ # looping through the columns
        y <- fft(sub.select[,i]) # selecting the values of the column   
        # 11 values, zero frequency + 5 positive frequencies + 5 negative frequencies
        y.tmp <- abs(y)
        # selecting the zero frequency + 5 positive frequencies 
        y.ampspec <- y.tmp[1:(length(y)/2 + 1)]
        # calculating the full frequency amplitude ([positive + negative frequencies])
        y.ampspec[2:(length(y)/2 + 1)] <- y.ampspec[2:(length(y)/2 + 1)] * 2
        f <- seq(from=0, to=fs/2, length=length(y)/2+1)
        y.ampspec <- y.ampspec[-1]
        f <- f[-1]
        f <- f[which.max(y.ampspec)]
        f.results <- c(f.results, f)
      }
      df <- rbind(df, c(freq, thr, gd, f.results[1], f.results[2], f.results[3], f.results[4], f.results[5]))
    }
  }
}
colnames(df) <- c("MFC_Hz", "thres", "Gd", "RT_Hz", "ACC_Hz", "v_Hz", "a_Hz", "Ter_Hz")
df <- data.frame(df)

p1 <- ggplot(df, aes(x = MFC_Hz, y = v_Hz)) + 
  geom_point(position = position_dodge(0.1)) +
  labs(title = 'Threshold = 4, Gamma Distance = 5', x='MFC Hz', y = "drift v Hz", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 5, label.y = 0.18, size = 4, color = "black") 

p2 <- ggplot(df, aes(x = MFC_Hz, y = RT_Hz)) + 
  geom_point(position = position_dodge(0.1)) +
  labs(title = '', x='MFC Hz', y = "RT Hz", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 5, label.y = 0.18, size = 4, color = "black") 

p3 <- ggplot(df, aes(x = MFC_Hz, y = ACC_Hz)) + 
  geom_point(position = position_dodge(0.1)) +
  labs(title = '', x='MFC Hz', y = "Accuracy Hz", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 5, label.y = 0.18, size = 4, color = "black") 

p4 <- ggplot(df, aes(x = MFC_Hz, y = a_Hz)) + 
  geom_point(position = position_dodge(0.1)) +
  labs(title = '', x='MFC Hz', y = "boundary a Hz", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 5, label.y = 0.18, size = 4, color = "black") 

p5 <- ggplot(df, aes(x = MFC_Hz, y = Ter_Hz)) + 
  geom_point(position = position_dodge(0.1)) +
  labs(title = '', x='MFC Hz', y = "Ter Hz", color = "Hz") +
  theme_linedraw() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 5, label.y = 0.18, size = 4, color = "black") 

p9 <- grid.arrange(p1, p2, p3, p4, p5)
p10 <- grid.arrange(p6, p7, p8, p9)
ggsave(paste0('~/Desktop/fft_super_plot','.jpg'),plot=p10, units = 'in', width=15, height =20)


write.csv(df,"~/Desktop/fft-analysis.csv", row.names = FALSE)


#######################################
## Plotting - figures for the paper  ##
#######################################
# Gamma Variability Plot 

# new labs for the plots
threshold.labs <- c("response threshold = 4")
names(threshold.labs) <- c(4)

sub1.1 <- subset(summary, threshold == 4)
p1.1 <- ggplot(sub1.1, aes(x=Hz, y = RT, color=Gd, group = Gd)) + 
  facet_wrap( ~ threshold, labeller = labeller(threshold = threshold.labs)) +
  labs(title = 'Mean response time (MRT)', y='MRT', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=RT-RT_ci, ymax=RT+RT_ci), width=.1) +
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub1.2 <- subset(summary, threshold == 4)
p1.2 <- ggplot(sub1.2, aes(x=Hz, y = ACC, color=Gd, group = Gd)) + 
  facet_wrap( ~ threshold, labeller = labeller(threshold = threshold.labs)) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right     geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1) + 
  labs(title = 'Probability correct (Pc)', y ="Pc",x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)+ 
  theme(text = element_text(size=20),
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub1.3 <- subset(summary, threshold == 4)
p1.3 <- ggplot(sub1.3, aes(x=Hz, y = v, color= Gd, group = Gd)) + 
  facet_grid(cols =vars(threshold), labeller = labeller(threshold = threshold.labs)) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Drift v', y='v',  x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=20),
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub1.4 <- subset(summary, threshold == 4)
p1.4 <- ggplot(sub1.4, aes(x= Hz, y = a, color= Gd, group = Gd)) + 
  facet_grid(cols = vars(threshold), labeller = labeller(threshold = threshold.labs)) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Boundary a', y='a',  x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.215, size = 4, color = "black") +
  theme(text = element_text(size=20),
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub1.5 <- subset(summary, threshold == 4)
p1.5 <- ggplot(sub1.5, aes(x=Hz, y = Ter, color=Gd, group = Gd)) + 
  facet_grid(cols = vars(threshold), labeller = labeller(threshold = threshold.labs)) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter', y='Ter', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  #stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 230, size = 4, color = "black") +
  theme(text = element_text(size=20),
        plot.margin=unit(c(1,1,1,1), "cm")) 

p1 <- ggarrange(p1.1, p1.2, p1.3, p1.4, p1.5, ncol = 2, nrow = 3)
  ggsave(paste0('/home/rien/repos/github_sync-model/Results/gv_behavioral_plot','.png'),plot=p1, units = 'in', width=12, height =12)

# Performance Bumps  Plot
sub2 <- subset(summary, threshold == 4)
p2.1 <- ggplot(sub2, aes(x=Hz, y = RT, color=Gd, group = Gd)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point(size = 4) + 
  labs(title = 'Mean response time (MRT), response threshold = 4', y='MRT', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  #geom_vline(data = sub1.1, aes(xintercept = as.numeric(Gd)), linetype = "longdash", size = 1) +
  geom_line(position = position_dodge(0.1), size = 2) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=RT-RT_ci, ymax=RT+RT_ci), width=1) +
  theme(text = element_text(size=20),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  scale_x_continuous(breaks = seq(from = 1, to = 9, by = 2))

p2.2 <- ggplot(sub2, aes(x=Hz, y = ACC, color=Gd, group = Gd)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point(size = 4) + 
  #geom_vline(data = sub5.1, aes(xintercept = as.numeric(Gd)), linetype = "longdash", size = 1) +
  geom_line(position = position_dodge(0.1), size = 2) + # errorbars overlap,  move them .05 to the left and right     geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1) + 
  labs(title = 'Probability correct (Pc), response threshold = 4', y='Pc', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=1)+ 
  theme(text = element_text(size=20),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  scale_x_continuous(breaks = seq(from = 1, to = 9, by = 2))

p2.3 <- ggplot(sub2, aes(x=Hz, y = v, color= Gd, group = Gd)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point(size = 4) + 
  #geom_vline(data = sub5.1, aes(xintercept = as.numeric(Gd)), linetype = "longdash", size = 1) +
  geom_line(position = position_dodge(0.1), size = 2) + 
  labs(title = 'Drift v, response threshold = 4', y='v', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  theme(text = element_text(size=20),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  scale_x_continuous(breaks = seq(from = 1, to = 9, by = 2))

p2.4 <- ggplot(sub2, aes(x= Hz, y = a, color= Gd, group = Gd)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point(size = 4) +
  #geom_vline(data = sub5.1, aes(xintercept = as.numeric(Gd)), linetype = "longdash", size = 1) +
  geom_line(position = position_dodge(0.1), size = 2) + 
  labs(title = 'Boundary a, response threshold = 4', y='a', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  theme(text = element_text(size=20),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  scale_x_continuous(breaks = seq(from = 1, to = 9, by = 2))

p2.5 <- ggplot(sub2, aes(x=Hz, y = Ter, color=Gd, group = Gd)) + 
  facet_grid(cols = vars(Gd), rows = vars(threshold), labeller = label_both) +
  geom_point(size = 4) +
  #geom_vline(data = sub5.1, aes(xintercept = as.numeric(Gd)), linetype = "longdash", size = 1) +
  geom_line(position = position_dodge(0.1), size = 2) + 
  labs(title = 'Nondecision time Ter, response threshold = 4', y='Ter', x = 'MFC Theta (Hz)', color = 'Gamma\nVariability\n(Hz)') +
  theme(text = element_text(size=20),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  scale_x_continuous(breaks = seq(from = 1, to = 9, by = 2))

p2 <- ggarrange(p2.1, p2.2, p2.3, p2.4, p2.5, ncol = 1, nrow = 5)
ggsave(paste0('/home/rien/repos/github_sync-model/Results/gv_overall_plot','.png'),plot=p2, units = 'in', width=12, height =13)

# Plot 3: threshold maps to boundary
myPalette <- colorRampPalette(brewer.pal(9, "Reds"))
sc <- scale_colour_gradientn(colours = myPalette(9), limits=c(0.18, 0.201))

sub3.1 <- subset(summary, threshold != 2 & Gd == 4  & Hz != 1 & Hz != 2)
p3.1 <- ggplot(sub3.1, aes(x= as.numeric(threshold), y = a, color = as.factor(Hz), group = Hz)) + 
  facet_grid(cols = vars(Hz), labeller = label_both) +
  geom_point(size = 4) + 
  labs(title = 'Boundary a', y='a', x = 'response threshold', color = "MFC\nTheta (Hz)") +
  theme(text = element_text(size=25)) +
  scale_y_continuous(breaks = seq(from = 0.1850, to = 0.2, by = 0.005)) +
  scale_x_continuous(breaks = seq(from = 3, to = 5, by = 1)) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 3, label.y =0.195, size =5, color = "black") +
  theme(text = element_text(size=20),
        panel.spacing = unit(1.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())
p3.1

# Plot 4: Gamma variability RT distributions
sub4.1 <- subset(data, thres == 3 & isi == 5)
p4.1 <- ggplot(sub4.1, aes(x = rt, fill = Gd)) + 
  facet_grid(cols = vars(Gd), labeller = label_both) +
  geom_density() +
  labs(title = 'Response time (RT) distribution', x='response time', fill = "Gamma\nVariability\n(Hz)") +
  theme(text = element_text(size=20),
        panel.spacing = unit(1.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())
p4.1

# Plot 5: Gamma variability RT distributions
gv.labs <- c("Gv = 1", "Gv = 2", "Gv = 3", "Gv = 4", "Gv = 5")
names(gv.labs) <- c(1,2,3,4,5)

sub5.1 <- subset(data, thres == 5 & isi == 3)
p5.1 <- ggplot(sub5.1, aes(x = rt, fill = Hz)) + 
  facet_grid(rows = vars(Hz), cols = vars(Gd), labeller = labeller(Gd = gv.labs)) +
  geom_density() +
  labs(title = 'Response time (RT) distribution', x='response time', fill = "MFC\ntheta\n(Hz)") +
  theme(text = element_text(size=20),
        panel.spacing = unit(1.5, "lines"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  scale_x_continuous(breaks = seq(from = 200, to = 600, by = 400)) +
  scale_y_continuous(breaks = seq(from = 0.002, to = 0.006, by = 0.004)) 
p5.1
ggsave(paste0('/home/rien/repos/github_sync-model/Results/figures_paper/figures-done/rt_distributions_gv_theta','.png'),plot=p5.1, units = 'in', width=12, height =17)


mfc.labs <- c("Theta = 1 Hz", "Theta = 2 Hz", "Theta = 3 Hz", "Theta = 4 Hz","Theta = 5 Hz", "Theta = 6 Hz", "Theta = 7 Hz", "Theta = 8 Hz", "Theta = 9 Hz")
names(mfc.labs) <- c(1,2,3,4,5,6, 7, 8, 9)

su6.1 <- subset(data, thres == 5 & isi == 3)
p6.1 <- ggplot(sub6.1, aes(x = rt, fill = Gd)) + 
  facet_grid(rows = vars(Hz), cols = vars(Gd), labeller = labeller(Hz =  mfc.labs)) +
  geom_density() +
  labs(title = 'Response time (RT) distribution', x='response time', fill = "Gamma\nVariability\n(Hz)") +
  theme(text = element_text(size=18),
        panel.spacing = unit(1.5, "lines"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_continuous(breaks = seq(from = 200, to = 600, by = 400)) +
  scale_y_continuous(breaks = seq(from = 0.002, to = 0.006, by = 0.004)) 
ggsave(paste0('/home/rien/repos/github_sync-model/Results/figures_paper/figures-done/rt_distributions_gv_theta','.png'),plot=p5.1, units = 'in', width=12, height =17)


# Plot 7: response threshold
thres.labs <- c("Resp. threshold = 3", "Resp. threshold = 4", "Resp. threshold = 5")
names(thres.labs) <- c(3,4,5)

sub7.1 <- subset(data, Gd == 3 & isi == 3 & thres != 2)
p7.1 <- ggplot(sub7.1, aes(x = rt)) + 
  facet_grid(rows = vars(thres), labeller = labeller(thres =  thres.labs)) +
  geom_density() +
  labs(title = 'Response time (RT) distribution', x='Response time (RT)', fill = "Gamma\nVariability\n(Hz)") +
  theme(text = element_text(size=18),
        panel.spacing = unit(1.5, "lines"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_continuous(breaks = seq(from = 200, to = 600, by = 400)) +
  scale_y_continuous(breaks = seq(from = 0.002, to = 0.006, by = 0.004)) 
ggsave(paste0('/home/rien/repos/github_sync-model/Results/figures_paper/figures-done/thres_distributions','.png'),plot=p7.1, units = 'in', width=12, height =12)

# Plot 8: response threshold plots
gamma.labs <- c("Gamma variability = 3")
names(gamma.labs) <- c(3)

sub8.1 <- subset(summary, Gd == 3 & threshold != 2)
p8.1 <- ggplot(sub8.1, aes(x=Hz, y = RT, color=threshold, group = threshold)) + 
  facet_wrap( ~ Gd, labeller = labeller(Gd = gamma.labs)) +
  labs(title = 'Mean response time (MRT)', y='MRT', x = 'MFC Theta (Hz)', color = 'Response\nthreshold') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=RT-RT_ci, ymax=RT+RT_ci), width=.1) +
  geom_point() + 
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm"))

sub8.2 <- subset(summary, Gd == 3 & threshold != 2)
p8.2 <- ggplot(sub8.2, aes(x=Hz, y = ACC, color=threshold, group = threshold)) + 
  facet_wrap( ~ Gd, labeller = labeller(Gd = gamma.labs)) +
  labs(title = 'Probability correct (Pc)', y ="Pc", x = 'MFC Theta (Hz)', color = 'Response\nthreshold') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)+ 
  geom_point() + 
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub8.3 <- subset(summary, Gd == 3 & threshold != 2)
p8.3 <- ggplot(sub8.3, aes(x=Hz, y = v, color=threshold, group = threshold)) + 
  facet_wrap( ~ Gd, labeller = labeller(Gd = gamma.labs)) +
  labs(title = 'Drift v', y='v', x = 'MFC Theta (Hz)', color = 'Response\nthreshold') +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_point() + 
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub8.4 <- subset(summary, Gd == 3 & threshold != 2)
p8.4 <- ggplot(sub8.4, aes(x=Hz, y = a, color=threshold, group = threshold)) + 
  facet_wrap( ~ Gd, labeller = labeller(Gd = gamma.labs)) +
  labs(title = 'Boundary a', y='a', x = 'MFC Theta (Hz)', color = 'Response\nthreshold') +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm")) 

sub8.5 <- subset(summary, Gd == 3 & threshold != 2)
p8.5 <- ggplot(sub8.5, aes(x=Hz, y = a, color=threshold, group = threshold)) + 
  facet_wrap( ~ Gd, labeller = labeller(Gd = gamma.labs)) +
  labs(title = 'Nondecision time Ter', y='Ter', x = 'MFC Theta (Hz)', color = 'Response\nthreshold') +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm")) 

p8 <- ggarrange(p8.1, p8.2, p8.3, p8.4, p8.5, ncol = 2, nrow = 3)
ggsave(paste0('/home/rien/repos/github_sync-model/Results/thres_behavioral_plot','.png'),plot=p8, units = 'in', width=12, height =12)

# Plot 9:
acc.labs <- c("error distribution", " correct distribution")
names(acc.labs) <- c(0, 1)
sub9.1 <- subset(data, thres == 5 & Gd == 3)
p9.1 <- ggplot(sub9.1, aes(x = rt, color= as.factor(Hz))) + 
  facet_grid(cols = vars(accuracy), labeller = labeller(accuracy = acc.labs)) +
  geom_density() +
  labs(title = 'Response time distribution', x='Response time (RT)', color = "MFC\ntheta\n(Hz)") +
  theme(text = element_text(size=20), 
        plot.margin=unit(c(1,1,1,1), "cm")) 
p9.1