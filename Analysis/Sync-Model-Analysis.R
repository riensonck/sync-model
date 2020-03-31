library(ggplot2) # library for plotting
library(ggpubr)
library(dplyr) # library to join dataframes
library(cowplot) # has the function background_grid
library(tidybayes)
library(RColorBrewer) # for color gradient
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
get.vaTer = function(theta, thres, drift, Pc, VRT, MRT, s = 0.1){
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
  
  return(c(theta, thres, drift,round(v,4), round(a, 4), round(Ter,4)))
}

## EZ Diffusion Model
## Code adapted from Eric-Jan Wagenmakers, Han L.J. van der Maas , and Raoul P. P. P. Grasman (2007) 
get.vaTer.isi = function(theta, thres, drift, isi, Pc, VRT, MRT, s = 0.1){
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
    Pc = 0.9999
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
  
  return(c(theta, thres, drift, isi, round(v,4), round(a, 4), round(Ter,4)))
}

#####################
## Importing Data  ##
#####################

setwd("~/Desktop/UGent/Thesis-Study/Data/Generated-Data/") # Setting the working directory 
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

#######################
## Summarizing data  ##
#######################
# Summarizing the RT data
RT_summary <- summarySE(data, measurevar = "rt", groupvars = c("theta", "drift", "thres"))
names(RT_summary) <- c("theta", "drift", "thres", "N", "RT", "RT_sd", "RT_se", "RT_ci")

# Summarizing the accuracy data
ACC_summary <- summarySE(data, measurevar = "accuracy", groupvars = c("theta", "drift", "thres"))
names(ACC_summary) <- c("theta", "drift", "thres", "N", "ACC", "ACC_sd", "ACC_se", "ACC_ci")

# Summarizing and transforming to drift v, boundary a, and nondecision-time Ter
EZ_data <- matrix(ncol=6, nrow = 640)
MRT <- ddply(data, .(theta, drift, thres) , summarize, mean = mean(rt)) # mean response time for correct and uncorrect decisions for each theta
Pc <- ddply(data, .(theta, drift, thres) , summarize, mean = mean(accuracy)) # Pc, proportion of correct decisions
VRT <- ddply(data, .(theta, drift, thres) , summarize, sd = sd(accuracy)) # VRT, variance of response times for correct decisions
for (i in 1:nrow(Pc)){
  vaTer <- get.vaTer(Pc$theta[i], Pc$drift[i], Pc$thres[i], Pc$mean[i], VRT$sd[i], MRT$mean[i])
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

# Plot 2: RT Distribution ~ theta, drift, threshold
data$theta <- factor(data$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
p2 <- ggplot(data, aes(x = rt, color= theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_density() +
  labs(title = ' RT Distribution ~ theta, drift, threshold', x='reaction time', color = "theta") 
p2 

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
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1) + 
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

isi <- subset(summary_isi, theta < 4)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
# Plot 9: accuracy - ISI
p9 <- ggplot(isi, aes(x=isi, y = ACC, color=theta, group = theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Delta)', y='accuracy', color ="delta") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p9

isi <- subset(summary_isi, theta < 8 & theta > 3)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
# Plot 9: accuracy - ISI
p9 <- ggplot(isi, aes(x=isi, y = ACC, color=theta, group = theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Theta)', y='accuracy', color ="theta") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p9

isi <- subset(summary_isi, theta < 13 & theta > 7)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
# Plot 9: accuracy - ISI
p9 <- ggplot(isi, aes(x=isi, y = ACC, color=theta, group = theta)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI (Alpha)', y='accuracy', color ="alpha") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=ACC-ACC_ci, ymax=ACC+ACC_ci), width=.1)
p9



isi_summary <- summarySE(data, measurevar = "accuracy", groupvars = c("isi", "drift", "thres"))

isi <- subset(summary_isi)
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
# Plot 9: accuracy - ISI
p9 <- ggplot(isi_summary, aes(x=isi, y = accuracy, color = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  labs(title = 'Accuracy ~ ISI', y='accuracy', color ="drift") +
  geom_line(position = position_dodge(0.1)) + # errorbars overlap,  move them .05 to the left and right 
  geom_errorbar(aes(ymin=accuracy-ci, ymax=accuracy+ci), width=.1)
p9

###############################
## Parameter Search - drift  ##
###############################
# Plot 10: Drift v ~ theta, drift, threshold
p10 <- ggplot(summary, aes(x=theta, y = v, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Drift v ~ theta, drift, threshold', y='drift v') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.18, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p10

# Plot 11: Boundary a ~ theta, drift, threshold
p11 <- ggplot(summary, aes(x=theta, y = a, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Boundary a ~ theta, drift, threshold', y='boundary a') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 0.215, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p11

# Plot 12: Nondecision time Ter ~ theta, drift, threshold3
sub12 <- subset(summary, thres == 3)
p12 <- ggplot(sub12, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 350, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p12

# Plot 13: Nondecision time Ter ~ theta, drift, threshold4
sub13 <- subset(summary, thres == 4)
p13 <- ggplot(sub13, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 450, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p13

# Plot 14: Nondecision time Ter ~ theta, drift, threshold5
sub14 <- subset(summary, thres == 5)
p14 <- ggplot(sub14, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y =550, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p14

# Plot 15: Nondecision time Ter ~ theta, drift, threshold6
sub15 <- subset(summary, thres == 6)
p15 <- ggplot(sub15, aes(x=theta, y = Ter, color=drift, group = drift)) + 
  facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() + 
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'Nondecision time Ter ~ theta, drift, threshold', y='Ter') +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 1, label.y = 650, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p15

##################################
## Parameter Correlation Plots  ##
##################################

# Plot 16: correlation Rt - accuracy ~ drfit, thres
p16 <- ggplot(summary, aes(x=RT, y = ACC, color=factor(theta), group = thres)) + 
  geom_point(shape = 1,colour = "black", size = 2, position = "jitter") +
  geom_point(size = 1, position = "jitter") + 
  labs(title = 'Correlation RT - Accuracy', y='accuracy', x = "RT", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 200, label.y = 0.4, size = 4, color = "black") +
  theme(text = element_text(size=15)) 
p16

#Plot 17: correlation RT-accuracy ~ theta 18
sub17 <- subset(summary, theta == 18)
p17 <- ggplot(sub17, aes(x=RT, y = ACC, color=factor(thres))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation RT - Accuracy ~ theta 18', y='accuracy', x = "RT", color = "threshold") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 300, label.y = 0.8, size = 4, color = "black") +
  theme(text = element_text(size=15)) +
  facet_wrap(~ drift, labeller = label_both )
 p17
 
 #Plot 18: correlation RT-accuracy ~ theta 6
 sub17 <- subset(summary, theta == 6)
 p18 <- ggplot(sub17, aes(x=RT, y = ACC, color=factor(thres))) + 
   geom_point(shape = 1,colour = "black", size = 2) +
   geom_point(size = 1) + 
   labs(title = 'Correlation RT - Accuracy ~ theta 6', y='accuracy', x = "RT", color = "threshold") +
   stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
   stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 300, label.y = 0.8, size = 4, color = "black") +
   theme(text = element_text(size=15)) +
   facet_wrap(~ drift, labeller = label_both )
 p18


# Plot 19: Correlation drift v - boundary a
p19 <- ggplot(summary, aes(x=v, y = a, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - boundary a', y='boundary a', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.05, label.y = 0.21, size = 5, color = "black") +
  theme(text = element_text(size=15)) 
p19

# Plot 20: Correlation drift v - boundary a
p20 <- ggplot(summary, aes(x=v, y = Ter, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Ter', y='Ter', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.15, label.y = 500, size = 5, color = "black") +
  theme(text = element_text(size=15)) 
p20

# Plot 21: Correlation drift v - boundary a ~ threshold
p21 <- ggplot(summary, aes(x=v, y = Ter, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Ter ~ threshold', y='Ter', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.15, label.y = 500, size = 5, color = "black") +
  theme(text = element_text(size=15))  + 
  facet_wrap(~ thres, labeller = label_both)
p21

# Plot 22: Correlation drift v - accuracy
p22 <- ggplot(summary, aes(x=v, y = ACC, color=factor(theta))) + 
  geom_point(shape = 1,colour = "black", size = 2) +
  geom_point(size = 1) + 
  labs(title = 'Correlation drift v - Accuracy', y='accuracy', x = "drift v", color = "theta") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, color = "black") +
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = 0.05, label.y = 1, size = 5, color = "black") +
  theme(text = element_text(size=15))  
p22

#####################################
## Needing new dataframe with ISI  ##
#####################################
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
  vaTer <- get.vaTer.isi(Pc$theta[i], Pc$drift[i], Pc$thres[i], Pc$isi[i], Pc$mean[i], VRT$sd[i], MRT$mean[i])
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

####################
## ISI Plotting  ##
###################
isi$theta <- factor(isi$theta, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))
# Plot 21: accuracy - ISI
sub21 <- subset(summary_isi, thres == 3 & drift == 1)
p21 <- ggplot(sub20, aes(x=isi, y = v, color=factor(theta), group = factor(theta))) + 
  #facet_grid(cols = vars(drift), rows = vars(thres), labeller = label_both) +
  geom_point() +
  facet_wrap(~ theta, labeller = label_both) +
  geom_line(position = position_dodge(0.1)) + 
  labs(title = 'drift v ~ ISI', y='drift v') 
p21


modS <- lm(v ~ thres + theta + thres*theta, data = summary)
f <- summary(modS)$fstatistic


modS <- glm(subj_choice ~ deltaev + deltaev*lottery_won,family=binomial(link='logit'),data=sf)
insert <- c(as.character(subjects_nonverbal[i]), 
            as.numeric(modS$coefficients[1]), 
            as.numeric(modS$coefficients[2]), 
            as.numeric(modS$coefficients[3]), 
            as.numeric(modS$coefficients[4]))
resmodS2_nonverbal[i,] <- insert
