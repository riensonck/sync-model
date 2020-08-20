# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-06-26
Description: This file plots the doubled edged RTs that I'm using in as input for fitting the sync model. 
"""

######################
# Importing packages #
######################
import pandas as pd
import numpy as np
from matplotlib import pyplot as pl
from Analyse_Functions import read_in_data, plot_double_RT, plot_RT_difference

########################
# Value initialization #
########################
# path to the folder containing the simulated data
path = '/home/rien/repos/github_sync-model/Data/1s-Response-Time-Data/'
# Manually pick the values that you want to plot
# simulated data
mfc_freq = np.arange(4, 9, 2)
thres = np.arange(2, 3, 1) 
# target data
mfc_freq_target = 8
thres_target = 2
gd_target = 2

###################
# Reading in data #
###################
# Create dataframes
pd_data = read_in_data(mfc_freq, thres, path) # method to read in multiple files
target = '{0}Behavioral_Data_simulation_sub0_thetaFreq{1}.00Hz_thresh{2}_drift{3}.0.csv'.format(path, mfc_freq_target, thres_target, gd_target)
pd_targetData = pd.read_csv(target)
pd_targetData['mfc'] = mfc_freq_target
pd_targetData['thres'] = thres_target

############
# Plotting #
############
rows = len(mfc_freq)
cols = len(thres)

############
# Figure 1 #
############
fig, axs = pl.subplots(rows, cols) # double-edged histograms
# if rows or cols onlycontain 1 value the notation is axs[i], if rows and cols have 
# multiple values the notation is axs[i,j]
if rows == 1 or cols == 1:
	for i in range(max(rows,cols)):
		if rows  == 1:
			pd_subset = pd_data[(pd_data["mfc"] == mfc_freq[0]) & (pd_data["thres"] == thres[i])]
		else: 
			pd_subset = pd_data[(pd_data["mfc"] == mfc_freq[i]) & (pd_data["thres"] == thres[0])]
		plot_double_RT(pd_subset, axs[i])
else: 
	for col in range(len(thres)):
		for row in range(len(mfc_freq)):
			pd_subset = pd_data[(pd_data["mfc"] == mfc_freq[row]) & (pd_data["thres"] == thres[col])]
			plot_double_RT(pd_subset, axs[row,col])
fig.tight_layout(h_pad=0)

############
# Figure 2 #
############
fig2, axs2 = pl.subplots(rows, cols) # double-edged difference histograms
if rows == 1 or cols == 1:
	for i in range(max(rows,cols)):
		if rows  == 1:
			pd_subset = pd_data[(pd_data["mfc"] == mfc_freq[0]) & (pd_data["thres"] == thres[i])]
		else: 
			pd_subset = pd_data[(pd_data["mfc"] == mfc_freq[i]) & (pd_data["thres"] == thres[0])]
		plot_RT_difference(pd_subset, pd_targetData, axs2[i])
else: 
	for col in range(len(thres)):
		for row in range(len(mfc_freq)):
			pd_subset = pd_data[(pd_data["mfc"] == mfc_freq[row]) & (pd_data["thres"] == thres[col])]
			plot_RT_difference(pd_subset, pd_targetData, axs2[row,col])
fig2.subplots_adjust(hspace = 1)
fig2.suptitle("Target: threshold: {0}, MFC freq: {1}".format(
	pd_targetData["thres"].unique()[0],
	pd_targetData["mfc"].unique()[0]))