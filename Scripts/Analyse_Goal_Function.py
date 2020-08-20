# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-07-05
Description: This file plots the goal function.
"""

######################
# Importing packages #
######################
import numpy as np
import pandas as pd
import time
import matplotlib
from matplotlib import pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Analyse_Functions import new_goal_func_from_files
# IMPORTANT: In Analyse_Functions.py, change def goal_func(x, *args) to goal_func(x, args) 
# for this file to work

########################
# Value initialization #
########################
# path to the folder containing the target data
path = '/home/rien/repos/github_sync-model/Data/Collapsing-Bounds-Data-New/'
# Manually pick the values that you want to plot
# simulated data
mfc_freq = np.arange(1, 10, 0.5)
thres = np.arange(3, 6, 1)
gds = np.arange(1, 6, 1) 
# target data
mfc_target = 8
thres_target = 3
gd_target = 4
target_file = '{0}Behavioral_Data_simulation_sub0_thetaFreq{1:.2f}Hz_thresh{2}_drift{3}.0.csv'.format(path, mfc_target, thres_target, gd_target)
pd_targetData = pd.read_csv(target_file)

######################################
# Still have to generate data files #
#####################################
results = []
for mfc in mfc_freq: 
  for gd in gds: 
    for thr in thres:
    	x = mfc
    	args = (gd, thr, pd_targetData)
    	t = time.time()
    	error = goal_func(x, args)
    	print('\ttime taken: %.2fmin' % ((time.time() - t) / 60.))
    	d = {"mfc": [mfc], "thres": [thr], "gd": [gd], "error": [error]}
    	data = pd.DataFrame(data=d)
    	results.append(data)
pd_results = pd.concat(results, axis=0, ignore_index=True)
pd_results.to_csv("Error_function_sub0_thetaFreq{0}.00Hz_thresh{1}_drift{2}".format(mfc_target, thres_target, gd_target))

######################################
# Using already generated data files #
######################################
results = []
for mfc in mfc_freq: 
  for gd in gds: 
    for thr in thres:
    	error = new_goal_func_from_files(mfc, thr, gd, target_file, path)
    	d = {"mfc": [mfc], "thres": [thr], "gd": [gd], "error": [error]}
    	data = pd.DataFrame(data=d)
    	results.append(data)
pd_results = pd.concat(results, axis=0, ignore_index=True)
pd_results.to_csv("Error_function_sub0_thetaFreq{0:.2f}Hz_thresh{1}_drift{2}".format(mfc_target, thres_target, gd_target))

#############
# Plotting #
############

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(pd_results['mfc'].values, pd_results['thres'].values, pd_results['error'].values)
ax.set_xlabel("mfc"); ax.set_ylabel("thres"); ax.set_zlabel("error")
ax.set_title("Target: mfc = {0}, thres = {1}, gd = {2}".format(mfc_target, thres_target, gd_target))

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(pd_results['mfc'].values, pd_results['gd'].values, pd_results['error'].values, cmap ="jet")
ax.set_xlabel("mfc"); ax.set_ylabel("gd"); ax.set_zlabel("error")
ax.set_title("Target: mfc = {0}, thres = {1}, gd = {2}".format(mfc_target, thres_target, gd_target))

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(pd_results['thres'].values, pd_results['gd'].values, pd_results['error'].values, cmap ="jet")
ax.set_xlabel("thres"); ax.set_ylabel("gd"); ax.set_zlabel("error")
ax.set_title("Target: mfc = {0}, thres = {1}, gd = {2}".format(mfc_target, thres_target, gd_target))