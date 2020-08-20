# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-07-03
Description: This file contains all the functions that I use to analyse the simulated data from the sync model. 
"""

######################
# Importing packages #
######################
import numpy as np
import pandas as pd
from matplotlib import pyplot as pl
import scipy.stats as sci

################
#  Main Code   #
################
path = "/home/rien/repos/github_sync-model/Results/model_recovery_results/"

# Response threshold recovery plot
#path to the fitting results file
file = "{0}mfc_recovery_goalfunction1.csv".format(path)
pd_data = pd.read_csv(file, index_col=None, header=0)
# calculating correlation
x = pd_data['mfc'].values
y = pd_data['mfc_result'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)

# plotting
fig, ax = pl.subplots()
ax.scatter(x, y)
ax.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax.set_xlabel("Target MFC theta", fontsize = 18); ax.set_ylabel("Recovered MFC theta", fontsize = 16)
ax.set_title("MFC theta recovery plot (n = {0})".format(pd_data.shape[0]), fontsize = 16)
ax.legend(loc='upper left', fontsize = 14)
pl.setp(ax.get_xticklabels(), fontsize=14)
pl.setp(ax.get_yticklabels(), fontsize=14)


# MFC theta recovery plot
#path to the fitting results file
file = "{0}mfc_thres_recovery_unimodal.csv".format(path)
pd_data = pd.read_csv(file, index_col=None, header=0)
x = pd_data['thres'].values
y = pd_data['thres_results2'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
# plotting
fig, ax = pl.subplots()
ax.scatter(x, y)
ax.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 3)))
ax.set_xlabel("Target response threshold", fontsize = 18); ax.set_ylabel("Recovered response threshold", fontsize = 16)
ax.set_title("Response threshold recovery plot (n = {0})".format(pd_data.shape[0]), fontsize = 16)
ax.legend(loc='best', fontsize = 14)
pl.setp(ax.get_xticklabels(), fontsize=14)
pl.setp(ax.get_yticklabels(), fontsize=14)


