# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-07-27
Description: This file analyses the participants fits 
"""
######################
# Importing packages #
######################
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy.stats as sci


########################
# Value initialization #
########################
# path to raw data
path = '/home/rien/repos/github_sync-model/Results/fitting_participants_results/'
file = '{0}participant_fits.csv'.format(path)
pd_data = pd.read_csv(file)
subj = 1
subj_data = pd_data[pd_data['subject']== subj]
insttxts = np.array(['LL', 'LR', 'RL', 'RR'])
inst_diff_order = np.array([3, 0, 2, 1])
insttxts_diff_order = np.array(['RR', 'LL', 'RL', 'LR'])

#############
# Plotting #
############
# plot 1
x = np.arange(4)
y = subj_data[['mfc_instr3','mfc_instr0', 'mfc_instr2', 'mfc_instr1']].values[0]
fig, ax = plt.subplots() 
ax.set_title("Subject: {0}".format(subj))
ax.plot(x, y)
ax.scatter(x, y, s = 50, color = "darkorange")
ax.set_xticks(x); ax.set_xticklabels(insttxts_diff_order)
ax.set_xlabel("instructions"); ax.set_ylabel("recovered MFC parameter")

# plot 2
x =  pd_data['v_instr0'].values
y = pd_data['mfc_instr0'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
fig2, ax2 = plt.subplots()
ax2.scatter(x, y, s = 50, color = "darkorange")
ax2.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax2.set_xlabel("Drift v", fontsize = 18); ax2.set_ylabel("Fitted MFC theta", fontsize = 18)
ax2.set_title("Correlation with drift v", fontsize = 18)
ax2.legend(loc='upper left', fontsize = 14)
plt.setp(ax2.get_xticklabels(), fontsize=14)
plt.setp(ax2.get_yticklabels(), fontsize=14)

# plot 3
x =  pd_data['a_instr0'].values
y = pd_data['mfc_instr0'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
fig2, ax2 = plt.subplots()
ax2.scatter(x, y, s = 50, color = "darkorange")
ax2.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax2.set_xlabel("Bound a", fontsize = 18); ax2.set_ylabel("Fitted MFC theta", fontsize = 18)
ax2.set_title("Correlation with bound a", fontsize = 18)
ax2.legend(loc='upper left', fontsize = 14)
plt.setp(ax2.get_xticklabels(), fontsize=14)
plt.setp(ax2.get_yticklabels(), fontsize=14)

# plot 4
x =  pd_data['acc_instr0'].values
y = pd_data['mfc_instr0'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
fig2, ax2 = plt.subplots()
ax2.scatter(x, y, s = 50, color = "darkorange")
ax2.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax2.set_xlabel("Pc", fontsize = 18); ax2.set_ylabel("Fitted MFC theta", fontsize = 18)
ax2.set_title("Correlation with probability correct (Pc)", fontsize = 18)
ax2.legend(loc='upper left', fontsize = 14)
plt.setp(ax2.get_xticklabels(), fontsize=14)
plt.setp(ax2.get_yticklabels(), fontsize=14)

# plot 4
x =  pd_data['t_instr0'].values
y = pd_data['mfc_instr0'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
fig2, ax2 = plt.subplots()
ax2.scatter(x, y, s = 50, color = "darkorange")
ax2.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax2.set_xlabel("Non-decision time", fontsize = 18); ax2.set_ylabel("Fitted MFC theta", fontsize = 18)
ax2.set_title("Correlation with non-decision time", fontsize = 18)
ax2.legend(loc='upper left', fontsize = 14)
plt.setp(ax2.get_xticklabels(), fontsize=14)
plt.setp(ax2.get_yticklabels(), fontsize=14)

# plot 4
x =  pd_data['rt_instr0'].values
y = pd_data['mfc_instr0'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
fig2, ax2 = plt.subplots()
ax2.scatter(x, y, s = 50, color = "darkorange")
ax2.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax2.set_xlabel("MRT", fontsize = 18); ax2.set_ylabel("Fitted MFC theta", fontsize = 18)
ax2.set_title("Correlation with mean response time (MRT)", fontsize = 18)
ax2.legend(loc='upper left', fontsize = 14)
plt.setp(ax2.get_xticklabels(), fontsize=14)
plt.setp(ax2.get_yticklabels(), fontsize=14)

# plot 4
x =  pd_data['eeg_instr0'].values
y = pd_data['mfc_instr0'].values
m, b = np.polyfit(x, y, 1)
r_spearman, p_value = sci.spearmanr(x, y)
fig2, ax2 = plt.subplots()
ax2.scatter(x, y, s = 50, color = "darkorange")
ax2.plot(x, m*x + b, label = "Spearman r = {0}, p-value = {1}".format(np.round(r_spearman,2), np.round(p_value, 2)))
ax2.set_xlabel("Peak frequency", fontsize = 18); ax2.set_ylabel("Fitted MFC theta", fontsize = 18)
ax2.set_title("Correlation with participants' peak frequency", fontsize = 18)
ax2.legend(loc='upper left', fontsize = 14)
plt.setp(ax2.get_xticklabels(), fontsize=14)
plt.setp(ax2.get_yticklabels(), fontsize=14)