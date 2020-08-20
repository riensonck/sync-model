# -*- coding: utf-8 -*-

"""
Author: Mehdi Senoussi and Rien Sonck
Date: 2020-07-25
Description: This file analyses the participants data
"""
######################
# Importing packages #
######################
import numpy as np
from matplotlib import pyplot as pl
import pandas as pd
import glob, os
from scipy import stats
from mne import parallel as par
from mne.stats import fdr_correction
from scipy import signal as sig
from statsmodels.stats.anova import AnovaRM
from statsmodels.formula.api import ols
from scipy import ndimage

###############
#  Functions  #
###############
def ezdiff(rt, correct, s = 1.0):
	""" 
    Simple EZ-diffusion model based on Wagenmakers et al, 2007

    Parameters
    ----------
    rt      : integer
              response time in seconds 
    correct : integer
              correct (1) or wrong (0) answer
    s       :
    """ 
	if sum(correct==1)==0:
		correct = np.zeros(len(correct))
		correct[0]=1
	logit = lambda p:np.log(p/(1-p))

	assert len(rt)>0
	assert len(rt)==len(correct)
	
	assert np.max(correct)<=1
	assert np.min(correct)>=0
	
	pc=np.mean(correct)
  
	assert pc > 0
	
	# subtract or add 1/2 an error to prevent division by zero
	if pc==1.0:
		pc=1 - 1/(2*len(correct))
	if pc==0.5:
		pc=0.5 + 1/(2*len(correct))
	MRT=np.mean(rt[correct==1])
	VRT=np.var(rt[correct==1])

	assert VRT > 0
	
	r=(logit(pc)*(((pc**2) * logit(pc)) - pc*logit(pc) + pc - 0.5))/VRT
	v=np.sign(pc-0.5)*s*(r)**0.25
	a=(s**2 * logit(pc))/v
	y=(-1*v*a)/(s**2)
	MDT=(a/(2*v))*((1-np.exp(y))/(1+np.exp(y)))
	t=MRT-MDT
	
	return([a,v,t])


########################
# Value initialization #
########################
# path to raw data
path = '/home/rien/repos/github_sync-model/Data/participant_data/'
file = '{0}behav_data_34obs.npy'.format(path)
insttxts = np.array(['LL', 'LR', 'RL', 'RR'])
fix_thresh_indeg = 1.5
n_t = 11 # there are 11 ISIs
step = 0.05
inst_diff_order = np.array([3, 0, 2, 1], dtype=np.int) # instructions order by accuracy (higher to lower)
data = []
obs_all = np.array([1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39])
n_obs = len(obs_all)

# data is a list of dictionnaries, each dictionnary is one participant
oad = np.load(file, allow_pickle=True)

# compute average accuracy and rts across whole experiment for each observer
corr_all = np.zeros(len(obs_all))
rts_med_all = np.zeros(len(obs_all))
for oind, obs_i in enumerate(obs_all):
	corr_all[oind] = oad[oind]['corr_resp'].mean()
	rts_med_all[oind] = np.median(oad[oind]['rts'])

# compute accuracy, RT and DDM params per instruction (across all ISIs)
a_inst = np.zeros([n_obs, 4])
v_inst = np.zeros([n_obs, 4])
t_inst = np.zeros([n_obs, 4])
corr_all_inst = np.zeros([n_obs, 4])
rts_med_all_inst = np.zeros([n_obs, 4])
for obs_i in np.arange(n_obs):
	for inst in np.arange(4):
		inst_mask = oad[obs_i]['instructs']==inst
		a_inst[obs_i, inst], v_inst[obs_i, inst], t_inst[obs_i, inst] =\
			ezdiff(oad[obs_i]['rts'][inst_mask], oad[obs_i]['corr_resp'][inst_mask], s = 1.0)
		corr_all_inst[obs_i, inst] = oad[obs_i]['corr_resp'][inst_mask].mean()
		rts_med_all_inst[obs_i, inst] = np.median(oad[obs_i]['rts'][inst_mask])

# make a dataframe of a, v, t, rt, acc across all instructions ans subjects
data = np.concatenate((rts_med_all_inst, corr_all_inst), axis = 1)
data = np.concatenate((data, a_inst ), axis = 1)
data = np.concatenate((data, v_inst ), axis = 1)
data = np.concatenate((data, t_inst ), axis = 1)
pd_data_inst = pd.DataFrame(data, columns = ('rt_inst1', 'rt_inst2', 'rt_inst3', 'rt_inst4',
											'acc_inst1', 'acc_inst2', 'acc_inst3', 'acc_inst4',
											'a_inst1', 'a_inst2', 'a_inst3', 'a_inst4',
											'v_inst1', 'v_inst2', 'v_inst3', 'v_inst4',
											't_inst1', 't_inst2', 't_inst3', 't_inst4'))
trial_count = []
for i in range(len(oad)):
	pd_data = pd.DataFrame(oad[i])
	count = len(pd_data[pd_data['instructs'] == 0])
	trial_count.append(pd.DataFrame({'subject': i, 'trials': count}, index = [i]))
pd_trials = pd.concat(trial_count, axis=0, ignore_index=True)
pd_trials.sort_values('trials', ascending=False)

####################################################
# 						PLOT 					   #
####################################################
# Plot 1: across instructions
titles = ['Bound', 'Drift Rate', 'Non-dec time', 'RT', 'Accuracy']
fig, axs = pl.subplots(2, 3)
for data_ind, data in enumerate([a_inst, v_inst, t_inst, rts_med_all_inst, corr_all_inst]):
	ax = axs.flatten()[data_ind]
	toplot = data[:, inst_diff_order]
	ax.errorbar(x = np.arange(4), y=toplot.mean(axis=0),
		yerr=toplot.std(axis=0)/(toplot.shape[0]**.5), color='k', fmt = 'o',
	linestyle='-', mfc = 'w', ms=8, mec='k', elinewidth=3, mew=2,
	capthick=0, ecolor ='k', linewidth=3, zorder=2)
	# 1way RM-ANOVA
	# prepare stuff for the ANOVA
	subjs = np.repeat(obs_all, 4)
	conds = np.tile(np.arange(1, 5), n_obs)
	df = pd.DataFrame({'obs':subjs, 'conds':conds,
				'perf':toplot.flatten()})
	aovrm = AnovaRM(data=df, depvar='perf', subject='obs', within=['conds'])
	res = aovrm.fit()
	ax.set_title('%s\nF=%.2f, pval=%.4f' %\
		(titles[data_ind], np.float(res.anova_table['F Value']), res.anova_table['Pr > F']),
		fontsize=10)
	ax.set_xticks(np.arange(4)); ax.set_xticklabels(insttxts[inst_diff_order])
	ax.grid()
pl.tight_layout()

# Plot 2: across subjects (choose one instruction)
titles = ['Accuracy', 'RT','Bound', 'Drift Rate', 'Non-dec time']
labels = ['accuracy (%)', 'response time (seconds)', 'bound a', 'drift v', 'nondecision time Ter']
subjects = np.arange(1, 35)
number_of_subplots = int(len(pd_data_inst.columns)/len(insttxts))
instr_index = np.arange(0, 20, 4) # select all the instructions 1
fig2, axs2 = pl.subplots(2, 3)

for i, ax in enumerate(axs2.reshape(-1)[:-1]):
	ax.scatter(subjects, pd_data_inst.iloc[:,instr_index[i]], s = 50)
	#axs.set_title("{0}".format(titles[i]), fontsize=15)
	ax.set_xlabel("subject", fontsize=18); ax.set_ylabel(labels[i], fontsize=18)
	ax.set_xticks(np.arange(0, 36, 5))
	ax.grid()
	pl.setp(ax.get_xticklabels(), fontsize=14)
	pl.setp(ax.get_yticklabels(), fontsize=14)
fig2.suptitle("Inter-individual Variability in Participant Performance", fontsize=20)
fig2.tight_layout(rect=[0, 0.03, 1, 0.95])

# Plot 3: RT distributions
bins = np.arange(0, 1000, 10)
instr = 1
fig3, axs3 = pl.subplots(6, 6)
for i, ax in enumerate(axs3.reshape(-1)[:-2]): # have 36 axs, only need 34
	theta_hist = np.zeros(shape=[2, 1, len(bins)-1])
	data = pd.DataFrame(oad[i])
	# get histograms
	theta_hist = np.histogram(data.loc[(data['instructs'] == instr) & (data['corr_resp'] == 1), 'rts'].values * 1000, bins=bins, normed = True)[0] # error histograms
	# run difference kernel
	theta_histDiff= -np.array([theta_hist[bin_n] - (theta_hist[bin_n-1] + theta_hist[bin_n+1]).mean() for bin_n in np.arange(1, len(bins)-2)]).T
	# smooth
	theta_histDiffSmooth = ndimage.gaussian_filter1d(theta_histDiff, axis = 0, sigma = 2)
	ax.plot(theta_histDiffSmooth)
	ax.set_xticks((0, 25, 50, 75))
	ax.set_xticklabels(labels = ["0", "250", "500", "750"])
	ax.set_xlabel("Response time (ms)", fontsize = 12)
	ax.set_title("Subject: {0}".format(i + 1))
pl.tight_layout()


# Plot 4: across subjects (choose one instruction)
titles = ['RT']
labels = ['response time (seconds)']
subjects = np.arange(1, 35)
number_of_subplots = int(len(pd_data_inst.columns)/len(insttxts))
instr_index = np.arange(0, 20, 4) # select all the instructions 1
fig2, axs2 = pl.subplots(1, 2)
instr = ['rt_inst1', 'rt_inst2']
titles = ['Condition 1', 'Condition 2']

for i, ax in enumerate(axs2):
	ax.scatter(subjects,pd_data_inst[instr[i]] , s = 50)
	ax.set_xlabel("subject number", fontsize=18); ax.set_ylabel(labels[0], fontsize=18)
	ax.set_xticks(np.arange(0, 36, 5))
	ax.set_yticks(np.arange(0.4, 0.6, 0.05))
	ax.axhline(pd_data_inst[instr[i]].mean(), linewidth = 3, label = "Average", color = "darkorange" )
	ax.grid()
	ax.legend()
	ax.set_title(titles[i], fontsize=18)
	pl.setp(ax.get_xticklabels(), fontsize=14)
	pl.setp(ax.get_yticklabels(), fontsize=14)
