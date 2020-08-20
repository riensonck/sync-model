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

# simple EZ-diffusion model based on Wagenmakers et al, 2007
def ezdiff(rt, correct, s = 1.0):
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

def temporal_surrogates(x, n_surrs, n_t):
	temp_surr = np.zeros(shape = [n_surrs, n_t])
	rand_ind = np.arange(n_t)
	for surr_n in range(n_surrs):
		np.random.shuffle(rand_ind)
		temp_surr[surr_n, :] = x[rand_ind]

	return np.abs(np.fft.fft(temp_surr, axis = 1)[:, 1:int((n_t/2 + 1))])

def temporal_padded_surrogates(x, n_surrs, added_samples_byside, n_t, n_t_pad, n_freq, freqpadmask, assymetry):
	temp_surr = np.zeros(shape = [n_surrs, n_t_pad+assymetry])
	rand_ind = np.arange(n_t)
	for surr_n in range(n_surrs):
		np.random.shuffle(rand_ind)
		rand_obs_vals = x[rand_ind]
		avgvals = rand_obs_vals.mean()
		temp_surr[surr_n, :] = np.hstack([np.repeat(avgvals, added_samples_byside+assymetry),
									rand_obs_vals,
									np.repeat(avgvals, added_samples_byside)])
	return np.abs(np.fft.fft(temp_surr, axis = 1)[:, freqpadmask])

# path to raw data
path = '/home/rien/repos/github_sync-model/Data/Participant_Data/'
file = '{0}behav_data_34obs.npy'.format(path)

insttxts = np.array(['LL', 'LR', 'RL', 'RR'])

fix_thresh_indeg = 1.5

# there were 11 ISIs
n_t = 11
step = 0.05
# instructions order by accuracy (higher to lower)
inst_diff_order = np.array([3, 0, 2, 1], dtype=np.int)


oad = []
obs_all = np.array([1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39])
n_obs = len(obs_all)

# oad is a list of dictionnaries, each dictionnary is one participant
oad = np.load(file, allow_pickle=True)

# Relevant keys of each dictionnary
# 	stim_combis: represents how the gratings were tilted (cw = clockwise, ccw = counter-clockwise):
#		 0 = cw (left stim) - cw (right stim)  /  1 = ccw - cw  /  2 = cw - ccw  /  3 = ccw - ccw
#
#	instructs: represents the instruction
#		0 = LL, 1 = LR, 2 = RL, 3 = RR


# We don't use it but you can convert the dictionnary to a pandas dataframe in case it's useful
pdoad = pd.DataFrame(oad)


# compute average accuracy and rts across whole experiment for each observer
corr_all = np.zeros(len(obs_all))
rts_med_all = np.zeros(len(obs_all))
for oind, obs_i in enumerate(obs_all):
	corr_all[oind] = oad[oind]['corr_resp'].mean()
	rts_med_all[oind] = np.median(oad[oind]['rts'])


a = np.zeros(n_obs)
v = np.zeros(n_obs)
t = np.zeros(n_obs)
corr_all = np.zeros(n_obs)
rts_med_all = np.zeros(n_obs)
for obs_i in np.arange(n_obs):
	a[obs_i], v[obs_i], t[obs_i] = ezdiff(oad[obs_i]['rts'], oad[obs_i]['corr_resp'], s = 1.0)
	corr_all[obs_i] = oad[obs_i]['corr_resp'].mean()
	rts_med_all[obs_i] = np.median(oad[obs_i]['rts'])



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




####################################################
# 						PLOT 					   #
####################################################

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






