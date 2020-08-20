# -*- coding: utf-8 -*-

"""
Author: Mehdi Senoussi and Rien SonckM
Date: 2020-07-26
Description: This file analyses the RT distribution peaks
"""
######################
# Importing packages #
######################
import numpy as np
import pandas as pd
from matplotlib import pyplot as pl
from scipy import ndimage

###############
#  Functions  #
###############
def thetas_distribution_smoothing(path, thetas, threshold, gd, bins, sigma):
	""" 
    Function that applies the difference kernel and then a Gaussian smoothing kernel to the error 
    and correct RT distributions corresponding to each MFC frequency in the model.  

    Parameters
    ----------
    path       : string 
                 path to the folder containing the data
    thetas     : numpy array
                 array containing all the mfc frequency values (here theta frequency)
    threshold  : integer
                 threshold parameter value
    bins       : numpy array
                 array containing the bins for the histograms
    sigma      : integer
                 smoothing of the Gaussian kernel
    """ 
	accuracy = [0, 1]
	thetas_hist_list = []
	all_thetas_hist = np.zeros(shape=[2, len(thetas), len(bins)-1])
	###################
	# Reading in data #
	###################
	for theta_ind, theta_freq in enumerate(thetas):
		data = np.loadtxt(path + 'Behavioral_Data_simulation_sub0_thetaFreq{0}.00Hz_thresh{1}_drift{2}.0.csv'.format(theta_freq, threshold, gd),
			skiprows=1, dtype=np.float, delimiter=',')
		data_err = data[data[:, 5]==0, :]
		data_cor = data[data[:, 5]==1, :]
		# creating histograms
		all_thetas_hist[0, theta_ind, :] = np.histogram(data_err[:, 6], bins=bins, normed = True)[0] # error histograms
		all_thetas_hist[1, theta_ind, :] = np.histogram(data_cor[:, 6], bins=bins, normed = True)[0] # correct histograms
	###################################
	# Difference Kernel and Smoothing #
	###################################
	# run difference kernel
	all_thetas_histDiff_err = -np.array([all_thetas_hist[0, :, bin_n] - (all_thetas_hist[0, :, bin_n-1] + all_thetas_hist[0, :, bin_n+1]).mean() for bin_n in np.arange(1, len(bins)-2)]).T
	all_thetas_histDiff_cor = -np.array([all_thetas_hist[1, :, bin_n] - (all_thetas_hist[1, :, bin_n-1] + all_thetas_hist[1, :, bin_n+1]).mean() for bin_n in np.arange(1, len(bins)-2)]).T
	# smooth
	all_thetas_histDiffSmooth_err = ndimage.gaussian_filter1d(all_thetas_histDiff_err, axis = 1, sigma = sigma)
	all_thetas_histDiffSmooth_cor = ndimage.gaussian_filter1d(all_thetas_histDiff_cor, axis = 1, sigma = sigma)
	return (all_thetas_histDiffSmooth_err, all_thetas_histDiffSmooth_cor)

########################
# Value initialization #
########################
path = '/home/rien/repos/github_sync-model/Data/Collapsing-Bounds-Data-New/'
thetas  = np.arange(4, 10)
thres = np.arange(3, 6)
bins = np.arange(0, 1000, 10)
sigmas = [3, 4, 3]
gds = np.arange(1, 6, 1)
data = [] 

for i in range(len(thres)): # or thres
	for j in range(len(gds)):
		# applying kernels to the error and correct RT distributions
		thetas_hist_list = thetas_distribution_smoothing(path, thetas, thres[i], gds[j], bins, sigmas[i])
		# finding where in the bins the peaks and troughs happen and mark it 
		time_step_change = np.diff(np.sign(np.diff(thetas_hist_list)))
		# extracting the index of the peak and trough bins based on the marker
		err_indices = [np.where(time_step_change[0][i] != 0) for i in range(len(time_step_change[0]))]
		cor_indices = [np.where(time_step_change[1][i] != 0) for i in range(len(time_step_change[1]))]
		# extracting the peaks and trough bin value
		peaks_troughs_err = [thetas_hist_list[0][i][err_indices[i]] for i in range(len(thetas_hist_list[0]))]
		peaks_troughs_cor = [thetas_hist_list[1][i][cor_indices[i]] for i in range(len(thetas_hist_list[1]))]
		# saving the data to a dataframe
		for k in range(len(thetas)):
			d = pd.DataFrame.from_dict({'thres': [thres[i]], 'theta': [thetas[k]], 'gd': [gds[j]], 
							'p1_timing_err': err_indices[k][0][0], 'p1_timing_cor': cor_indices[k][0][0],
							'p1_err': [peaks_troughs_err[k][0]], 'p1_cor': [peaks_troughs_cor[k][0]]})
							#'p2_timing_err': err_indices[k][0][2],
							#'p2_timing_cor': cor_indices[k][0][2],
							#'p2_err': [peaks_troughs_err[k][2]],
							#'p2_cor': [peaks_troughs_cor[k][2]]})
			data.append(d)
pd_data = pd.concat(data, axis=0, ignore_index=True)
# adding the relative height of the peaks (peak 1 - peak 2)
pd_data['rel_height_err_min'] = pd_data['p1_err'] - pd_data['p2_err']
pd_data['rel_height_cor_min'] = pd_data['p1_cor'] - pd_data['p2_cor']
# adding the relative height of the peaks (peak 2 / peak 1)
pd_data['rel_height_err_div'] = pd_data['p2_err'] / pd_data['p1_err']
pd_data['rel_height_cor_div'] = pd_data['p2_cor'] / pd_data['p1_cor']
# adding the relative timing of the peaks (peak 2 timing / peak 1 timing)
pd_data['rel_timing_err_div'] = pd_data['p2_timing_err'] / pd_data['p1_timing_err']
pd_data['rel_timing_cor_div'] = pd_data['p2_timing_cor'] / pd_data['p1_timing_cor']
# adding the relative timing of the peaks (peak 2 timing - peak 1 timing)
pd_data['rel_timing_err_min'] = pd_data['p2_timing_err'] - pd_data['p1_timing_err']
pd_data['rel_timing_cor_min'] = pd_data['p2_timing_cor'] - pd_data['p1_timing_cor']
# adding the relative timing of the peaks ((peak 2 timing - peak 1 timing / peak 1 timing)
pd_data['rel_timing_err'] = (pd_data['p2_timing_err'] - pd_data['p1_timing_err']) / pd_data['p2_timing_err']
pd_data['rel_timing_cor'] = (pd_data['p2_timing_cor'] - pd_data['p1_timing_cor']) / pd_data['p2_timing_cor']

# adding the relative height of the peaks (peak 1 - peak 2)
pd_data['rel_height_1stpeak_div'] = pd_data['p1_err'] / pd_data['p1_cor']
pd_data['rel_height_1stpeak_min'] = pd_data['p1_err'] - pd_data['p1_cor']

############
# Plotting #
############
fig, axs = pl.subplots(1, 2)
for theta_ind, theta_freq in enumerate(thetas):
	axs[0].plot(bins[:-1][1:-1], thetas_hist_list[0][theta_ind], label=theta_freq)
	axs[1].plot(bins[:-1][1:-1], thetas_hist_list[1][theta_ind], label=theta_freq)
axs[0].legend(); axs[1].legend()
axs[0].set_xlabel("Reaction Time"); axs[1].set_xlabel("Reaction Time")
axs[0].set_title("Error Distribution"); axs[1].set_title("Correct Distribution")

# Plotting the absolute height
fig2, axs2 = pl.subplots(1, 2)
for thr in thres: 
	axs2[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "p1_err"].values, label="threshold:" + str(thr))
	axs2[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "p2_cor"].values, label="threshold:" + str(thr))
	axs2[0].set_title("Absolute Height: Error Peak 1")
	axs2[1].set_title("Absoulte Height: Correct Peak 2")
	axs2[0].set_ylim(-.0020, 0.0075); axs2[1].set_ylim(0, 0.0040)
	axs2[0].set_xlabel("Hz"); axs2[1].set_xlabel("Hz")
	axs2[0].legend(); axs2[1].legend() 
fig2.tight_layout()

# Plotting the relative height
fig3, axs3 = pl.subplots(1, 2)
for thr in thres: 
	axs3[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_height_err_min"].values, label="threshold:" + str(thr))
	axs3[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_height_cor_min"].values, label="threshold:" + str(thr))
	axs3[0].set_title("Relative Height: Error Peak 1 - Error Peak 2")
	axs3[1].set_title("Relative Height: Correct Peak 1 - Correct Peak 2")
	#axs3[0].set_ylim(-.0025, 0.0040); axs3[1].set_ylim(0.0010, 0.0040)
	axs3[0].set_xlabel("Hz"); axs3[1].set_xlabel("Hz")
	axs3[0].legend(); axs3[1].legend() 
fig3.tight_layout()

# Plotting the relative height
fig4, axs4 = pl.subplots(1, 2)
for thr in thres: 
	axs4[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_height_err_div"].values, label="threshold:" + str(thr))
	axs4[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_height_cor_div"].values, label="threshold:" + str(thr))
	axs4[0].set_title("Relative Height: Error Peak 2 / Error Peak 1")
	axs4[1].set_title("Relative Height: Correct Peak 2 / Correct Peak 1")
	#axs4[0].set_ylim(-.0025, 0.0040); axs4[1].set_ylim(0.0010, 0.0040)
	axs4[0].set_xlabel("Hz"); axs4[1].set_xlabel("Hz")
	axs4[0].legend(); axs4[1].legend() 
fig4.tight_layout()


# Plotting the timing of the second peak
fig5, axs5 = pl.subplots(1, 2)
for thr in thres: 
	axs5[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "p2_timing_err"].values, label="threshold:" + str(thr))
	axs5[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "p2_timing_cor"].values, label="threshold:" + str(thr))
	axs5[0].set_title("Absolute Timing: Error Peak 2 ")
	axs5[1].set_title("Absolute Timing: Correct Peak 2")
	#axs5[0].set_ylim(-.0025, 0.0040); axs5[1].set_ylim(0.0010, 0.0040)
	axs5[0].set_xlabel("Hz"); axs5[1].set_xlabel("Hz")
	axs5[0].legend(); axs5[1].legend() 
fig5.tight_layout()


# Plotting the relative timing
fig6, axs6 = pl.subplots(1, 2)
for thr in thres: 
	axs6[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_timing_err_div"].values, label="threshold:" + str(thr))
	axs6[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_timing_cor_div"].values, label="threshold:" + str(thr))
	axs6[0].set_title("Relative Timing: Error Peak 2 / Error Peak 1 ")
	axs6[1].set_title("Relative Timing: Correct Peak 2 / Correct Peak 1")
	#axs6[0].set_ylim(-.0025, 0.0040); axs6[1].set_ylim(0.0010, 0.0040)
	axs6[0].set_xlabel("Hz"); axs6[1].set_xlabel("Hz")
	axs6[0].legend(); axs6[1].legend() 
fig6.tight_layout()


# Plotting the relative timing
fig7, axs7 = pl.subplots(1, 2)
for thr in thres: 
	axs7[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_timing_err_min"].values, label="threshold:" + str(thr))
	axs7[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_timing_cor_min"].values, label="threshold:" + str(thr))
	axs7[0].set_title("Relative Timing: Error Peak 2 - Error Peak 1 ")
	axs7[1].set_title("Relative Timing: Correct Peak 2 - Correct Peak 1")
	#axs7[0].set_ylim(-.0025, 0.0040); axs7[1].set_ylim(0.0010, 0.0040)
	axs7[0].set_xlabel("Hz"); axs7[1].set_xlabel("Hz")
	axs7[0].legend(); axs7[1].legend() 
fig7.tight_layout()

# Plotting the relative timing
fig8, axs8 = pl.subplots(1, 2)
for thr in thres: 
	axs8[0].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_timing_err"].values, label="threshold:" + str(thr))
	axs8[1].scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_timing_cor"].values, label="threshold:" + str(thr))
	axs8[0].set_title("Relative Timing: (Error Peak 2 - Error Peak 1) / Error Peak 1")
	axs8[1].set_title("Relative Timing: (Correct Peak 2 - Correct Peak 1) / Error Peak 1")
	#axs8[0].set_ylim(-.0025, 0.0040); axs8[1].set_ylim(0.0010, 0.0040)
	axs8[0].set_xlabel("Hz"); axs8[1].set_xlabel("Hz")
	axs8[0].legend(); axs8[1].legend() 
fig8.tight_layout()

# Plotting the relative height
fig9, axs9 = pl.subplots()
for thr in thres: 
	axs9.scatter(thetas, pd_data.loc[(pd_data["thres"] == thr), "rel_height_1stpeak_div"].values, label="threshold:" + str(thr))
	axs9.set_title("Relative Height: (Error Peak 1 / Correct Peak 1)")
	#axs8[0].set_ylim(-.0025, 0.0040); axs8[1].set_ylim(0.0010, 0.0040)
	axs9.set_xlabel("Hz")
	axs9.legend()

# Plotting the absolute height
fig10, axs10 = pl.subplots(1, 2)
for theta in thetas: 
	axs10[0].scatter(gds, pd_data.loc[(pd_data["theta"] == theta) & (pd_data['thres'] == 3), "p1_err"].values, label="Theta:" + str(theta))
	axs10[1].scatter(gds, pd_data.loc[(pd_data["theta"] == theta) & (pd_data['thres'] == 3), "p1_cor"].values, label="Theta:" + str(theta))
	axs10[0].set_title("Absolute Height: Error Peak 1")
	axs10[1].set_title("Absoulte Height: Correct Peak 2")
	#axs10[0].set_ylim(-.0020, 0.0075); axs10[1].set_ylim(0, 0.0040)
	axs10[0].set_xlabel("Gd"); axs10[1].set_xlabel("Gd")
	axs10[0].legend(); axs10[1].legend() 
fig10.tight_layout()

# Plotting the absolute height
fig11, axs11 = pl.subplots()
for theta in thetas: 
	axs11.scatter(gds, pd_data.loc[(pd_data["theta"] == theta) & (pd_data['thres'] == 3), "rel_height_1stpeak_min"].values, label="Theta:" + str(theta))
	axs11.set_title("Relative Height: (Error Peak 1 / Correct Peak 1)")
	#axs10[0].set_ylim(-.0020, 0.0075); axs10[1].set_ylim(0, 0.0040)
	axs11.set_xlabel("Gd"); 
	axs11.legend();
fig11.tight_layout()