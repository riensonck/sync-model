#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: Mehdi Senoussi and Rien Sonck
@description: This file generates simulations using the sync-model
@date: 2020/06/28
"""

######################
# Importing packages #
######################
import time
import numpy as np
from Sync_Model import model_sim
from mne import parallel as par

########################
# Value initialization #
########################
# pick the model parameter values which you want to generate simulations with. 
gds = np.arange(4, 5, 1)        # gamma distance parameter
threshs = np.arange(2, 3, 1)    # threshold parameter
mfcs = np.arange(4, 10, 1)      # mfc frequency parameter, often the theta frequency range is used (4Hz-7Hz)

###################
# Data generation #
###################
# running the simulations parallel
parallel, my_cvstime, _ = par.parallel_func(model_sim, n_jobs = -1, verbose = 40)

for gd in gds:
    for thr in threshs:
        print('thresh %i, drift %.1f' % (thr, gd))
        t = time.time()
        parallel(my_cvstime(Threshold = thr, drift = gd, Nsubjects=1, theta_freq = mfc, save_behavioral = True, save_eeg = False, return_dataframe = False) for mfc in mfcs)
        print('\ttime taken: %.2fmin' % ((time.time() - t) / 60.))

# Work in progress
"""
thr = 4
d = 5
Cg_1 = np.concatenate([np.arange(0.2, 0.5, 0.05), np.arange(0.5, 1.01, 0.1)])
Cg_2 = np.concatenate([np.arange(0.25, 0.55, 0.05), np.arange(0.65, 1.15, 0.1)])
for i in range(len(Cg_1)):
    print('Cg_1: {0:.2f}, Cg_2: {1:.2f}'.format(Cg_1[i], Cg_2[i]))
    t = time.time()
    parallel(my_cvstime(Cg_1[i], Cg_2[i], Threshold = thr, drift = d, Nsubjects=1, theta_freq = theta, save_eeg = False) for theta in np.arange(2, 16))
    print('\ttime taken: %.2fmin' % ((time.time() - t) / 60.))
"""