# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-04-03
Description: This file compares the MFC amplitude of two frequencies (you can add more or less). 
			 The reson for checking them is that we have to check if the damping parameter in the sync-model.py is 
			 keeping the amplitude of all frequencies between -0.5 and 0.5, if not this results in a confound when
			 analysing the impact of MFC frequencies on task performance. 
"""

######################
# Importing packages #
######################
# 3rd party packages
import numpy as np
import matplotlib.pyplot as pl

# Some issue with the newer version of numpy to load .npz 
# Solution source: https://stackoverflow.com/questions/55890813/how-to-fix-object-arrays-cannot-be-loaded-when-allow-pickle-false-for-imdb-loa
# save np.load to a local variable 
np_load_old = np.load
# modify the default parameters of np.load
np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)

##############
# Functions #
#############
 
def loading_data(freq, drift, thres):
	# creating the filename, this file is an output of the sync-model.py
	file_name = "EEG_Data_simulation_sub0_thetaFreq" + str(freq) + ".00Hz_thresh" + str(thres) + "_drift" + str(drift) + ".0_256Hz.npz"
	# load in the simulated EEG data 
	data = np.load(file_name)
	# unpack the file 
	files = data.files
	for file in files: 
		EEG = data[file]
	# EEG is a numpy.ndarray, containing a dictionary (size is 1)
	EEG = EEG.item() # getting the dictionary out of the ndarray
	MFC = EEG['MFC'] # getting the MFC activation
	return MFC


###################
# Loading in Data #
###################
MFC1 = loading_data(1, 2, 4) # 1 Hz, drift 2, thres 4
MFC20 = loading_data(20, 2, 4) # 20 Hz, drift 2, thres 4
MFCs = [MFC1, MFC20]

###############################
# Plotting the MFC activation #
###############################
# frequencies laying far away from each other on the frequency spectrum are best to be plotted separatly
for MFC in MFCs:
	pl.figure() # drawing the canvas
	pl.xlabel("time")
	pl.ylabel("amplitude")
	pl.plot(MFC[:, 0]) # plotting the first excitatory unit