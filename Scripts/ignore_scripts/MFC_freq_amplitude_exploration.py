#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:38:24 2019

@authors: Mehdi Senoussi (main author) and Rien Sonck
"""

#########################
#    Importing Modules   #
#########################
# Third party modules
import numpy as np
import matplotlib.pyplot as pl
from mne import parallel as par
from scipy import signal as sig

##################
#    Functions   #
##################

def phase_updating(Neurons=[], Radius=1, Damp=0.3, Coupling=0.3, multiple=True):
    """
    Returns updated value of the inhibitory (I) and excitatory (E) phase neurons
        @Param Neurons, list containing the current activation of the phase code units
        @Param Radius,  bound of maximum amplitude
        @Param Damp, strength of the attraction towards the Radius, e.g. OLM cells reduce  pyramidal cell activation
        @Param Coupling, frequency*(2*pi)/sampling rate
        @Param multiple, True or false statement whether the Neurons array includes more than one node or not
        
        Formula (2) and (3) from Verguts (2017) 
    """
    # updating the phase neurons from the processing module
    if multiple:
        Phase = np.zeros ((len(Neurons[:,0]),2)) # Creating variable to hold the updated phase values ( A zero for each E and I phase neuron of each phase code unit)
        r2 = np.sum(Neurons * Neurons, axis = 1) # calculating the amplitude depending on the activation of the E and I phase neurons
        # updating the E phase neurons, the higher the value of the I neurons (Neurons[:, 1]), the lower the value of the E neurons
        Phase[:,0] = Neurons[:,0] -Coupling * Neurons[:,1] - Damp *((r2>Radius).astype(int)) * Neurons[:,0] 
        # updating the I phase neurons, the higher the value of the E neurons (Neuron[:, 0]), the higher the value of the I neurons
        Phase[:,1] = Neurons[:,1] +Coupling * Neurons[:,0] - Damp * ((r2>Radius).astype(int)) * Neurons[:,1]
    # updating the phase neurons of the MFC
    else:
        Phase = np.zeros((2)) 
        r2 = np.sum(Neurons[0] * Neurons[1], axis = 0) 
        Phase[0] = Neurons[0] -Coupling*Neurons[1] - Damp * ((r2>Radius).astype(int)) * Neurons[0]
        Phase[1] = Neurons[1] +Coupling*Neurons[0] - Damp * ((r2>Radius).astype(int)) * Neurons[1]    
        
    return Phase

def hertz_to_damp(freq):
    """
    transformation from hertz (the MFC frequency) to MFC
    damping parameter to keep a relatively similar amplitude across frequencies
    (tested from 2 to 15 Hz)
    """
    return .0184*freq + .0084

#################################
#    Timing of the Experiment   #
#################################
srate = 500                                               # sampling rate per second (500 sample points each second, so 500 Hz sampling)
Preinstr_time = int(.2 * srate)                           # pre-instruction time 1000 sample points (0.2s)
Instr_time = int(.2 * srate)                              # instruction presentation 1000 sample points (0.2s)
Prep_time = (np.arange(2.7,3.3,.05) * srate).astype(int)  # ISI ranging from 2.7s to 3.25s, multiplying it by the sampling rate to get the amount of sample points we have for each ISI
Stim_time = int(.05 * srate)                              # Stimulus presentation of 25 timepoints (0.05s)
Resp_time = .7 * srate                                    # max response time 350 sample points (0.7s)
FB_time = int(.1 * srate)                                 # Feedback presentation 50 sample points (0.1s)
ITI = (np.arange(0,.1,.025) * srate).astype(int)          # Inter-trial interval (ITI) ranging from 0s to 0.075s (from 0 to 37 sample points)
Response_deadline = .7 * srate                            # Response deadline 350 sample points (0.7s)

# The maximum possible trial sample points, 2286 sample points
total_trial_samples = (Preinstr_time + Instr_time + max(Prep_time) + Stim_time + Resp_time + FB_time + max(ITI)).astype(int)

##########################
#    Experiment Setup   #
#########################
nInstr = 4                                        # number of instructions
nTilts = 2                                        # number of tilt directions
nSides = 2                                        # number of stimuli locations
nStim = nTilts * nSides                           # number of stimuli in total
nResp = 4                                         # number of responses
nReps = 10                                        # number of replications
UniqueTrials = nInstr * nStim * len(Prep_time)    # number of different unique trials, 192 unique trials
total_trial_amount = UniqueTrials * nReps         # Total amount of trials, 1920 trials

fig, axs = pl.subplots(2, 1) # create two subplots

# Will loop through different frequencies to check their amplitude
for theta_freq_ind, theta_freq in enumerate(np.arange(2, 16, 1)):
    Threshold=5; drift=2; Nsubjects=1;  # keeping these constant

    ###########################
    #    Processing Module   #
    ##########################
    nNodes = nStim + nResp                                                # total model nodes in the processing module = stimulus nodes + response nodes
    
    # Each node has 2 phase neurons and one rate neuron
    Phase = np.zeros((nNodes,2, total_trial_samples, total_trial_amount)) # phase neurons: we will update all the sample points of each trial for both phase neurons
    Rate = np.zeros((nNodes, total_trial_samples, total_trial_amount))    # rate neurons

    #######################
    #    Control Module   #
    #######################
    r2_MFC = 1                                       # radius MFC
    Ct = (theta_freq/srate) * 2 * np.pi              # coupling theta waves, how many times does the frequency couple the gamma frequencies 
    # 15Hz, 15 cycles per second. One second is sampled 500 times so 15Hz is sampled 
    # how many times we are sample the frequency peaks eahc second. 
    # formula for angular velocity 
    # how many number of rotations in a period of time are sampled
    # have 15Hz cycles and sample rate is 500Hz, 0.03 

    damp_MFC = hertz_to_damp(theta_freq)  # damping parameter MFC
    MFC_slope = 10                        # steepness of the MFC burst threshold

    MFC = np.zeros((2,total_trial_samples, total_trial_amount))       # MFC phase units, two phase neurons
    MFC[:,0,0] = [.1,.1] # set the first timepoint of the first trial to .1 for both MFC phase neurons
    preparatory_period = np.zeros((total_trial_amount))  

    ############################################
    #            Preinstruction  Period        #
    ###########################################
    time = 0
    trial = 0    
    # FIRST STEP: copying over the phase values of the previous trial to the current trial
    if trial > 0:                                       # starting points of new trials are end points of the previous trial
        MFC[:,0,trial] = MFC[:,time,trial-1]            # Taking both phase neurons of the MFC of the current trial and setting it equal to the phase neurons of the previous trial 
    # SECOND STEP: updating phase code units each sample point in the preinstruction period
    for time in range(Preinstr_time): # looping across the sample points of the pre-instruction time
        MFC[:, time+1, trial] = phase_updating(Neurons = MFC[:, time, trial], Radius = r2_MFC, Damp = damp_MFC, Coupling = Ct, multiple=False)                                ### updating the MFC node
    # THIRD STEP: Showing the instructions --> Phase reset
    t = time                    # setting the current sample point of the preinstruction time on the current trial
    # phase reset of the MFC phase neurons due to the instruction
    MFC[: , t, trial] = np.ones((2)) * r2_MFC 

    ##########################################
    #            Preparatory Period          #
    ##########################################
    # Instruction presentation and preparation period
    preparatory_period[trial] = Prep_time[10] ### Set the preparatory period (ISI), for this script we use a constant ISI
    for time in range(t , int(t + Instr_time + int(preparatory_period[trial]))): ### looping of the sample points of the ISI + instruction time period
        MFC[:, time+1, trial] = phase_updating(Neurons=MFC[:, time, trial], Radius=r2_MFC, Damp=damp_MFC, Coupling=Ct, multiple=False)                             ### updating the MFC node
    t = time            
    # Response period: syncing bursts and rate code stimulation
    while time < t + Response_deadline: # Add the amount of sample points of the response period to the current sample point we are at
        time += 1
        MFC[:,time+1,trial]=phase_updating(Neurons=MFC[:,time,trial], Radius=r2_MFC, Damp=damp_MFC, Coupling=Ct, multiple=False)

    ##########################################
    #            Plotting        #
    ##########################################
    axs[0].plot(MFC[0, :, 0], label='%iHz - %.2fd' % (theta_freq, damp_MFC)) # plotting all the sample points of one trial of one MFC phase neuron
    axs[0].set_ylabel("amplitude")
    fftamp = np.abs(np.fft.fft(MFC[0, :, 0])) # fast fourier transform
    axs[1].plot(fftamp, label='%iHz' % theta_freq)
    axs[1].set_ylabel("magnitude")

axs[0].legend()
