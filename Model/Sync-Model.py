#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:38:24 2019
@authors: Pieter Verbeke, Mehdi Senoussi, and Rien Sonck
"""
import numpy as np
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

# transformation from hertz (the MFC frequency) to MFC
# damping parameter to keep a relatively similar amplitude
# (tested from 2 to 15 Hz)
def hertz_to_damp(freq):
    return .0184*freq + .0084


# Model function
def Model_sim(Threshold=5, drift=2, Nsubjects=1, theta_freq = 5, save_eeg = False):
    # Threshold = 5
    # drift = 2
    # Nsubjects = 1
    """
    Model simulation that writes away two files
        @Param Threshold, the response threshold
        @Param drift, drift rate of Neurons in different layers, Note! this is actually the inverse of drift rate in the drift diffusion model
        @Param Nsubjects, how many subjects you want to simulate data for
        
        @File csv-file, should ressemble the behavioral data file you get after testing a participant
        @File npy-file, contains simulated EEG data
    """

    # timing of the experiment
    srate = 500                                               # sampling rate per second (500 sample points each second)
    Preinstr_time = int(.2 * srate)                           # pre-instruction time 1000 sample points (0.2s)
    Instr_time = int(.2 * srate)                              # instruction presentation 1000 sample points (0.2s)
    Prep_time = (np.arange(1.7,2.2,.05) * srate).astype(int)  # ISI ranging from 1.7s to 2.2s, multiplying it by the sampling rate to get the amount of sample points we have for each ISI
    Stim_time = int(.05 * srate)                              # Stimulus presentation of 25 timepoints (0.05s)
    Resp_time = 1 * srate                                     # max response time 350 sample points (1s)
    #Resp_time = .7 * srate                                   # max response time ... sample points (0.7s)
    FB_time = int(.1 * srate)                                 # Feedback presentation 50 sample points (0.1s)
    ITI=(np.arange(1,1.9,.25) * srate).astype(int)            # Inter-trial interval (ITI) ranging from 0s to 0.075s (from 0 to 37 sample points)
    Response_deadline = 1 * srate                             # Response deadline 350 sample points (0.7s)
    # Response_deadline = 0.7 * srate

    # The maximum possible trial sample points 
    total_trial_samples = (Preinstr_time + Instr_time + max(Prep_time) + Stim_time + Resp_time + FB_time + max(ITI)).astype(int)  

    # variables for randomization
    nInstr = 4                                        # number of instructions
    nTilts = 2                                        # number of tilt directions
    nSides = 2                                        # number of stimuli locations
    nStim = nTilts * nSides                           # number of stimuli in total
    nResp = 4                                         # number of responses
    nReps = 5                                         # number of replications
    UniqueTrials = nInstr * nStim * len(Prep_time)    # number of different unique trials
    total_trial_amount = UniqueTrials * nReps         # Total amount of trials

    ###########################
    #    Processing Module   #
    ##########################
    nNodes = nStim + nResp                          # total model nodes = stimulus nodes + response nodes
    r2max = 1                                       # max amplitude
    Cg_1 = (30/srate) * 2 * np.pi                   # Coupling gamma waves, for the stimulus nodes
    Cg_2 = Cg_1 + (drift/srate) * 2 * np.pi         # Coupling gamma waves with frequency difference of 2 Hz, for the response nodes
    damp = 0.3                                      # damping parameter, e.g. OLM cells that damp the gamma amplitude
    decay = 0.9                                     # decay parameter
    noise = 0.05                                    # noise parameter

    Phase = np.zeros((nNodes,2,total_trial_samples,total_trial_amount))             # phase neurons, each node has two phase neurons, we update it each timestep, based on the sample rate of each trial
    Rate = np.zeros((nNodes, total_trial_samples, total_trial_amount))               # rate neurons, each node has one rate neuron

    # Weights initialization
    W = np.ones((nStim,nResp))*0.5
    W[(0,2),1] = 0.1
    W[(0,2),3] = 0.1
    W[(1,3),0] = 0.1
    W[(1,3),2] = 0.1

    #########################
    #    Integrator Module  #
    #########################
    Integr = np.zeros(shape = [nResp, total_trial_samples, total_trial_amount]);              # inhibitory weights inducing competition
    inh = np.ones((nResp,nResp))*-0.01
    for i in range(nResp):
        inh[i,i] = 0

    cumul = 1
    # setting collapsing bound function from [Palestro, Weichart, Sederberg, & Turner (2018)]
    t_thresh = np.linspace(0, 1, 1000) # there are 500 timepoints 
    #a_thresh = Threshold; k_thresh = 2; a_prime_thresh = .0; lamb_thresh = .35 # see eq. 1 of Palestro et al., 2018
    #Threshold_byTime = a_thresh - (1 - np.exp(-(t_thresh/lamb_thresh)**k_thresh)) * ((a_thresh/2.) - a_prime_thresh)
    Threshold_byTime =  np.ones(len(t_thresh)) * Threshold # without lapsing bounds

    #######################
    #    Control Module   #
    ######################
    # theta_freq = 5
    r2_MFC=1                                        #radius MFC
    Ct=(theta_freq/srate)*2*np.pi                   #coupling theta waves
    damp_MFC=.1                                     #damping parameter MFC
    acc_slope=10                                    #MFC slope parameter, is set to -5 in equation (7) of Verguts (2017)
                                                    #(steepness of burst threshold)
                                                    
    MFC = np.zeros((2,total_trial_samples,total_trial_amount))      # MFC phase units, two phase neurons
    Be=0                                            				#bernoulli (rate code MFC)

    LFC = np.zeros((nInstr,total_trial_amount))     # LFC stores information for each instruction for each trial
    LFC_sync = np.zeros((nInstr,4))
    LFC_sync[0,:]=[0,1,4,5]                         # LL  sync left stimulus nodes with left hand nodes
    LFC_sync[1,:]=[2,3,6,7]                         # RR  sync right stimulus nodes with right hand nodes
    LFC_sync[2,:]=[0,1,6,7]                         # LR  sync left stimulus nodes with right hand nodes
    LFC_sync[3,:]=[2,3,4,5]                         # RL  sync right stimulus nodes with left hand nodes


    tiltrate=.1                                     		  # mean tilt ~1.8 degrees =.2*90/10 
    #Instr_activation=np.diag(np.ones((4)))         		  #Instruction activation matrix
    Stim_activation=np.zeros((nStim,nResp))         		  #Stimulus activation matrix
    Stim_activation[0,:]=np.array([1,0,1,0])*tiltrate         #Activate 2 stimuli with left tilt (LL)
    Stim_activation[1,:]=np.array([0,1,0,1])*tiltrate         #Activate 2 stimuli with right tilt(RR)
    Stim_activation[2,:]=np.array([1,0,0,1])*tiltrate         #Activate left stimulus with left tilt and right with right tilt
    Stim_activation[3,:]=np.array([0,1,1,0])*tiltrate         #Activate left stimulus with right tilt and right stimulus with left tilt


    for sub in range(Nsubjects): 
                
        # Randomization for instructions, tilt of stimuli and ISI's
        # let's say: 1 = LL       left stim, left resp   | two times left tilt
        #            2 = RR       right stim, right resp | two time right tilt
        #            3 = LR       left stim, right resp  | left tilt (left) and right tilt (right)
        #            4 = RL       right stim, left resp  | right tilt (left) and left tilt (right)
        
        
        ##################################
        #    Create a factorial design  #
        #################################
        Instr = np.repeat(range(nInstr), nStim * len(Prep_time)) # Repeat the instructions (nInstr: 0-4) for the ISI's of each stimulus and put it into an array
        Stim = np.tile(range(nStim), nInstr * len(Prep_time)) # Repeat the stimuli for each instruction, total amount of stimuli
        Preparation = np.floor(np.array(range(UniqueTrials))/(nStim))%len(Prep_time) # Preparation Period, 11 levels 
        Design = np.column_stack([Instr, Stim, Preparation]) # Create an array that has a stack of lists, each list contains instruction, stimulus and a preparation period
        Design = np.tile(Design,(nReps,1)) # Repeat the design nReps
        np.random.shuffle(Design) # shuffle the design making it have a random order
        Design = Design.astype(int)
        
        #####################################################
        #    Oscillations start point of the phase neurons  #
        #####################################################
        start = np.random.random((nNodes,2))          # Draw random starting points for the two phase neurons of each node
        start_MFC = np.random.random((2))             # MFC phase neurons starting point
        # assign starting values
        Phase[:,:,0,0] = start
        MFC[:,0,0] = start_MFC

        #################################
        #            Records           #
        ################################

        Hit = np.zeros((total_trial_samples,total_trial_amount))      # Hit record, check for the sampeling points of each trial
        RT = np.zeros((total_trial_amount))                           # RT record, 
        accuracy = np.zeros((total_trial_amount))                     # Accuracy record
        Instruct_lock = np.zeros((total_trial_amount))                # Instruction onset record
        Stim_lock = np.zeros((total_trial_amount))                    # Stimulus onset record
        Response_lock = np.zeros((total_trial_amount))                # Response onset record 
        resp = np.ones((total_trial_amount)) * -1                     # Response record
        preparatory_period = np.zeros((total_trial_amount))  
        sync = np.zeros((nStim, nResp, total_trial_amount))           # Sync record between the stimuli and the responses on each trial               


        ############################################
        #            Preinstruction  Period        #
        ###########################################
        
        time = 0
        for trial in range(total_trial_amount): # for every trial in total amount of trials
            
            # FIRST STEP: copying over the phase values of the previous trial to the current trial
            if trial > 0:                                       ### index 0 of the phase neurons are already assigned random starting values, starting points are end points of previous trials
                Phase[:,:,0,trial] = Phase[:,:,time,trial-1]    ### Taking both phase neurons of all the nodes of the current trial and setting it equal to the phase neurons of the previous triaL
                MFC[:,0,trial] = MFC[:,time,trial-1]            ### Taking both phase neurons of the MFC of the current trial and setting it equal to the phase neurons of the previous trial 
        
            # SECOND STEP: updating phase code units each sample point in the preinstruction period
            ## Pre-instruction time = no stimulation and no bursts
            for time in range(Preinstr_time): # looping across the sample points of the pre-instruction time
                ## Cg_1 and Cg_2 are the oscillating gamma frequencies, stimulus and response nodes have different gamma frequencies!
                ## Ct the oscillating frequency of the MFC phase neurons
                ## r2max is the radius (amplitude) of a pair of inhibitory and excitatory neurons
                ## damp is the damping value acting on the excitatory (e.g. OLM cells)
                Phase[0:nStim, : , time + 1, trial] = phase_updating (Neurons=Phase[0:nStim, :, time, trial], Radius=r2max, Damp=damp, Coupling=Cg_1, multiple=True)            ### updating the stimulus nodes
                Phase[nStim:nNodes, : , time + 1, trial]  = phase_updating(Neurons=Phase[nStim:nNodes, : , time,trial], Radius=r2max, Damp=damp, Coupling=Cg_2, multiple=True)  ### updating the response nodes
                MFC[:, time+1, trial] = phase_updating(Neurons = MFC[:, time, trial], Radius=r2_MFC, Damp=damp_MFC, Coupling=Ct, multiple=False)                                ### updating the MFC node
            
            # THIRD STEP: Showing the instructions --> Phase reset
            t = time                    ### setting the current sample point of the preinstruction time on the current trial
            Instruct_lock[trial] = t    ### setting the instruction onset of the current trial
            
            ## phase reset of the MFC phase neurons due to the instruction
            MFC[: , t, trial] = np.ones((2)) * r2_MFC 
            
        ##########################################
        #            Preparatory Period          #
        ##########################################
            
            ## Instruction presentation and preparation period
            ## start syncing but no stimulation yet
            preparatory_period[trial] = Prep_time[Design[trial,2]] ### Set the preparatory period (ISI), use the Preparation variable (which is randomized and 11 levels) to select a Prep_time
            
            for time in range(t , int(t + Instr_time + int(preparatory_period[trial]))): ### looping of the sample points of the ISI + instruction time period
                
                # FIRST Step: set the LFC for the current trial
                LFC[Design[trial,0], trial] = 1 ### set the LFC to 1 for the instruction that is shown on this trial
                
                # SECOND STEP: updating phase code units each sample point in the preparatory period
                Phase[0:nStim, :, time + 1, trial] = phase_updating(Neurons=Phase[0:nStim, :, time, trial], Radius=r2max, Damp=damp, Coupling=Cg_1, multiple=True)         ### updating the stimulus nodes
                Phase[nStim:nNodes, :, time + 1, trial]=phase_updating(Neurons=Phase[nStim:nNodes, :, time, trial], Radius=r2max, Damp=damp, Coupling=Cg_2, multiple=True) ### updating the response nodes
                MFC[:, time+1, trial] = phase_updating(Neurons=MFC[:, time, trial], Radius=r2_MFC, Damp=damp_MFC, Coupling=Ct, multiple=False)                             ### updating the MFC node
            
                # THIRD STEP: Rate code MFC neuron activation is calculated by a bernoulli process, start syncing
                Be = 1 / (1 + np.exp(-acc_slope * (MFC[0,time,trial]-1))) ### Equation (7) in Verguts (2017)
                prob = np.random.random()
                    
                # FOURTH STEP: 
                if prob < Be:       
                    
                    Hit[time, trial] = 1
                    Gaussian = np.random.normal(size=[1,2])
                    for Ins in range(nInstr):
                        if LFC[Ins,trial]: ### Checks which of the 4 instruction is set to 1 in the LFC
                            for nodes in LFC_sync[Ins,:]: ### take the 4 nodes associated with the current instruction
                                Phase[int(nodes), :, time + 1, trial] = decay * Phase[int(nodes), :, time, trial] + Gaussian ### Update the nodes that the LFC selected
            
            t=time
            Stim_lock[trial]=t
            
            ##########################################
            #            Responds Period             #
            ##########################################
            
            # Response period: syncing bursts and rate code stimulation
            
            while resp[trial] == -1 and time < t + Response_deadline: ### while the response of this trial is still equal to -1 (no answer has been given)
                
                time += 1
                
                # FIRST STEP: updating phase code units of processing module
                Phase[0:nStim,:,time+1,trial]=phase_updating(Neurons=Phase[0:nStim,:,time,trial], Radius=r2max, Damp=damp, Coupling=Cg_1, multiple=True)
                Phase[nStim:nNodes,:,time+1,trial]=phase_updating(Neurons=Phase[nStim:nNodes,:,time,trial], Radius=r2max, Damp=damp, Coupling=Cg_2, multiple=True)
                MFC[:,time+1,trial]=phase_updating(Neurons=MFC[:,time,trial], Radius=r2_MFC, Damp=damp_MFC, Coupling=Ct, multiple=False)
                
                # SECOND STEP: bernoulli process in MFC rate
                Be = 1/(1+np.exp(-acc_slope*(MFC[0,time,trial]-1)))
                prob = np.random.random()
                    
                # THIRD STEP: Burst 
                if prob<Be:       
                    
                    Hit[time,trial]=1;
                    Gaussian=np.random.normal(size=[1,2])
                    for Ins in range(nInstr):
                        if LFC[Ins,trial]:
                            for nodes in LFC_sync[Ins,:]:               
                                Phase[int(nodes),:,time+1,trial] = decay * Phase[int(nodes), :, time, trial] + Gaussian
                    
                    
                # FOURTH STEP: updating rate code units
                Rate[0:nStim, time, trial] = Stim_activation[Design[trial,1],:]*(1/(1+np.exp(-5*Phase[0:nStim,0,time,trial]-0.6))) ### Updating ratecode units for the stimulus nodes
                Rate[nStim:nNodes, time, trial] = np.matmul(Rate[0:nStim, time, trial],W)*(1/(1+np.exp(-5*Phase[nStim:nNodes,0,time,trial]-0.6))) ### Updating ratecode units for the response nodes
                Integr[:, time+1, trial] = np.maximum(0, Integr[:, time, trial]+cumul*Rate[nStim:nNodes, time, trial]+np.matmul(inh,Integr[:, time, trial]))+noise*np.random.random((nResp))
                
                
                for i in range(nResp):
                    if Integr[i, time+1, trial] > Threshold_byTime[time-t]: # collapsing bounds
                        resp[trial]=i
                        Integr[:, time+1, trial] = np.zeros((nResp))
                    
            RT[trial]=(time-t)*(1000/srate)
            t=time
            Response_lock[trial]=t
            
            if Design[trial,0]==0:
                if (Design[trial,1]==0 or Design[trial,1]==2 ) and resp[trial]==0:
                    accuracy[trial]=1
                elif (Design[trial,1]==1 or Design[trial,1]==3 ) and resp[trial]==1:
                    accuracy[trial]=1
                else:
                    accuracy[trial]=0
        
            if Design[trial,0]==1:
                if (Design[trial,1]==0 or Design[trial,1]==3 ) and resp[trial]==2:
                    accuracy[trial]=1
                elif (Design[trial,1]==1 or Design[trial,1]==2 ) and resp[trial]==3:
                    accuracy[trial]=1
                else:
                    accuracy[trial]=0
        
            if Design[trial,0]==2:
                if (Design[trial,1]==0 or Design[trial,1]==2 ) and resp[trial]==2:
                    accuracy[trial]=1
                elif (Design[trial,1]==1 or Design[trial,1]==3 ) and resp[trial]==3:
                    accuracy[trial]=1
                else:
                    accuracy[trial]=0
        
            if Design[trial,0]==3:
                if (Design[trial,1]==0 or Design[trial,1]==3 ) and resp[trial]==0:
                    accuracy[trial]=1
                elif (Design[trial,1]==1 or Design[trial,1]==2 ) and resp[trial]==1:
                    accuracy[trial]=1
                else:
                    accuracy[trial]=0
        
            
            for time in range(t, t+ FB_time+ ITI[int(np.round(np.random.random()*3))]):
                
                #updating phase code units of processing module
                Phase[0:nStim,:,time+1,trial]=phase_updating(Neurons=Phase[0:nStim,:,time,trial], Radius=r2max, Damp=damp, Coupling=Cg_1, multiple=True)
                Phase[nStim:nNodes,:,time+1,trial]=phase_updating(Neurons=Phase[nStim:nNodes,:,time,trial], Radius=r2max, Damp=damp, Coupling=Cg_2, multiple=True)
                MFC[:,time+1,trial]=phase_updating(Neurons=MFC[:,time,trial], Radius=r2_MFC, Damp=damp_MFC, Coupling=Ct, multiple=False)
                    
            for st in range(nStim):
                for rs in range(nResp):
                    sync[st,rs, trial]=np.corrcoef(Phase[st,0,int(Stim_lock[trial]):int(Response_lock[trial]),trial],Phase[nStim+rs,0,int(Stim_lock[trial]):int(Response_lock[trial]),trial])[0,1]   
        
        Trials=np.arange(total_trial_amount)
        Design=np.column_stack((Trials, Design, resp, accuracy, RT, Instruct_lock, Stim_lock, Response_lock))
        Column_list='trial,instr,stim,isi,response,accuracy,rt,instr_onset,stim_onset,resp_onset'
        #Column_list_2='Visual 1, Visual 2, Visual 3, Visual 4, Motor 1, Motor 2, Motor 3, Motor 4, MFC'
        filename_behavioral='Behavioral_Data_simulation_sub%i_thetaFreq%.2fHz_thresh%i_drift%.1f' % (sub, theta_freq, Threshold, drift)
        np.savetxt(filename_behavioral+'.csv', Design, header=Column_list, delimiter=',',fmt='%.2f')
        
        if save_eeg:
	        Phase_ds = sig.resample(Phase, int( total_trial_samples/2.), axis = 2)
	        MFC_ds = sig.resample(MFC, int( total_trial_samples/2.), axis = 1)
	        Rate_ds = sig.resample(Rate, int( total_trial_samples/2.), axis = 1)
	        Integr_ds = sig.resample(Integr, int( total_trial_samples/2.), axis = 1)

	        EEG_data = {'Phase':Phase_ds[:,0,:,:], 'MFC':MFC_ds[0,:,:], 'Rate':Rate_ds, 'Integr':Integr_ds}
	        filename_EEG='EEG_Data_simulation_sub%i_thetaFreq%.2fHz_thresh%i_drift%.1f_256Hz' % (sub, theta_freq, Threshold, drift)
	        np.savez(filename_EEG+'.npz', EEG_data)

    #return np.mean(accuracy)


# compute surrogates
import time

###########################
# Generate Model Behavior #
###########################


drifts = np.arange(1, 5)
threshs = np.arange(3, 5)
parallel, my_cvstime, _ = par.parallel_func(Model_sim, n_jobs = -1, verbose = 40)

for d in drifts:
    for thr in threshs:
        print('thresh %i, drift %.1f' % (thr, d))
        t = time.time()
        parallel(my_cvstime(Threshold = thr, drift = d, Nsubjects=1, theta_freq = theta, save_eeg = False) for theta in np.arange(2, 16))
        print('\ttime taken: %.2fmin' % ((time.time() - t) / 60.))
