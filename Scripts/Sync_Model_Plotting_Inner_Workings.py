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


Threshold=3; drift=3; Nsubjects=1; theta_freq = 4; save_eeg = False; sim_path = './'

srate = 500                                               # sampling rate per second
Preinstr_time = int(.2 * srate)                                 # pre-instruction time (1s)
Instr_time = int(.2 * srate)                              #  instruction presentation (200 ms)
Prep_time = (np.arange(1.7,2.2,.05) * srate).astype(int)  # ISI ranging from 1700 to 2200 ms, multiplying it by the sampling rate to get how many samples we have for each ISI
Stim_time = int(.05 * srate)                              # Stimulus presentation of 50 ms
Resp_time = .7 * srate                                     # max response time of 1s
FB_time = int(.1 * srate)                                 # Feedback presentation of 500 ms
ITI = (np.arange(0,.1,.025) * srate).astype(int)          # ITI ranging from 1000 ms to 1900 ms
Response_deadline = .7 * srate                             # Response deadline

# max trial time
TotT = (Preinstr_time + Instr_time + max(Prep_time) + Stim_time + Resp_time + FB_time + max(ITI)).astype(int)  

# variables for randomization
nInstr = 4                                        # number of instructions
nTilts = 2                                        # number of tilt directions
nSides = 2                                        # number of stimuli locations
nStim = nTilts * nSides                           # number of stimuli in total
nResp = 4                                         # number of responses
nReps = 20                                        # number of replications
UniqueTrials = nInstr * nStim * len(Prep_time)    # number of different unique trials
Tr = UniqueTrials * nReps                         # Total amount of trials

###########################
#    Processing Module   #
##########################
nNodes = nStim + nResp                          # total model nodes = stimulus nodes + response nodes
r2max = 1                                       # max amplitude
Cg_1 = 0.38                                     # Coupling gamma waves, for the stimulus nodes
Cg_2 = Cg_1 + drift/100                         # Coupling gamma waves with frequency difference of 2 Hz, for the response nodes
damp = 0.3                                      # damping parameter, e.g. OLM cells that damp the gamma amplitude
decay = 0.9                                     # decay parameter
noise = 0.05                                    # noise parameter

Phase = np.zeros((nNodes,2,TotT,Tr))             # phase neurons, each node has two phase neurons, we update it each timestep, based on the sample rate of each trial
Rate = np.zeros((nNodes, TotT, Tr))                        # rate neurons, each node has one rate neuron

# Weights initialization
W = np.ones((nStim,nResp))*0.5
W[(0,2),1] = 0.1
W[(0,2),3] = 0.1
W[(1,3),0] = 0.1
W[(1,3),2] = 0.1

#########################
#    Integrator Module  #
#########################
Integr = np.zeros(shape = [nResp, TotT, Tr]);              # inhibitory weights inducing competition
inh = np.ones((nResp,nResp))*-0.01
for i in range(nResp):
    inh[i,i] = 0

cumul = 1
#Threshold=4

#######################
#    Control Module   #
######################
# theta_freq = 5
r2_MFC=1                                        #radius MFC
Ct=(theta_freq/srate)*2*np.pi                            #coupling theta waves
damp_MFC=.03                                    #damping parameter MFC
acc_slope=10                                    #MFC slope parameter, is set to -5 in equation (7) of Verguts (2017)
                                                #(steepness of burst threshold)
                                                
MFC = np.zeros((2,TotT,Tr))                     # MFC phase units, two phase neurons
Be=0                                            #bernoulli (rate code MFC)

LFC = np.zeros((nInstr,Tr))                     # LFC stores information for each instruction for each trial
LFC_sync = np.zeros((nInstr,4))
LFC_sync[0,:]=[0,1,4,5]                         # LL  sync left stimulus nodes with left hand nodes
LFC_sync[1,:]=[2,3,6,7]                         # RR  sync right stimulus nodes with right hand nodes
LFC_sync[2,:]=[0,1,6,7]                         # LR  sync left stimulus nodes with right hand nodes
LFC_sync[3,:]=[2,3,4,5]                         # RL  sync right stimulus nodes with left hand nodes


tiltrate=.1                                     # mean tilt ~1.8 degrees =.2*90/10 
#Instr_activation=np.diag(np.ones((4)))         #Instruction activation matrix
Stim_activation=np.zeros((nStim,nResp))         #Stimulus activation matrix
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
    """
    TODO: No idea what preparation is doing
    """
    Preparation = np.floor(np.array(range(UniqueTrials))/(nStim))%len(Prep_time) # Preparation Period, 11 levels 
    Design = np.column_stack([Instr, Stim, Preparation]) # Create an array that has a stack of lists, each list contains instruction, stimulus and a preparation period
    Design = np.column_stack([np.tile(Design,(nReps,1)), np.repeat(np.arange(nReps), UniqueTrials)]) # Repeat the design nReps
    np.random.shuffle(Design) # shuffle the design making it have a random order

    Design = Design.astype(int)
    
    #####################################################
    #    Oscillations start point of the phase neurons  #
    #####################################################
    start = np.random.random((nNodes,2))          # Draw random starting points for the two phase neurons of each node
    """
    # TODO: MFC = ACC? Or ACC is a part of the MFC?  
    """
    start_MFC = np.random.random((2))             # Acc phase neurons starting point
    # assign starting values
    Phase[:,:,0,0] = start
    MFC[:,0,0] = start_MFC

    #################################
    #            Records           #
    ################################

    Hit = np.zeros((TotT,Tr))                     # Hit record, check for the sampeling points of each trial
    RT = np.zeros((Tr))                           # RT record, 
    accuracy = np.zeros((Tr))                     # Accuracy record
    Instruct_lock = np.zeros((Tr))                # Instruction onset record
    Stim_lock = np.zeros((Tr))                    # Stimulus onset record
    Response_lock = np.zeros((Tr))                # Response onset record 
    resp = np.ones((Tr)) * -1                     # Response record
    preparatory_period = np.zeros((Tr))  
    sync = np.zeros((nStim, nResp, Tr))           # Sync record between the stimuli and the responses on each trial               


    ############################################
    #            Preinstruction  Period        #
    ###########################################
    
    time = 0
    Tr = 10
    for trial in range(Tr): # for every trial in total amount of trials
        
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
                if Integr[i, time+1, trial]>Threshold:
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




import os
import time as tm
from matplotlib import pyplot as pl


def plot_shadederr(subp, data, error_measure = 'avg_sem', percs = [], color = 'blue',
    x = None, xlim = None, ylim = None, label = None, linestyle = '-', alpha = 1, linewidth = 1):
    error_val_up = np.zeros(data.shape[-1]); error_val_low = error_val_up
    if error_measure != None:
        if error_measure[:3]=='avg':
            curve_val = data.mean(axis=0)
        elif error_measure[:3]=='med':
            curve_val = np.median(data, axis=0)

        if error_measure[3:]=='sem':
            error_val_low = curve_val - data.std(axis=0) #/ np.sqrt(data.shape[0])
            error_val_up = error_val_low + 2*curve_val
        elif error_measure[-4:]=='perc':
            if len(percs)==0: percs = [25,75]
            error_val_low = np.percentile(data, percs[0], axis=0)
            error_val_up = np.percentile(data, percs[1], axis=0)

    subp.plot(x, curve_val, color = color, label = label, alpha = alpha, linestyle = linestyle, linewidth = linewidth)
    if error_measure != None:
        subp.fill_between(x, error_val_up, error_val_low, color = color, alpha = .2)
    pl.grid(); pl.ylim(ylim); pl.xlim(xlim)
    if label != None: pl.legend()


insttxts = np.array(['LL', 'RR', 'LR', 'RL'])

# timing of the experiment
srate = 500                                               # sampling rate per second
Preinstr_time = int(.2 * srate)                                 # pre-instruction time (1s)
Instr_time = int(.2 * srate)                              #  instruction presentation (200 ms)
Prep_time = (np.arange(1.7,2.2,.05) * srate).astype(int)  # ISI ranging from 1700 to 2200 ms, multiplying it by the sampling rate to get how many samples we have for each ISI
Stim_time = int(.05 * srate)                              # Stimulus presentation of 50 ms
Resp_time = .7 * srate                                     # max response time of 1s
FB_time = int(.1 * srate)                                 # Feedback presentation of 500 ms
ITI = (np.arange(0,.1,.025) * srate).astype(int)          # ITI ranging from 1000 ms to 1900 ms
Response_deadline = .7 * srate                             # Response deadline

# max trial time
TotT = (Preinstr_time + Instr_time + max(Prep_time) + Stim_time + Resp_time + FB_time + max(ITI)).astype(int)  


x = np.arange(Phase.shape[-2])/srate - .2

fig, axs = pl.subplots(4, 1)
for ind, i in enumerate(np.arange(4)):
    # task relevant units
    axs[ind].plot(x, Phase[LFC_sync[Design[i,0], :].astype(np.int), 0, :, i].T)#, color=['red', 'orange', 'purple', 'cyan'])
    #axs[ind].plot(x, Phase[LFC_sync[Design[i,0], 0].astype(np.int), 0, :, i])#, color=['red', 'orange', 'purple', 'cyan'])
    axs[ind].plot(x, Integr[:, :, i].T)
    axs[ind].plot(x, MFC[0, :, i].T, 'r--')
    # task irrelevant units
    # task_irr_units_mask = np.ones(8, dtype=np.bool)
    # task_irr_units_mask[LFC_sync[Design[i,0], :].astype(np.int)] = False
    # axs[i].plot(x, Phase[task_irr_units_mask, 0, :, i].T, color=[0, 0, .5])
    axs[ind].axvline(x=0, color='k')
    axs[ind].axvline(x=.200, color='k')
    axs[ind].axvline(x=(Prep_time[Design[i, 2]] + Instr_time)/srate, color='k')
    axs[ind].axvline(x=(Prep_time[Design[i, 2]] + Instr_time + Resp_time)/srate, color='k')
    axs[ind].axvline(x=RT[i]/1000 + (Prep_time[Design[i, 2]] + Instr_time)/srate, linestyle='--',
        color=['r', 'g'][int(accuracy[i])])
    axs[ind].set_ylabel(insttxts[Design[i,0]])


# plot excitatory and inhibitory units

x = np.arange(Phase.shape[-2])/srate - .2

fig, axs = pl.subplots(4, 1)
for ind, i in enumerate(np.arange(4)):
    rule_relevant_Snode = LFC_sync[Design[i,0], 0].astype(np.int)
    axs[ind].plot(x, Rate[rule_relevant_Snode, :, i]*10, color='k')#, color=['red', 'orange', 'purple', 'cyan'])
    axs[ind].plot(x, MFC[0, :, i].T, 'r--')

    axs[ind].plot(x, Phase[rule_relevant_Snode, 0, :, i], color='r')#, color=['red', 'orange', 'purple', 'cyan'])
    axs[ind].plot(x, Phase[rule_relevant_Snode, 1, :, i], color='b')#, color=['red', 'orange', 'purple', 'cyan'])
    axs[ind].plot(x, Integr[:, :, i].T)
    axs[ind].axvline(x=(Prep_time[Design[i, 2]] + Instr_time)/srate, color='k')
    axs[ind].axvline(x=(Prep_time[Design[i, 2]] + Instr_time + Resp_time)/srate, color='k')
    axs[ind].axvline(x=RT[i]/1000 + (Prep_time[Design[i, 2]] + Instr_time)/srate, linestyle='--',
        color=['r', 'g'][int(accuracy[i])])
    axs[ind].set_ylabel(insttxts[Design[i,0]])


    #### Plotting Rien #########
x = np.arange(Phase.shape[-2])/srate - .2
insttxts = np.array(['LL', 'RR', 'LR', 'RL'])
fig, axs = pl.subplots()
instr = 3
# task relevant units
axs.plot(x, Phase[LFC_sync[Design[instr ,0], 0].astype(np.int), 0, :, instr ].T, label = "cortical column")#, color=['red', 'orange', 'purple', 'cyan'])
axs.plot(x, Phase[LFC_sync[Design[instr ,0], 2].astype(np.int), 0, :, instr ].T, label = "cortical column")#, color=['red', 'orange', 'purple', 'cyan'])
#axs[ind].plot(x, Phase[LFC_sync[Design[i,0], 0].astype(np.int), 0, :, i])#, color=['red', 'orange', 'purple', 'cyan'])
#axs.plot(x, Integr[:, :, i].T)
axs.plot(x, MFC[0, :, instr ].T, 'b--', label ="MFC theta")
# task irrelevant units
# task_irr_units_mask = np.ones(8, dtype=np.bool)
# task_irr_units_mask[LFC_sync[Design[i,0], :].astype(np.int)] = False
# axs[i].plot(x, Phase[task_irr_units_mask, 0, :, i].T, color=[0, 0, .5])
#axs.axvline(x=0, color='k')
#axs.axvline(x=.200, color='k')
#axs.axvline(x=(Prep_time[Design[instr , 2]] + Instr_time)/srate, color='k')
#axs.axvline(x=(Prep_time[Design[instr , 2]] + Instr_time + Resp_time)/srate, color='k')
#axs.axvline(x=RT[instr ]/1000 + (Prep_time[Design[instr , 2]] + Instr_time)/srate, linestyle='--',
#        color=['r', 'g'][int(accuracy[instr ])])
axs.set_ylabel(insttxts[Design[instr ,0]])
axs.legend() 

x = np.arange(Phase.shape[-2])/srate - .2
insttxts = np.array(['LL', 'RR', 'LR', 'RL'])
fig, axs = pl.subplots()
instr = 3
# task relevant unit
axs.plot(x, Phase[LFC_sync[Design[instr ,0], 0].astype(np.int), 0, :, instr ].T, label = "Phase code neurons of cortical column 5 (40 Hz)", color = 'dodgerblue', linewidth = 3)#, color=['red', 'orange', 'purple', 'cyan'])
axs.plot(x, Phase[LFC_sync[Design[instr ,0], 2].astype(np.int), 0, :, instr ].T, label = "Phase code neurons of cortical column 1 (38 Hz)", color = 'lightpink', linewidth = 3)#, color=['red', 'orange', 'purple', 'cyan'])
axs.plot(x, MFC[0, :, instr ].T, '--', label ="MFC theta (4 Hz)",  linewidth = 5, color = 'gold')
axs.set_ylabel("Neural activity (μV)", fontsize = 30); axs.set_xlabel("Time (s)", fontsize = 25)
pl.setp(axs.get_xticklabels(), fontsize=25)
pl.setp(axs.get_yticklabels(), fontsize=25)
axs.legend(fontsize = 30, loc = "upper right") 



# Plot task-relevant units
x = np.arange(Phase.shape[-2])/srate - .2
insttxts = np.array(['LL', 'RR', 'LR', 'RL'])
fig, axs = pl.subplots()
instr = 3
# task relevant unit
axs.plot(x, Phase[LFC_sync[Design[instr ,0], 2].astype(np.int), 0, :, instr ].T, label = "", color = 'dodgerblue', linewidth = 3)#, color=['red', 'orange', 'purple', 'cyan'])
axs.plot(x, Phase[LFC_sync[Design[instr ,0], 0].astype(np.int), 0, :, instr ].T, label = "", color = 'lightpink', linewidth = 3)#, color=['red', 'orange', 'purple', 'cyan'])
axs.plot(x, MFC[0, :, instr ].T, '--', label ="",  linewidth = 5, color = 'gold')
axs.set_ylabel("Neural activity (μV)", fontsize = 30); axs.set_xlabel("Time (s)", fontsize = 30)
axs.plot(x, Integr[0, :, instr].T, linewidth = 5, label = "Left press")
axs.plot(x, Integr[2, :, instr].T, linewidth = 5, label = "Right press")
#axs.set_title("Response threshold = 5", fontsize = 30)
pl.setp(axs.get_xticklabels(), fontsize=25)
pl.setp(axs.get_yticklabels(), fontsize=25)
#axs.axvline(x=(Prep_time[Design[instr , 2]] + Instr_time)/srate, color='k')
#axs.axvline(x=RT[instr ]/1000 + (Prep_time[Design[instr , 2]] + Instr_time)/srate, linestyle='--', color=['r', 'g'][int(accuracy[instr ])])
axs.legend(fontsize = 25, loc = "upper right") 