import numpy as np
import matplotlib.pyplot as pl
from scipy import signal as sig
from mne import parallel as par


def temporal_surrogates(x, n_surrs, n_t):
	temp_surr = np.zeros(shape = [n_surrs, n_t])
	rand_ind = np.arange(n_t)
	for surr_n in range(n_surrs):
		np.random.shuffle(rand_ind)
		temp_surr[surr_n, :] = x[rand_ind]

	return np.abs(np.fft.fft(temp_surr, axis = 1)[:, 1:int((n_t/2 + 1))])

################################################################################################################################################
################################################################################################################################################

# Simulate one sinusoid and our sampling procedure

n_time_points = 1500

# simulate a sinusoid
freq = 6
start_t = 0
end_t = 1.5
x = np.linspace(start_t, end_t, n_time_points)
sine = np.sin(x*freq*np.pi*2)

# add noise
avg_noise = 1
std_noise = .5
noise = (np.random.randn(n_time_points) + avg_noise) * std_noise
sine_noise = sine + noise

# sample from it as we do in the instructWM protocol
n_samples = 11
start_t_sample = .5
end_t_sample = 1
samples_t = np.linspace(start_t_sample, end_t_sample, n_samples)
samples_inds = np.array([np.where(x <= st)[0][-1]+1 for st in samples_t])
sine_noise_sampled = sine_noise[samples_inds]

# compute frequency spectrum of sampled data
fft_sampled = np.fft.fft(sine_noise_sampled)
amp_spectrum_sampled = np.abs(fft_sampled)
sampling_interval = (end_t_sample - start_t_sample)/(n_samples-1)
sampling_freq = 1/sampling_interval
freqs_sampled = np.fft.fftfreq(n_samples, sampling_interval)
mask_good_freqs = np.arange(1, int(n_samples/2. + 1))
freqs_sampled_good = np.abs(freqs_sampled[mask_good_freqs])

# plot

fig, axs = pl.subplots(4,1)
axs[0].plot(x, sine)
axs[1].plot(x, sine_noise)
axs[2].plot(samples_t, sine_noise_sampled)
axs[2].set_xlim(start_t, end_t)

axs[3].plot(freqs_sampled_good, amp_spectrum_sampled[mask_good_freqs], 'o-')



################################################################################################################################################
################################################################################################################################################

###### Simulate for a cohort of observers ######
n_obs = 50

n_time_points = 1500

# simulate a sinusoid
freq = 6
start_t = 0
end_t = 1.5
x = np.linspace(start_t, end_t, n_time_points)

# noise params
avg_noise = 0
std_noise = .25

# sample params
n_samples = 11
start_t_sample = .5
end_t_sample = 1
samples_t = np.linspace(start_t_sample, end_t_sample, n_samples)
samples_inds = np.array([np.where(x <= st)[0][-1]+1 for st in samples_t])

sine_all = np.zeros(shape = [n_obs, n_time_points])
sine_noise_all = np.zeros(shape = [n_obs, n_time_points])
sampled_all_obs = np.zeros(shape = [n_obs, n_samples])
amp_spectrum_sampled_all = np.zeros(shape = [n_obs, len(freqs_sampled_good)])

# individual variability parameters
phase_var = np.random.rand(n_obs) * np.pi * 2
phase_var = np.linspace(0, 2*np.pi, n_obs)

for obs_i in np.arange(n_obs):
	sine_all[obs_i, :] = np.sin(phase_var[obs_i] + (x*freq*np.pi*2))/4 + .75

	# add noise
	noise = (np.random.randn(n_time_points) + avg_noise) * std_noise
	sine_noise_all[obs_i, :] = sine_all[obs_i, :] + noise

	sampled_all_obs[obs_i, :] = sine_noise_all[obs_i, samples_inds]

	# compute spectrum
	fft_sampled = np.fft.fft(sampled_all_obs[obs_i, :])
	amp_spectrum_sampled = np.abs(fft_sampled)
	amp_spectrum_sampled_all[obs_i, :] = amp_spectrum_sampled[mask_good_freqs]



fig, axs = pl.subplots(4,1)
fig.set_figheight(8)
# Underlying pure oscillation
for obs_i in np.arange(n_obs):
	axs[0].plot(x, sine_all[obs_i, :])
axs[0].plot(x, sine_all.mean(axis = 0), 'k', linewidth=2)
axs[0].set_title('Hypothetical underlying oscillation', fontsize=8)
axs[0].set_ylabel('Performance')

# Noisy
axs[1].plot(x, sine_noise_all.mean(axis = 0), 'k')
axs[1].set_title('Adding NOISE', fontsize=8)
axs[1].set_ylabel('Performance', fontsize=8)

# Sampled time points
for obs_i in np.arange(n_obs):
	axs[2].plot(samples_t, sampled_all_obs[obs_i, :], linewidth=.5)
axs[2].plot(samples_t, sampled_all_obs.mean(axis = 0), 'k', linewidth=2)
axs[2].set_xlim(start_t_sample, end_t_sample)
axs[2].set_title('Noisy oscillation sampled with %s time points (at %.1f Hz)' %\
	(n_samples, sampling_freq), fontsize=8)
axs[2].set_xlabel('Time (s)', fontsize=8)
axs[2].set_ylabel('Performance', fontsize=8)

# Frequency spectrum
axs[3].errorbar(freqs_sampled_good, amp_spectrum_sampled_all.mean(axis=0),
	amp_spectrum_sampled_all.std(axis=0)/np.sqrt(n_obs))
axs[3].set_xlabel('Frequency (Hz)', fontsize=8)
axs[3].set_ylabel('Amplitude', fontsize=8)

# pl.subplots_adjust(hspace=.5)

##########################################
# compute surrogates
n_surrs = 100000
# compute surrogates
parallel, my_cvstime, _ = par.parallel_func(temporal_surrogates, n_jobs = -1, verbose = 40)
surr_sampled = np.rollaxis(np.array(parallel(my_cvstime(x = x, n_surrs = n_surrs, n_t = n_samples) for x in sampled_all_obs)), 1)

# plot results + surrogate
ci = 1 - .05 / len(freqs_sampled_good)
surr_sampled_avg = surr_sampled.mean(axis = 1)
allftimefunction = np.sort(surr_sampled_avg, 0)
percentile = np.int(np.floor(ci * n_surrs))
upperlim = allftimefunction[percentile, :]
expected = surr_sampled_avg.mean(axis = 0)

acc_spectra = amp_spectrum_sampled_all.mean(axis=0)
pvals = np.zeros(len(freqs_sampled_good))
for ind, sp in enumerate(acc_spectra):
	temp = np.where(sp < allftimefunction[:, ind])[0]
	if temp.size:
		pvals[ind] = 1 - temp[0] / float(n_surrs)

print('pvalues for %s (Hz):\n\t%s\n' % (np.str(freqs_sampled_good), np.str(pvals)))
#####################


# plot surrogates
axs[3].plot(freqs_sampled_good, upperlim, 'k--')
axs[3].plot(freqs_sampled_good, expected, 'k.-')


pl.tight_layout()
pl.subplots_adjust(hspace=.5)














