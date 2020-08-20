# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-07-31
Description: This file make splots to explain neural oscillations
"""

######################
# Importing packages #
######################
import numpy as np
import pandas as pd
from scipy import ndimage
from matplotlib import pyplot as pl


x = np.linspace(-np.pi * 2, np.pi * 2, 1000)
# plot 1: synchronized and asynchronized
fig, axs = pl.subplots(1, 2)
title = ['Synchronized phase', 'Asynchronized phase']
phase = [0.05, 4]
label = ['Signal 1', 'Signal 2']
for i, ax in enumerate(axs):
	ax.set_title(title[i], fontsize = 35)
	ax.plot(x + phase[i], np.sin(x * 4), linewidth = 4, label = label[0])
	ax.plot(x, np.sin(x * 4), linewidth = 4, label = label[1])
	ax.set_xlim(-2, 2)
	ax.legend(loc='upper right', fontsize = 35)
	ax.set_ylabel("Neural activity (μV)", fontsize=35); ax.set_xlabel("Time (s)", fontsize=30)
	ax.tick_params(
	axis='both',        # changes apply to the x-axis
	which='both',      	# both major and minor ticks are affected
	bottom=False,      	# ticks along the bottom edge are off
	top=False,         	# ticks along the top edge are off
	left=False,
	labelbottom=False,	# labels along the bottom edge are off
	labelleft=False) 	# labels along the left edge are off


# plot 2: frequency bands
fig2, axs2 = pl.subplots(5, 1)
title = ['delta band oscillations', 'theta band oscillations', 'alpha band oscillations', 'beta band oscillations', 'gamma band oscillations']
freq = [2, 6, 10, 21, 65]
x = np.linspace(-np.pi, np.pi, 1000)
for i, ax in enumerate(axs2):
	noise = np.random.normal(0 , 0.2, 1000)
	ax.set_title(title[i], fontsize = 20)
	ax.plot(x, noise + np.sin(x * freq[i]), linewidth = 4)
	ax.tick_params(
	axis='both',        # changes apply to the x-axis
	which='both',      	# both major and minor ticks are affected
	bottom=False,      	# ticks along the bottom edge are off
	top=False,         	# ticks along the top edge are off
	left=False,
	labelbottom=False,	# labels along the bottom edge are off
	labelleft=False) 	# labels along the left edge are off
	ax.axis('off')
fig2.tight_layout() 

# plot 3: 
x = np.linspace(-np.pi, np.pi, 1000)
# plot 1: synchronized and asynchronized
fig, ax = pl.subplots(1, 1)
ax.plot(x, np.sin(x * 2), linewidth = 4)
ax.set_ylabel("Neural activity (μV)", fontsize=20); ax.set_xlabel("Time (s)", fontsize=20)
ax.tick_params(
axis='both',        # changes apply to the x-axis
which='both',      	# both major and minor ticks are affected
bottom=False,      	# ticks along the bottom edge are off
top=False,         	# ticks along the top edge are off
left=False,
labelbottom=False,	# labels along the bottom edge are off
labelleft=False) 	# labels along the left edge are off

# plot 3: 
x = np.linspace(-np.pi, np.pi, 1000)
# plot 1: synchronized and asynchronized
fig, ax = pl.subplots(1, 1)
ax.plot(x, np.sin(x * 6), linewidth = 4, label ="Higher Frequency (6 Hz)")
ax.plot(x, np.sin(x * 2), linewidth = 4, label ="Lower Frequency (2 Hz)")
ax.set_ylabel("Neural activity (μV)", fontsize=35); ax.set_xlabel("Time (s)", fontsize=35)
ax.legend(loc='upper right', fontsize = 30)
ax.tick_params(
axis='both',        # changes apply to the x-axis
which='both',      	# both major and minor ticks are affected
bottom=False,      	# ticks along the bottom edge are off
top=False,         	# ticks along the top edge are off
left=False,
labelbottom=False,	# labels along the bottom edge are off
labelleft=False) 	# labels along the left edge are off

# plot 2: frequency bands
fig2, ax = pl.subplots(1, 1)
title = 'neural oscillation'
freq = 4
x = np.linspace(-np.pi, np.pi, 1000)
noise = np.random.normal(0 , 0.1, 1000)
ax.set_title(title, fontsize = 20)
ax.plot(x,  np.sin(x * freq), linewidth = 4)
ax.tick_params(
axis='both',        # changes apply to the x-axis
which='both',      	# both major and minor ticks are affected
bottom=False,      	# ticks along the bottom edge are off
top=False,         	# ticks along the top edge are off
left=False,
labelbottom=False,	# labels along the bottom edge are off
labelleft=False) 	# labels along the left edge are off
fig2.tight_layout() 

ax.axis('off')


Gaussian = np.random.normal(size=[1,])
noise = np.random.normal(0.1, size = [250,])

# plot 3: 
x = np.linspace(-np.pi, np.pi, 1000)
sin_x = np.sin(x * 38)
sin_y = np.sin(x * 42)
sin_z = np.sin(x * 5)

# plot 1: synchronized and asynchronized
fig, ax = pl.subplots(1, 1)
ax.plot(x, sin_x, linewidth = 4)
ax.plot(x, sin_y, linewidth = 4)
ax.plot(x, sin_z, linewidth = 4)
ax.tick_params(
axis='both',        # changes apply to the x-axis
which='both',      	# both major and minor ticks are affected
bottom=False,      	# ticks along the bottom edge are off
top=False,         	# ticks along the top edge are off
left=False,
labelbottom=False,	# labels along the bottom edge are off
labelleft=False) 	# labels along the left edge are off
ax.axis('off')