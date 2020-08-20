# -*- coding: utf-8 -*-

"""
Author: Rien Sonck
Date: 2020-08-08
Description: This file plots a figure demostrating parameter space
"""
######################
# Importing packages #
######################
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl
import random
import numpy as np

############
# Plotting #
############

# create new figure
figure = plt.figure()
# create third axis
axis = mpl.axes3d.Axes3D(figure)
# 3 lists of coordinates
x_coordinates = np.random.rand(100) * 100
y_coordinates = np.random.rand(100) * 100
random.shuffle(x_coordinates)
z_coordinates = (x_coordinates  ** 2) + np.random.rand(100) * 100
# shuffle coordinates
axis.scatter(x_coordinates, y_coordinates, z_coordinates, s = 50)
#axis.set_xlabel("Parameter 1", fontsize = 18); axis.set_ylabel("Parameter 2", fontsize = 18); axis.set_zlabel("Behavioral outcome", fontsize = 18)
axis.set_yticklabels([])
axis.set_xticklabels([])
axis.set_zticklabels([])