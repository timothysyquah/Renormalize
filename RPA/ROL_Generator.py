#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 10:50:50 2022

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt
from RPA_functions import *
from scipy.interpolate import CubicSpline

#Sweep over f, c, chi, N
fsweep = np.linspace(0.5,0.5,1)
csweep = np.array([1.0,0.1])
chisweep = np.linspace(0,10.5,100)
Nsweep = np.array([16,32,64,128,256])

#array ordering
sweeparray = [fsweep,csweep,chisweep,Nsweep]
label = ['f','c','chiN','N']
#use methods in RPA Functions
dictionary = generator_ROL_emperical(sweeparray,label)
array,header = convert_dictionary_to_npy(dictionary)
np.savetxt("ROL_dataset.dat",array,delimiter=" ",header = header)

