#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 18:00:21 2022

@author: tquah
"""

import pickle 
import numpy as np
import matplotlib.pyplot as plt


file = open("dataC1.pkl",'rb')

dictionary = pickle.load(file)
header = list(dictionary)

#N=100 chi=0.0

plt.close('all')
S = dictionary['N64.0_chi0.0']['SKAA']

plt.plot(S[:,0]*4,S[:,1])
plt.ylim(0,3)
# plt.yscale('log')