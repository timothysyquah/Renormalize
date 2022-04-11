#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:55:46 2022

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f

# plt.close('all')
plt.figure()

dfn, dfd = 29, 18
mean, var, skew, kurt = f.stats(dfn, dfd, moments='mvsk')   
x = np.linspace(0,30,100)


dfn = 10    
dfd = 10
c =2
scale = 10
plt.plot(x, (f.pdf(x, dfn, dfd,loc = 2)-c)*scale,
       'r-', lw=5, alpha=0.6, label='f pdf')