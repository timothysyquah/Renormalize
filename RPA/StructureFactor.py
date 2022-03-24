#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:07:16 2022

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
from RPA_functions import *


plt.close('all')
plt.figure()

data = np.loadtxt('QinMorse.dat')
data_sort = np.argsort(data[:,0])

plt.plot(data[data_sort,0],data[data_sort,1],'ok',label = 'Morse-Qin-2011')


N = 64
b = 1
x = np.linspace(0,4,1000)
chi = np.array([0.])
# chi = np.round(np.arange(10,10.51,0.1),2)
for i in range(len(chi)):

    isk = inverse_Sk_diblock(x*x*N/6,N,chi[i]) # units are in Rg
    C = np.power(N*b/6,3/2)/N
    plt.plot(x*np.sqrt(N/6),1/isk/N,label = f'$\chi= {chi[i]}$')

plt.ylim(0,0.05)
plt.xlim(0,12)
plt.legend()
plt.xlabel(r'$q$ $R_{g,0}$')
plt.ylabel('$S(q)$/$CN$')