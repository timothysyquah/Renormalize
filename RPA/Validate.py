#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:07:16 2022

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
from RPA_functions import *
from scipy.interpolate import CubicSpline

plt.close('all')
plt.figure()

data = np.loadtxt('QinMorse.dat')
data_sort = np.argsort(data[:,0])

plt.plot(data[data_sort,0],data[data_sort,1],'ok',label = 'Morse-Qin-2011')


N = 10
b = 1
x = np.linspace(0,4,1000)
chi = np.array([0.,6.0,8.0,10.0])/N
# chi = np.round(np.arange(10,10.51,0.1),2)
for i in range(len(chi)):

    isk = inverse_Sk_diblock(x*x*N/6,N,chi[i]) # units are in Rg
    C = np.power(N*b/6,3/2)/N
    plt.plot(x*np.sqrt(N/6),1/isk/N,label = f'$\chi N= {N*chi[i]}$')

# plt.ylim(0,0.05)
plt.xlim(0,12)
plt.legend()
plt.xlabel(r'$q$ $R_{g,0}$')
plt.ylabel('$S(q)$/$CN$')

x = np.linspace(0.1,4,1000)

minx = np.min(x)
maxx = np.max(x)
chilist = np.linspace(0,10.5,1000)
plt.figure()
Speaklist = []
for i in range(len(chilist)):
    chi = chilist[i]/N
    isk = inverse_Sk_diblock(x*x*N/6,N,chi) # units are in Rg
    qRg= x*np.sqrt(N/6)
    sk = 1/isk/N
    cs = CubicSpline(qRg, sk)
    loc = cs.derivative().roots()[-1]
    Speaklist.append(cs(loc))
    
    
    # plt.plot(x,cs(x))
    # plt.plot(loc,cs(loc),'ok')
Speak = np.vstack(Speaklist)
plt.plot(chilist,1/Speak/2,'--k')
plt.xlim(0,12)
plt.ylim(0,12)

plt.figure()
Nbar = 1000
Speaklist = []
Speak_ROL_list = []
for i in range(len(chilist)):
    chi = chilist[i]/N
    isk = inverse_Sk_diblock(x*x*N/6,N,chi) # units are in Rg
    qRg= x*np.sqrt(N/6)
    sk = 1/isk/N
    cs = CubicSpline(qRg, sk)
    loc = cs.derivative().roots()[-1]
    Speaklist.append(cs(loc))
    extra = inverse_SK_diblock_ROL(0.5,N,chi,Nbar)
    Speak_ROL_list.append(np.array([chilist[i]*N,extra,1/cs(loc)]))
    # plt.plot(x,cs(x))
    # plt.plot(loc,cs(loc),'ok')
Speak = np.vstack(Speaklist).reshape((len(Speaklist)))
Speak_ROL = np.vstack(Speak_ROL_list)
# plt.plot(chilist,Speak/2,'--k')
plt.plot(chilist,Speak_ROL[:,1]*np.sqrt(Nbar),'--r')
plt.xlim(0,10.5)
plt.ylim(-100,1000)

plt.close('all')

plt.figure()


plt.plot(chilist+Speak_ROL[:,1]/2,Speak_ROL[:,1]/2,'--r')
maping = CubicSpline(chilist[:-1]+Speak_ROL[:-1,1]/2,Speak_ROL[:-1,1])
Sinv_Rnorm = maping(chilist[:-1])+Speak_ROL[:-1,2]

plt.xlabel('$\chi_e$')
plt.ylabel('$cN\delta S^{-1}(q^* R_g)$')
plt.xlim(0,10.5)
plt.ylim(-2,2)

plt.figure()


plt.plot(chilist[:-1],Sinv_Rnorm/2,'--r')
plt.plot(chilist,1/Speak/2,'--k')


dataset = np.vstack((np.array(chilist),1/Speak/2))

np.savetxt('RPA_Sinv.dat',dataset)
