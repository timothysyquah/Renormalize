#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 13:33:18 2022

@author: tquah
"""
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,leastsq
from scipy.interpolate import interp1d,CubicSpline
import numpy as np
import re
import os

plt.close('all')
def ExtractNumber(line):
    return re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)

def Pade_Approximation(x,a,b,c):
    num = a*x+b*x**2
    den = 1+c*x
    return num/den


# def Cost_Function(params,x,y):
    
    
    
    


Nbarfiles = glob.glob('ROL*dat')

Renorm = dict()
CL_data = np.loadtxt('../PeakExtraction/CL_dataset_large.dat')


fig1 = plt.figure()
ax1 = fig1.gca()
fig2 = plt.figure()
ax2 = fig2.gca()
logic = True

chimatchlist = []


for i in range(0,len(Nbarfiles)):

    N = float(ExtractNumber(Nbarfiles[i])[0])
    if N not in CL_data[:,2]:
        continue
    
    
    Renorm[i] = np.loadtxt(Nbarfiles[i]).transpose()

    temp=np.loadtxt(Nbarfiles[i]).transpose()
    tempfunc = interp1d(temp[:,1],temp[:,0])
    minInterp = np.min(temp[:,1])
    maxInterp = np.max(temp[:,1])
    
    
    ax1.plot(temp[:,0],temp[:,1],label = f'N = {N}')
    row = np.where(CL_data[:,2]==N)[0]
    chib = CL_data[row,0]
    fmchib = CL_data[row,1] 

    Stemp = CL_data[row,2]/CL_data[row,4]/2
    
    
    
    #remove some stuff
    loca  = np.where((Stemp>minInterp)&(Stemp<maxInterp))[0]
    # print(Stemp)
    # print(minInterp)
    # print(maxInterp)
    Stemp = Stemp[loca]
    chi_e = tempfunc(Stemp)
    
    ax1.plot(chi_e,Stemp,'ok')
    ax2.scatter(chib[loca]/N,chi_e/N,label = f'N = {N}')
    
    # param,pcov = curve_fit(Pade_Approximation,chib[loca]/N,chi_e/N)
    # ax1.plot(Pade_Approximation(chib[loca]/N,param[0],param[1],param[2])*N,Stemp,'^b')
    chimatchlist.append(np.vstack((chib[loca],chi_e,Stemp,N*np.ones(len(chi_e)))).transpose())
ax1.legend()
ax2.legend()


chimatch = np.vstack(chimatchlist)
param,pcov = curve_fit(Pade_Approximation,chimatch[:,0]/chimatch[:,3],chimatch[:,1]/chimatch[:,3])

ax1.plot(Pade_Approximation(chimatch[:,0]/chimatch[:,3],param[0],param[1],param[2])*chimatch[:,3],chimatch[:,2],'sr')


x = np.linspace(0,0.8,100)
ax2.plot(x,Pade_Approximation(x,param[0],param[1],param[2]),'--k')
    # ax2.scatter(fmchib[loca],chi_e)

    # plt.figure(figsize = (6,6))
    # # plt.title('$N=1000$')
    # plt.plot(CL_data[row,0],chi_e/CL_data[row,2],'ok',label ='Smeared')
    # plt.plot(CL_data[row,1],chi_e/CL_data[row,2],'^g',label = 'Unsmeared-FM')
    # plt.ylabel('$\chi_e N$')    
    # plt.xlabel('$\chi_b N$')
    # xdata = np.linspace(0,np.max(CL_data[row,0]),100)
    # param1,pcov1 = curve_fit(Pade_Approximation,CL_data[row,0],chi_e/CL_data[row,2])
    # param2,pcov2 = curve_fit(Pade_Approximation,CL_data[row,2],chi_e/CL_data[row,2])
    # plt.plot(xdata,Pade_Approximation(xdata,param1[0],param1[1],param1[2],param1[3]),label='Pade-Fit')
    # plt.legend()
    # # plt.savefig('/home/tquah/Figures/N1000Fit.pdf', dpi = 300)
    
# plt.figure(figsize = (6,6))

# uniqueN = np.unique(CL_data[:,2])

# for i in range(0,len(uniqueN)):
#     row = np.where(CL_data[:,2]==uniqueN[i])[0]
#     chie = Pade_Approximation(CL_data[row,0],param1[0],param1[1],param1[2],param1[3])
#     plt.scatter(chie*CL_data[row,2],CL_data[row,2]/CL_data[row,4]/2)
#     x = p    