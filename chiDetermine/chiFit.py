#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 13:33:18 2022

@author: tquah
"""
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d,CubicSpline
import numpy as np
import re
import os


def ExtractNumber(line):
    return re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)

def Pade_Approximation(x,a,b,c,d):
    num = a*x+b*x**2
    den = 1+c*x+d*x**2
    return num/den

Nbarfiles = glob.glob('ROL*')

Renorm = dict()
CL_data = np.loadtxt('CL_dataset.dat')

for i in range(0,len(Nbarfiles)):

    N = float(ExtractNumber(Nbarfiles[i])[0])
    Renorm[i] = np.loadtxt(Nbarfiles[i]).transpose()

    if N<1000:
        continue
    temp=np.loadtxt(Nbarfiles[i]).transpose()
    tempfunc = interp1d(temp[:,1],temp[:,0])
    
    row = np.where(CL_data[:,2]==N)[0]
    
    Stemp = CL_data[row,2]/CL_data[row,4]/2
    chi_e = tempfunc(Stemp)
    
    plt.figure(figsize = (6,6))
    plt.title('$N=1000$')
    plt.plot(CL_data[row,0],chi_e/CL_data[row,2],'ok',label ='Smeared')
    plt.plot(CL_data[row,1],chi_e/CL_data[row,2],'^g',label = 'Unsmeared-FM')
    plt.ylabel('$\chi_e N$')    
    plt.xlabel('$\chi_b N$')
    xdata = np.linspace(0,np.max(CL_data[row,0]),100)
    param1,pcov1 = curve_fit(Pade_Approximation,CL_data[row,0],chi_e/CL_data[row,2])
    param2,pcov2 = curve_fit(Pade_Approximation,CL_data[row,2],chi_e/CL_data[row,2])
    plt.plot(xdata,Pade_Approximation(xdata,param1[0],param1[1],param1[2],param1[3]),label='Pade-Fit')
    plt.legend()
    plt.savefig('/home/tquah/Figures/N1000Fit.pdf', dpi = 300)
    
plt.figure(figsize = (6,6))

uniqueN = np.unique(CL_data[:,2])

for i in range(0,len(uniqueN)):
    row = np.where(CL_data[:,2]==uniqueN[i])[0]
    chie = Pade_Approximation(CL_data[row,0],param1[0],param1[1],param1[2],param1[3])
    plt.scatter(chie*CL_data[row,2],CL_data[row,2]/CL_data[row,4]/2)
#     x = p    