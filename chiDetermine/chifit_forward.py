#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 15:38:49 2022

@author: tquah
"""
import sys
sys.path.append('../RPA/')
from RPA_functions import *
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,least_squares
from scipy.interpolate import interp1d,CubicSpline
import numpy as np
import re
import os

def Pade_Approximation(_x,_a,_b,_c):
    _num = _a*_x+_b*np.square(_x)
    _den = 1+_c*_x
    return _num/_den

class CL_Dataset():
    def __init__(self,_chib,_Sinv2,_N,_dict,_f,_c):
        self._chib = _chib/_N
        self._Sinv2 = _Sinv2
        self._N = _N
        self._dict = _dict
        self._f = _f
        self._c = _c

    def CostFunction(self,_x):
        #we want to tweak x until S and SROL is minimized
        #x should be our parameters that map chib to chie
        #may need to add constraint on chi
        _chie = Pade_Approximation(self._chib,_x[0],_x[1],_x[2])
        _chieN = _chie*self._N
        _chic, _scalc,loc = self._dict[(self._f,self._c,self._N)].Evaluate_Chie(_chieN)
        print(self._Sinv2[loc])
        print(np.abs(_scalc-self._Sinv2[loc]) )
        return np.abs(_scalc-self._Sinv2[loc]) 
    
    
    



plt.close('all')



ROL_Interp = Build_ROL('../RPA/ROL_dataset.dat')
header = list(ROL_Interp)

#validate mathods
#blind points
ROLREF = np.loadtxt('../RPA/ROL_dataset.dat')
plt.plot(ROLREF[:,-2],ROLREF[:,-1]/2,'ok',markersize=1)

plt.xlim(-1,15)
plt.ylim(0,12)

for i in range(len(header)):
    minmax,sminmax = ROL_Interp[header[i]].MinMax()
    chie = np.linspace(minmax[0],minmax[1],1000)
    chiea,Sinva,loc = ROL_Interp[header[i]].Evaluate_Chie(chie)
    plt.plot(chiea,Sinva/2)

# Ok First need to write in the CL dataset

CL_Data = np.loadtxt("../PeakExtraction/CL_dataset_large.dat")
Ncutoff = 100
loc = np.where(CL_Data[:,2]>Ncutoff)[0]
CL_Data_Safe = CL_Data[loc,:]
uniqueN = np.unique(CL_Data_Safe[:,2])
f = 0.5 # need to supply this
c = 1.0 # need to supply this

# need to remove anydata we cannot probe using emperical ROL

# Newlist = []
# for i in range(0,len(uniqueN)):
#     if (f,c,uniqueN[i]) not in header:
#         print(f"Missing f={f} c={c} N={uniqueN[i]}")
#         continue
#     loc = header.index((f,c,uniqueN[i]))
#     chiBounds, Sbounds = ROL_Interp[header[loc]].MinMax()
    
#     #remove anydata using mask
    
#     loc_cl = np.where(CL_Data_Safe[:,2]==uniqueN[i])[0]
#     chib = CL_Data_Safe[loc_cl,0]
#     S = CL_Data_Safe[loc_cl,4]
#     Sinv2 = uniqueN[i]/S/2
    
#     #apply masks
    
#     m1 = np.ma.masked_where(chib>=chiBounds[0],chib).mask # not sure if chib mask makes sense need to test
#     m2 = np.ma.masked_where(chib<=chiBounds[1],chib).mask
#     m3 = np.ma.masked_where(Sinv2>=Sbounds[0],Sinv2).mask
#     m4 = np.ma.masked_where(Sinv2<=Sbounds[1],Sinv2).mask
#     m = m1*m2*m3*m4
#     finloc = np.where(m==True)[0]
    
#     Newlist.append(np.array([uniqueN[i]*np.ones(len(finloc)),chib[finloc],Sinv2[finloc]]))
# CleanData = np.hstack(Newlist).transpose()
saveparam = []
for i in range(0,len(uniqueN)):
    if (f,c,uniqueN[i]) not in header:
        print(f"Missing f={f} c={c} N={uniqueN[i]}")
        continue
    loc = header.index((f,c,uniqueN[i]))
    # loc_cl = np.where(CleanData[:,0]==uniqueN[i])[0]
    loc_cl = np.where(CL_Data_Safe[:,2]==uniqueN[i])[0]
    # Sinv2 = uniqueN[i]/S/2
       
    fitinit = CL_Dataset(CL_Data_Safe[loc_cl,0],uniqueN[i]/CL_Data_Safe[loc_cl,4]/2,uniqueN[i],ROL_Interp,f,c)
    param = least_squares(fitinit.CostFunction,x0 = np.ones(3)*0.5,xtol = 1e-12,ftol = 1e-12)
    print(param)
    saveparam.append(param.x)
    chie = Pade_Approximation(CL_Data_Safe[loc_cl,0],param.x[0],param.x[1],param.x[2])
    plt.plot(chie,uniqueN[i]/CL_Data_Safe[loc_cl,4]/2,'ok')