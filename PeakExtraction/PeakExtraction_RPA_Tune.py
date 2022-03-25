#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:07:40 2022

@author: tquah
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import json
import pickle
import re
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
sys.path.append('../RPA')
from RPA_functions import *
from scipy.signal import find_peaks
plt.close('all')
hackdistance = 860
#Function to get all directories
def ExtractNumber(line):
    return re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)
def Prune(Dlist,Desired_Depth):
    returndire = []
    for dire in Dlist:
        split = dire.split('/')
        if Desired_Depth+1==len(split):
            returndire.append(dire)
    return returndire

#general path 

#change directory

Working_Directory = os.getcwd()
file = open("data.pkl",'rb')

dictionary = pickle.load(file)
header = list(dictionary)
datalist = []
for i in range(0,len(header),1):
    # i=3
    N = dictionary[header[i]]['N']
    chi_a = dictionary[header[i]]['Chi']

    S = (dictionary[header[i]]['SKAA'][:,1]+ dictionary[header[i]]['SKAA'][:,1]-2* dictionary[header[i]]['SKAB'][:,1])##/np.power(N,1.5)
    SAA = dictionary[header[i]]['SKAA'][:,1]#/np.power(N,1.5)
    x = dictionary[header[i]]['SKAA'][:,0]*np.sqrt(N)
    xmainplot = np.logspace(-10,2,1000)

    chi = np.mean(dictionary[header[i]]['Chi_I'][-1000:,1]) #set this by default in future should use Kris' Stats stuff to get the actual one
    
    isk_1 = inverse_Sk_diblock(x*x,N,chi) # units are in Rg
    isk_2 = inverse_Sk_diblock(x*x,N,chi_a) # units are in Rg
    F = Fitting_Functions(0.5,N)
    loc = np.where(x/np.sqrt(N)<1.2)[0]
    
    param,pcov = curve_fit(F.SK_Diblock, x[loc], S[loc])
    param2,pcov2 = curve_fit(F.SK_Diblock_Add, x[loc], S[loc],p0 = [chi/N,0,0,0,0])

    Sfit = F.SK_Diblock(x,param)
    Sfit_add = F.SK_Diblock_Add(x, param2[0],param2[1], param2[2], param2[3],param2[4])
    Sfit_add_plot = F.SK_Diblock_Add(xmainplot, param2[0],param2[1], param2[2], param2[3],param2[4])
    logic =np.isnan(Sfit_add_plot)
    nonnan = np.where(logic==0)
    cs = CubicSpline(xmainplot[nonnan],Sfit_add_plot[nonnan])
    locus = cs.derivative().roots()
    peak_loc = find_peaks(Sfit_add_plot,prominence=0.5, width=10)[0]
    distance = np.abs(peak_loc-hackdistance)
    pick_peak = np.where(distance ==np.min(distance))
    dictionary[header[i]]['PeakStuff'] = np.array([chi_a,N, xmainplot[peak_loc[pick_peak]][0],\
                                                   Sfit_add_plot[peak_loc[pick_peak]][0]])
    datalist.append(dictionary[header[i]]['PeakStuff'] )
    # if len(peak_loc)==0:
    # plt.figure()

    # plt.plot(x,S,'r^',label = r'$S(q)/CN$')
    # plt.plot(x,SAA,'ok',label = r'$S_{AA}(q)/CN$')
    
    # plt.ylabel('$S(q)/CN$')
    # plt.xlabel('$q \ (R^{-1}_{g,0})$')
    # # plt.plot(x,1/isk_2,label =fr'$\chi N = {chi_a*N}$')
    # # plt.plot(x,1/isk_1,label =fr'$\chi N = {chi*N:0.2f}$')
    # # plt.plot(x,1/isk_1,label =fr'$\chi N = {chi*N:0.2f}$')
    # plt.plot(x,Sfit,label =fr'$\chi N = {param[0]*N:0.2f}-RPA-Fit$')
    # plt.plot(x,Sfit_add,label =fr'$\chi N = {param2[0]*N:0.2f}-RPA-AugFit$')
    # plt.plot(xmainplot,Sfit_add_plot,label =fr'$\chi N = {param2[0]*N:0.2f}-RPA-AugFit$')
    # plt.plot(xmainplot,cs(xmainplot),label =fr'$\chi N = {param2[0]*N:0.2f}-RPA-AugFit$')
    # # plt.plot(locus,cs(locus),'sb')
    # plt.plot(xmainplot[peak_loc],Sfit_add_plot[peak_loc],'Dg')

    # plt.xlim(0,np.max(x))
    # plt.ylim(np.min(S[loc]),np.max(S[loc]))
    # plt.legend()
    # break
data = np.vstack(datalist)
data[:,3] = data[:,3]/data[:,1]
Nunique = np.unique(data[:,1])
plt.figure()
for i in range(0,len(Nunique)):
    loc = np.where(data[:,1]==Nunique[i])[0]
    plt.scatter(data[loc,0]*data[loc,1],1/data[loc,3]/2,label = fr"$N={Nunique[i]}$")
    
    
rpa = np.loadtxt('../RPA/RPA_Sinv.dat').transpose()
plt.plot(rpa[:,0],rpa[:,1])
plt.legend()