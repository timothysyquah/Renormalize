#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:22:33 2022

@author: tquah
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import json
import pickle
from scipy.stats import f,beta
import re
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
sys.path.append('../RPA')
from RPA_functions import *
from scipy.signal import find_peaks
plt.close('all')
hackdistance = 1000
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

def Gaussian(x,a,b,c):
    return a*np.exp(-(np.square(x-b)/(2*np.square(c))))
def fdist(x,a,b,c,d):
    return c/beta.pdf(x,a,b,loc = d)

#general path 

#change directory
cutoff = 1
Working_Directory = os.getcwd()
file = open("dataC1.pkl",'rb')

dictionary = pickle.load(file)
header = list(dictionary)
datalist = []       
plt.close('all')
plt.figure()
datalist = []
for i in range(0,len(header),1):
    # i=3
    N = dictionary[header[i]]['N']
    chi_a = dictionary[header[i]]['Chi']

    S = (dictionary[header[i]]['SKAA'][:,1]+ dictionary[header[i]]['SKAA'][:,1]-2* dictionary[header[i]]['SKAB'][:,1])
    SAA = dictionary[header[i]]['SKAA'][:,1]#/np.power(N,1.5)
    x = dictionary[header[i]]['SKAA'][:,0]*np.sqrt(N)
    xmainplot = np.logspace(-10,2,1000)

    chi = 2/np.mean(dictionary[header[i]]['Chi_I'][-1000:,1]) #set this by default in future should use Kris' Stats stuff to get the actual one
    print(chi*N)
    # print(chi/np.sqrt(N))
    # print(chi_a)
    loc = np.where(x/np.sqrt(N)<cutoff)[0]
    Speak = np.max(S[loc]) #this is a naiive way to pick out the peak
    peakloc = np.where(Speak==S)[0]
    plt.plot(x[loc]/np.sqrt(N),S[loc]/np.power(N,1.5),label = rf'$N = {N}$, $\chi = {chi_a*N}$')
    plt.plot(x[peakloc]/np.sqrt(N),Speak/np.power(N,1.5),'ok')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    #first remove noisy section
    plt.ylim(0.01,10.0)
    plt.yscale('log')

    plt.tight_layout()
    datalist.append(np.array([chi_a,chi,N,x[peakloc][0],Speak]))
    
data = np.vstack(datalist)



np.savetxt('CL_dataset_large.dat',data) 
uniqueN = np.unique(data[:,2])
plt.figure()
for i in range(0,len(uniqueN)):
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    #first remove noisy section
    plt.ylim(0.01,10.0)
    plt.yscale('log')

    plt.tight_layout()
    datalist.append(np.array([chi_a,chi,N,x[peakloc][0],Speak]))
    
data = np.vstack(datalist)



np.savetxt('CL_dataset_large.dat',data) 
uniqueN = np.unique(data[:,2])
plt.figure()
for i in range(0,len(uniqueN)):
    loc = np.where(data[:,2]==uniqueN[i])[0]
    plt.scatter(data[loc,2]*data[loc,0],(data[loc,2])/data[loc,4]/2,label = fr'$N = {uniqueN[i]}$')
plt.legend()    
plt.xlabel(r'$\alpha N$')
plt.ylabel(r'$CNS^{-1}(q)/2$')
plt.savefig('/home/tquah/Figures/invS.pdf',dpi = 300)


# plt.figure()
# for i in range(0,len(uniqueN)):
#     loc = np.where(data[:,2]==uniqueN[i])[0]
#     plt.scatter(data[loc,2]*data[loc,0],data[loc,4]/np.square(data[loc,2])*np.sqrt(data[loc,2]),label = fr'$N = {uniqueN[i]}$')
# plt.legend()    

    # plt.yscale('log')