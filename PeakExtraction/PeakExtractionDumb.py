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
file = open("dataC1_zeta10.pkl",'rb')

dictionary = pickle.load(file)
header = list(dictionary)
datalist = []       
plt.close('all')
plt.figure()
datalist = []
for i in range(0,len(header),1):
    N = dictionary[header[i]]['N']
    chi_a = dictionary[header[i]]['Chi']
    if chi_a*N < 2.0:
        continue
    
    if N==256:
        continue
    if N==100:
        continue
    S = (dictionary[header[i]]['SKAA'][:,1]+ dictionary[header[i]]['SKBB'][:,1]-2* dictionary[header[i]]['SKAB'][:,1])
    x = dictionary[header[i]]['SKAA'][:,0]*np.sqrt(N)
    xmainplot = np.logspace(-10,2,1000)

    chi = 2/np.mean(dictionary[header[i]]['Chi_I'][-1000:,1]) #set this by default in future should use Kris' Stats stuff to get the actual one
    loc = np.where(x/np.sqrt(N)<cutoff)[0]
    Speak = np.max(S[loc]) #this is a naiive way to pick out the peak
    peakloc = np.where(Speak==S)[0]
    datalist.append(np.array([chi_a*N,chi*N,N,x[peakloc][0],Speak/4]))    
data = np.vstack(datalist)
np.savetxt('CL_dataset_large.dat',data) 
