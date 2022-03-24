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
sys.path.append('../RPA')
from RPA_functions import *

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

#general path 

#change directory

Working_Directory = os.getcwd()
file = open("data.pkl",'rb')

dictionary = pickle.load(file)
header = list(dictionary)
for i in range(0,len(header),1):
    N = dictionary[header[i]]['N']
    chi_a = dictionary[header[i]]['Chi']

    S = (dictionary[header[i]]['SKAA'][:,1]+ dictionary[header[i]]['SKAA'][:,1]-2* dictionary[header[i]]['SKAB'][:,1])##/np.power(N,1.5)
    SAA = dictionary[header[i]]['SKAA'][:,1]#/np.power(N,1.5)
    plt.figure()
    x = dictionary[header[i]]['SKAA'][:,0]*np.sqrt(N)
    chi = np.mean(dictionary[header[i]]['Chi_I'][-1000:,1]) #set this by default in future should use Kris' Stats stuff to get the actual one
    
    isk_1 = inverse_Sk_diblock(x*x,N,chi) # units are in Rg
    isk_2 = inverse_Sk_diblock(x*x,N,chi_a) # units are in Rg
    F = Fitting_Functions(0.5,N)
    
    param,pcov = curve_fit(F.SK_Diblock, x, S)
    param2,pcov2 = curve_fit(F.SK_Diblock_Add, x, S)

    Sfit = F.SK_Diblock(x,param)
    Sfit_add = F.SK_Diblock_Add(x, param2[0],param2[1], param2[2], param2[3])

    plt.plot(x,S,'r^',label = r'$S(q)/CN$')
    plt.plot(x,SAA,'ok',label = r'$S_{AA}(q)/CN$')
    
    plt.ylabel('$S(q)/CN$')
    plt.xlabel('$q \ (R^{-1}_{g,0})$')
    # plt.plot(x,1/isk_2,label =fr'$\chi N = {chi_a*N}$')
    # plt.plot(x,1/isk_1,label =fr'$\chi N = {chi*N:0.2f}$')
    # plt.plot(x,1/isk_1,label =fr'$\chi N = {chi*N:0.2f}$')
    plt.plot(x,Sfit,label =fr'$\chi N = {param[0]*N:0.2f}-RPA-Fit$')
    plt.plot(x,Sfit_add,label =fr'$\chi N = {param2[0]*N:0.2f}-RPA-AugFit$')

    plt.xlim(0,5)
    plt.ylim(0,5)
    plt.legend()
    if i==10:
        break