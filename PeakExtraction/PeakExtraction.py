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
Path = '/media/tquah/TIMMY/Projects/CL_CG/FORCEMATCH-Smeared_Compressible-Unsmeared-Compressible'

#change directory
os.chdir(Path)

Working_Directory = os.getcwd()
file = open("data.pkl",'rb')

dictionary = pickle.load(file)
header = list(dictionary)
for i in range(0,len(header),1):
    
    S = dictionary[header[i]]['SKAA'][:,1]+ dictionary[header[i]]['SKAA'][:,1]-2* dictionary[header[i]]['SKAB'][:,1]
    plt.figure()
    plt.plot(dictionary[header[i]]['SKAA'][:,0],S,'r^')
    plt.plot(dictionary[header[i]]['SKAA'][:,0],dictionary[header[i]]['SKAA'][:,1],'ok')
    plt.xlim(0,5)
    plt.ylim(0,5)
    break    