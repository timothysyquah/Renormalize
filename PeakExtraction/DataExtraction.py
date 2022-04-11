#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:07:40 2022

@author: tquah
"""


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
#Path = '/media/tquah/TIMMY/Projects/CL_CG/FORCEMATCH-Smeared_Compressible-Unsmeared-Compressible'
Path = "/media/tquah/TIMMY/Projects/PodProjects/Projects/CG/CG"
#change directory
os.chdir(Path)

Working_Directory = os.getcwd()

Directories = glob.glob('./**/chi*',recursive=True)


#remove .dat stuff using depth
Directories = Prune(Directories,2)

Data_Dictionary = dict()

for i in range(0,len(Directories)):
    os.chdir(Directories[i])
    #extract 
    temp_details = Directories[i].split('/')
    N = float(ExtractNumber(temp_details[1])[0])
    chi = float(ExtractNumber(temp_details[2])[0])
    
    if np.loadtxt('STATUS')==1:
        os.chdir(Working_Directory)

        continue
    
    
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]] = dict()

    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['Chi'] = chi/N
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['N'] = N
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['Operators'] = np.loadtxt('model1_operators.dat')

    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['SKAA'] = np.loadtxt('model1_SKAA.dat')
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['SKAB'] = np.loadtxt('model1_SKAB.dat')
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['SKBB'] = np.loadtxt('model1_SKBB.dat')
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['Chi_I'] = np.loadtxt('chiN.dat')

    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['Status'] = np.loadtxt('STATUS')
    Data_Dictionary[temp_details[1]+'_'+temp_details[2]]['OP'] = np.loadtxt('model1_OrientationPersistenceOP.dat')
    os.chdir(Working_Directory)
    
f = open('data.pkl',"wb")
pickle.dump(Data_Dictionary,f)
f.close()
    
