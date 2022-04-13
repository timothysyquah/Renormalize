#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 09:08:41 2022

@author: tquah
"""
import matplotlib.pyplot as plt
import numpy as np
import pickle
import re
plt.close('all')
def extract_number(line):
    return re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)

def extract_numberlist(lst,pos):
    rlst = []
    for ele in lst:
        rlst.append(float(extract_number(ele)[pos]))
    
    return rlst


chi_column = 0
fm_chi_column = 1
n_column = 2
qstar_column = 3
structure_column = 4
chicutoff = 2.0
C = 1.0
k_cutoff = 2.0 #rg
dataname = 'CL_dataset_large.dat'
cl_data = np.loadtxt(dataname)


uniquen = np.unique(cl_data[:,n_column])



plt.figure()
for i in range(0,len(uniquen)):
    loc = np.where(cl_data[:,n_column]==uniquen[i])[0]
    tempdata = cl_data[loc,:]
    chisort = np.argsort(tempdata[:,chi_column])
    tempdata = tempdata[chisort,:]
    plt.scatter(tempdata[:,chi_column],np.power(tempdata[:,n_column],1.)/tempdata[:,structure_column]/2)


#plot structure factor by N

file = open("dataC1_zeta10.pkl",'rb')
dictionary = pickle.load(file)
header = list(dictionary)
nlist = extract_numberlist(header, 0)

for i in range(0,len(uniquen)):
    loc = np.where(nlist==uniquen[i])[0]
    plt.figure()
    plt.title(f'N = {uniquen[i]}')

    for j in range(0,len(loc)):
        chi = dictionary[header[loc[j]]]['Chi']
        
        if chi*uniquen[i]<chicutoff:
            continue
        
        SKAA = dictionary[header[loc[j]]]['SKAA']
        SKBB = dictionary[header[loc[j]]]['SKBB']
        SKAB = dictionary[header[loc[j]]]['SKAB']
        kcutoff = np.where(SKAA[:,0]<2.0)[0]

        x = SKAA[kcutoff,0]*np.sqrt(uniquen[i])
        S = SKAA[:,1]+SKBB[:,1]-2*SKAB[:,1]
        plt.plot(x,S[kcutoff]/uniquen[i]/4,label = f'$\chi N ={chi * uniquen[i]}$')
        plt.plot(x,SKAA[kcutoff,1]/uniquen[i],label = f'$\chi N ={chi * uniquen[i]}$')

    plt.legend()
    plt.xlabel('$k R_g$')
    plt.ylabel("$S(k)/cN$")
    plt.savefig(f"/home/tquah/Figures/C_{C}_N_{uniquen[i]}_S.pdf",dpi = 300)

#try to fit to ROL








