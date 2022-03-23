#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:07:16 2022

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt


def inverse_Sk(k,f,b,N,chi):
    num = 1+(1/18)*(k**2)*(b**2)*N
    den = N*f*(1-f)
    return num/den-2*chi

# import numpy as np
def h(x,s):
    return (1-np.exp(-x*s))/x
#x=k^2b^2 N/6 or (Rg^2)(k^2)

def g(x,s):
    return 2*(np.exp(-x*s)+s*x-1)/x/x

def StructureFactorMatrix(x,f):
    SAA = g(x,f)
    SBB = g(x,1-f)
    SAB = h(x,f)*h(x,1-f)
    return SAA,SBB,SAB

def inverse_Sk_diblock(x,N,chi):
    S_AA,S_BB,S_AB = StructureFactorMatrix(x,0.5)
    detS = (S_AA*S_BB)-(S_AB*S_AB)
    return (g(x,1)/N/detS)-2*chi

N = 1
b = 1
f = 0.5
k  = np.linspace(0,3,100)
chi =np.array([1,1.5,1.8])

plt.close('all')
plt.figure()



for i in range(len(chi)):
    isk = inverse_Sk(k,f,b,N,chi[i])
    plt.plot(k*b*np.sqrt(N),1/isk/N,label = f'$\chi= {chi[i]}$')

plt.ylim(0,3)
plt.legend()
plt.xlabel('$kbN^{1/2}$')
plt.ylabel('$S(k)/N$')
plt.tight_layout()


def h(x,s):
    return (1-np.exp(-x*s))/x
#x=k^2b^2 N/6 or (Rg^2)(k^2)

def g(x,s):
    return 2*(np.exp(-x*s)+s*x-1)/x/x

def StructureFactorMatrix(x,f):
    SAA = g(x,f)
    SBB = g(x,1-f)
    SAB = h(x,f)*h(x,1-f)
    return SAA,SBB,SAB

def inverse_Sk_diblock(x,N,chi):
    S_AA,S_BB,S_AB = StructureFactorMatrix(x,0.5)
    detS = (S_AA*S_BB)-(S_AB*S_AB)
    return (g(x,1)/N/detS)-2*chi

plt.close('all')
plt.figure()



x = np.linspace(0,100,1000)
chi = np.array([10.])
# chi = np.round(np.arange(10,10.51,0.1),2)
for i in range(len(chi)):

    isk = inverse_Sk_diblock(x,N,chi[i])
    plt.plot(x*((N/6)**(1/2)),1/isk/N/np.sqrt(N)*np.sqrt(6)*6,label = f'$\chi= {chi[i]}$')

# plt.ylim(0,0.05)
plt.xlim(0,12)
plt.legend()
plt.xlabel('$ka N^{1/2}$')
plt.ylabel('S(k)/N')