#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 11:43:30 2022

@author: tquah
"""
import numpy as np
def inverse_Sk(k,f,b,N,chi):
    num = 1+(1/18)*(k**2)*(b**2)*N
    den = N*f*(1-f)
    return num/den-2*chi
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
