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

class Fitting_Functions():
    def __init__(self,f,N):
        self.f = f
        self.N = N
    def SK_Diblock(self,x,a):
        xx =x*x
        S_AA,S_BB,S_AB = StructureFactorMatrix(xx,self.f)
        detS = (S_AA*S_BB)-(S_AB*S_AB)
        return 1/((g(xx,1)/self.N/detS)-2*a)

    def SK_Diblock_Add(self,x,a,b,c,d,e):
        xx =x*x
        S_AA,S_BB,S_AB = StructureFactorMatrix(xx,self.f)
        detS = (S_AA*S_BB)-(S_AB*S_AB)
        return (1/(((g(xx,1)/self.N/detS)-2*a)+b+c*np.power(x,2)+d*np.power(x,4)))+e






def inverse_SK_diblock_ROL(f,N,chi,Nbar):
    def e_fun(f,a,b):
        return a+b*np.power(1/2-f,2)
    a = e_fun(f, 0.652, -0.799)
    b = e_fun(f, 0.0, -18.2)
    c = e_fun(f, 2.5, 0.857)
    d = e_fun(f, -7.58, -1460)
    chisN = 10.495
    den = chisN-(chi*N)
    
    pc1 = a*np.power((chi*N),c)/np.sqrt(den)
    pc2 = b*np.power((chi*N),c)/den
    pc3 = d
    
    return 1/np.sqrt(Nbar)*(pc1+pc2+pc3)