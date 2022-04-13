#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 11:43:30 2022

@author: tquah
"""
import numpy as np
from scipy.interpolate import interp1d
import re

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


class PeakSK():
    def __init__(self,_f,_N,_chi,_Nbar,_chies=10.495):
        self._f = _f
        self._N = _N
        self._chi = _chi
        self._Nbar = _Nbar
        self._chies = _chies
    # from Qin et al. (JCP 2011) 
    def deltaS_inv_ROL(self):
        def e_fun(_f,_a,_b):
            return _a+_b*np.power(1/2-_f,2)
        _a = e_fun(self._f, 0.652, -0.799)
        _b = e_fun(self._f, 0.0, -18.2)
        _c = e_fun(self._f, 2.5, 0.857)
        _d = e_fun(self._f, -7.58, -1460)
        _den = self._chies-(self._chi*self._N)    
        _pc1 = _a*np.power((self._chi*self._N),_c)/np.sqrt(_den)
        _pc2 = _b*np.power((self._chi*self._N),_c)/_den
        _pc3 = _d
        return 1/np.sqrt(self._Nbar)*(_pc1+_pc2+_pc3)
    def S_inv_RPA(self):
        return 2*(self._chies-(self._chi*self._N))
def MinMax(_v):
    return np.min(_v),np.max(_v)
   
class ROL_Object():
    def __init__(self,_chie,_Sinv,_N):
        self._chie = _chie 
        self._Sinv = _Sinv
        self._f = interp1d(self._chie, self._Sinv)
        self._finv = interp1d(self._Sinv, self._chie)
        self._chimax = np.max(_chie)
        self._chimin = np.min(_chie)
        self._Smax = np.max(_Sinv)
        self._Smin = np.min(_Sinv)
    def Evaluate_Chie(self,_inputchi):
        _loc = np.where((_inputchi>=self._chimin) & (_inputchi<=self._chimax))
        if len(_loc)>0:
            _chi_in = _inputchi[_loc]
            return _inputchi[_loc], self._f(_chi_in),_loc
        
    def Evaluate_Sinv(self,_inputS):
        _loc = np.where((_inputS>=self._Smin) & (_inputS<=self._Smax))
        if len(_loc)>0:
            S_in = _inputS[_loc]
            return _inputS[_loc], self._finv(S_in),_loc

    def MinMax(self):
        return MinMax(self._chie),MinMax(self._Sinv)
        

def inverse_SK_diblock_ROL(f,N,chi,Nbar):
    def e_fun(f,a,b):
        return a+b*np.power(1/2-f,2)
    a = e_fun(f, 0.652, -0.799)
    b = e_fun(f, 0.0, -18.2)
    c = e_fun(f, 2.5, 0.857)
    d = e_fun(f, -7.58, -1460)
    chisN = 10.495
    den = chisN-(chi*N)    
    
    if np.abs(den)<1e-6:
        pc1 = 0
        pc2 = 0
        pc3 = d
    if den<0:
        pc1 = 0
        pc2 = b*np.power((chi*N),c)/den
        pc3 = d
    else:
        pc1 = a*np.power((chi*N),c)/np.sqrt(den)
        pc2 = b*np.power((chi*N),c)/den
        pc3 = d
    
    
    return 1/np.sqrt(Nbar)*(pc1+pc2+pc3)

def create_nmap(_list,_label):
    _map = np.meshgrid(*_list)

    _dict = dict()
    for i in range(len(_map)):
        _dict[_label[i]] = _map[i].flatten()
    return _dict
def generator_ROL_emperical(_lst_array,_label):

    _sweepmap = create_nmap(_lst_array,_label)
    _PeakS = PeakSK(_sweepmap['f'],_sweepmap['N'],_sweepmap['chiN']/_sweepmap['N'],\
           np.power(6,3)*_sweepmap['N']*np.power(_sweepmap['c'],2))

    _dS_ROL = _PeakS.deltaS_inv_ROL()
    _S_RPA = _PeakS.S_inv_RPA()
    _chie = _sweepmap['chiN']+_dS_ROL/2
    _Sinv = _dS_ROL+_S_RPA
    _sweepmap['chie'] = _chie
    _sweepmap['Sinv'] = _Sinv
    return _sweepmap
def convert_dictionary_to_npy(_dict):
    #requires everything to have the same length
    _keys = list(_dict)
    _list = []
    _str = ""

    for _ele in _keys:
        _list.append(_dict[_ele])
        _str+=f"{_ele} "
    _array = np.vstack(_list).transpose()
    return _array,_str
    
def mc_list(_list,_target):
    #prechecks
    assert len(_list)==len(_target), "both objects should be equal length"
    for _ele in _list:
        assert len(_list[0])==len(_ele), "All Arrays in list should be equal length"
    
    _maskarray = np.ones_like(_list[0])
    for i in range(0,len(_list)):
        _tempmask = np.ma.masked_values(_list[i],_target[i])  
        _maskarray*=_tempmask.mask
    return np.where(_maskarray==True)

def Build_ROL(_ROLTXT):
    def RemoveNan(_array):
        _logic = np.isnan(_array)
        return  np.where(_logic==False)
    
    _dictionary = dict()
    _data = np.loadtxt(_ROLTXT)
    # in this dataset 
    _f = open(_ROLTXT)
    _header = _f.readline()
    _label = _header[2:-2].split(' ')
    _floc = _label.index('f')
    _cloc = _label.index('c')
    _Nloc = _label.index('N')
    _chieloc = _label.index('chie')
    _Sinvloc = _label.index('Sinv')

    _f = _data[:,_floc]
    _c = _data[:,_cloc]
    _N = _data[:,_Nloc]
    _chie = _data[:,_chieloc]
    _Sinv = _data[:,_Sinvloc]
    # _Nbar = np.power(6,3)*np.power(_c,2)*_N
    for _i in range(0,len(_data)):
        if [_f[_i],_c[_i],_N[_i]] in list(_dictionary):
            continue
        _fullloc = mc_list([_f,_c,_N],[_f[_i],_c[_i],_N[_i]])
        _chitemp = _chie[_fullloc]
        _Sinvtemp = _Sinv[_fullloc]

        _loc = RemoveNan(_Sinvtemp)
        _chitemp = _chitemp[_loc]
        _Sinvtemp = _Sinvtemp[_loc]

        
        if len(_loc)==0:
            continue
        _dictionary[_f[_i],_c[_i],_N[_i]] = ROL_Object(_chitemp,_Sinvtemp,_N[_i])
    return _dictionary
def ExtractNumber(line):
    return re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)


def extract_numberlist(lst,pos):
    rlst = []
    for ele in lst:
        rlst.append(float(ExtractNumber(ele)[pos]))
    
    return rlst
