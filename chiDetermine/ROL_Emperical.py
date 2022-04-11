#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 14:16:10 2022

@author: tquah
"""
import numpy as np

def E_ROL(f,chi,N,Nbar):
    a = 0.652-0.799*(1/2-f)**2
    b = -18.2*(1/2-f)**2 
    