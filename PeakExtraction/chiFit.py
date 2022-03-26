#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 13:33:18 2022

@author: tquah
"""
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import re
import os


def ExtractNumber(line):
    return re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)

Nbarfiles = glob.glob('ROL')