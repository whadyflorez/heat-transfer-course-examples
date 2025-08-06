#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 13:10:45 2025

@author: whadymacbook2016
"""
import numpy as np
from scipy.optimize import fsolve

qg=10
h=10
k=0.1 
e=0.8 
Tm=25+273.15 
Ts=40+273.15
d=2e-3 
L=5e-2
A=L**2
s=5.67e-8 

def f(x):
    T=x[0]
    Rconv=1/(h*A)
    Rrad=1/(e*s*(T+Tm)*(T**2+Tm**2)*A)
    Rcond=d/(k*A)
    desbalance=qg-(T-Tm)/Rconv-(T-Tm)/Rrad-(T-Ts)/Rcond
    return desbalance

T=[100+273.15]
sol=fsolve(f,T)

T=80+273
Rconv=1/(h*A)
Rrad=1/(e*s*(T+Tm)*(T**2+Tm**2)*A)
Rcond=d/(k*A)
qg=(T-Tm)/Rconv+(T-Tm)/Rrad+(T-Ts)/Rcond




