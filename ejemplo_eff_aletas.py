#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 08:51:20 2026

@author: whadymacbook2016
"""
import numpy as np
from scipy.special import k0,k1,i0,i1

r1=0.01
r2=0.03
d=2e-3
e=5e-3
h=10
Tf=25
Tb=180
thetab=Tb-Tf
N=10
k=400
Af=2*np.pi*(r2**2-r1**2)+2*np.pi*r2*d

m=np.sqrt(2*h/(k*d))

qf=2*np.pi*k*r1*d*thetab*m*(k1(m*r1)*i1(m*r2)-i1(m*r1)*k1(m*r2))/\
(k0(m*r1)*i1(m*r2)+i0(m*r1)*k1(m*r2))

qmax=h*Af*thetab

nu_def=qf/qmax

nu_for=(2*r1/(m*(r2**2-r1**2)))*(k1(m*r1)*i1(m*r2)-i1(m*r1)*k1(m*r2))/\
(k0(m*r1)*i1(m*r2)+i0(m*r1)*k1(m*r2))

qe=h*(2*np.pi*r1*e)*thetab*(N-1)

qtot=N*qf+qe

eff=qf/(h*2*np.pi*r1*d*thetab)
    













