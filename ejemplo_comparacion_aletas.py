#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:46:01 2025

@author: whadymacbook2016
"""
import numpy as np
from scipy.special import k1,k0,i1,i0
k=205 
Tb=100 
Tf=25 
h=50 
thetab=Tb-Tf


#alecta rectangular
L=50e-3 
t=2e-3 
b=20e-3 
P=2*(b+t)
Ac=t*b
As=2*L*b+2*L*t
m=np.sqrt(h*P/(k*Ac))
M=np.sqrt(h*P*k*Ac)*thetab
ql=M*np.tanh(m*L)
effl=ql/(h*As*thetab)
efectl=ql/(h*Ac*thetab)


#aleta anular
r1=10e-3 
r2=25e-3 
Asr=np.pi*(r2**2-r1**2)*2
m=np.sqrt(2*h/(k*t))
qr=2*np.pi*k*r1*t*thetab*m*(k1(m*r1)*i1(m*r2)-i1(m*r1)*k1(m*r2))/\
   (k0(m*r1)*i1(m*r2)+i0(m*r1)*k1(m*r2)) 
effr=qr/(h*Asr*thetab)
efectir=qr/(2*np.pi*r1*t*h*thetab)





















