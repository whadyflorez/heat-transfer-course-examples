#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 07:51:30 2026

@author: whadymacbook2016
"""
import numpy as np

rho=700.0 
Lc=1.0e-3/2 
cp=2400.0 
k=0.34 
Ti=20.0 
Tf=220.0 
Tinf=300.0
h=55.0  
Lh=3.0 

Bi=h*Lc/k 
alpha=k/(rho*cp)
theta=(Tf-Tinf)/(Ti-Tinf)

Fo=-np.log(theta)/Bi
t=Fo*Lc**2/alpha
V=Lh/t

Fo_tau=-np.log(0.368)/Bi
tau=Fo_tau*Lc**2/alpha



