#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:41:50 2026

@author: whadymacbook2016
"""
import numpy as np

rho=2700.0 
c=950.0 
k=240.0
Ti=25.0 
Tg=300.0 
h=75.0 
R=75e-3/2 
As=4.0*np.pi*R**2 
V=(4.0/3.0)*np.pi*R**3 
Lc=(1.0/3.0)*R

Bi=h*Lc/k
t=-rho*V*c/(h*As)*np.log(0.1)

Q=-rho*V*c*(Ti-Tg)*(np.exp(-h*As*t/(rho*V*c))-1)






 
