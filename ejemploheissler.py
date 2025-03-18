#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:13:37 2025
ejemplo pared de ladrillo transitoria
@author: whadymacbook2016
"""
import numpy as np

k=0.72
rho=1920.0 
cp=835.0 
alfa=k/(rho*cp)

h=50.0 
Tm=300.0 
Ti=20.0 
L=0.15 

Bi=h*L/k

t=3600.0*3 
x1=0.0/L
x2=7.5e-2/L
x3=15.0e-2/L 

Fo=alfa*t/L**2 

epsilon1=1.4289 
c1=4*np.sin(epsilon1)/(2*epsilon1+np.sin(2*epsilon1))

theta1=c1*np.exp(-epsilon1**2*Fo)*np.cos(epsilon1*x1)

T1=theta1*(Ti-Tm)+Tm

theta2=c1*np.exp(-epsilon1**2*Fo)*np.cos(epsilon1*x2)

T2=theta2*(Ti-Tm)+Tm


theta3=c1*np.exp(-epsilon1**2*Fo)*np.cos(epsilon1*x3)

T3=theta3*(Ti-Tm)+Tm

Q0=rho*cp*1*(2*L)*(Ti-Tm)
Q=Q0*(1-np.sin(epsilon1)/epsilon1*theta1)

Fo120=(-1/epsilon1**2)*np.log((1/c1)*(120-Tm)/(Ti-Tm))
t120=Fo120*L**2/alfa








