#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 08:10:32 2026

@author: whadymacbook2016
"""
import numpy as np

k=237
t=1e-3 
L=12e-3  
d=100e-3 #profundidad
a=200e-3 #ancho
nf=50
s=a/(nf-1) #longitud de los espacios
T0=400
TL=350
Tm=300
h=150 
At=t*d
As=2*L*d
Ae=s*d  #espacios entre aletas
P=2*(t+d)
thetab=T0-Tm
thetaL=TL-Tm

M=np.sqrt(h*P*k*At)*thetab
m=np.sqrt((h*P)/(k*At))
qf=M*(np.cosh(m*L)-thetaL/thetab)/np.sinh(m*L) #1 aleta
Qf=nf*qf #calor total de las 50 aletas

#calor de los espacios
qtop=h*Ae*(T0-Tm)
qbottom=h*Ae*(TL-Tm)
Qtop=qtop*(nf-1)
Qbottom=qbottom*(nf-1)

Qtot=Qf+Qbottom+Qtop

qmax=h*As*(T0-Tm)
effic=qf/qmax





















