#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 08:25:40 2026

@author: whadymacbook2016
"""
import numpy as np

k=186 
hg=35
hc=20 
Tg=457
Tc=80
L=0.6 
At=6e-3*2 
Af=0.6*2*2
d=6e-3 
e=0.05 
Lc=L+0.5*d
Ap=Lc*d
parametro_x=Lc**(3.0/2.0)*(hg/(k*Ap))**(1.0/2.0)
P=2*2+2*6e-3
m=np.sqrt(hg*P/(k*At))
eff=np.tanh(m*Lc)/(m*Lc)
Rf=1/(eff*hg*Af)
N=1/0.05-1
A_sf=1*2-N*Af
A_espacio=2*0.05
Abase=2*1

Rtotf=1/(eff*hg*Af*N)
Rtotsf=1/((N+1)*hg*A_espacio)
Rtot_paralelo=1/((1/Rtotf)+(1/Rtotsf))

Rconv=1/(hc*Abase)

Req=Rconv+Rtot_paralelo

q=(Tg-Tc)/Req

qsf=(Tg-Tc)/((1/(hc*Abase))+(1/(hg*Abase)))




















