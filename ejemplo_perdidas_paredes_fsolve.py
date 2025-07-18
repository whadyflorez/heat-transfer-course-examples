#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 08:46:44 2021

@author: whadyimac
"""
#EJEMPLO PERDIDAS EN PAREDES DE HORNOS CON CONVECCION Y RADIACION BALANCES SIMPLES

import numpy as np
from scipy.optimize import fsolve

k=0.1 #W/(m.K)
L=0.5 #m
depth=1.0
height=1.0
s=5.67e-8
e=0.8
As=height*depth
At=As
hin=100.0 #W/(m^2.K)
hout=10.0
Tin=800.0+273.15
Tout=25+273.15

def balances(T):
    B=np.zeros(2)
    B[0]=hin*As*(Tin-T[0])+s*e*As*(Tin**4-T[0]**4)-k*At*(T[0]-T[1])/L
    B[1]=k*At*(T[0]-T[1])/L-hout*As*(T[1]-Tout)-s*e*As*(T[1]**4-Tout**4)
    return B

Ts=fsolve(balances,[750.0+273.15,30.0+273.15])

qconv_int=hin*As*(Tin-Ts[0])
qrad_int=s*e*As*(Tin**4-Ts[0]**4)
qcond=k*At*(Ts[0]-Ts[1])/L
qconv_out=hout*As*(Ts[1]-Tout)
qrad_ext=s*e*As*(Ts[1]**4-Tout**4)


print('Ts=',Ts)

    




