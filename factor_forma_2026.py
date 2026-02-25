#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 08:19:16 2026

@author: whadymacbook2016
"""
import numpy as np
z=0.5 
L=2.0 
T1=25
m=0.01
Cp=4180.0
D=0.025
k=1.5
T2=80.0

F=3*D/2
S=2*np.pi*L/np.log(4*z/D)

q=S*k*(T2-T1)
Tout=q/(m*Cp)+T1

#correccion
Ttubo=0.5*(T1+Tout)
qcorregido=S*k*(T2-Ttubo)
Tout_corregida=qcorregido/(m*Cp)+T1



