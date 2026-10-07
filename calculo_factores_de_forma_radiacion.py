#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 08:07:51 2026

@author: whadymacbook2016
"""
#ejemplo factores de forma entre 2 discos paralelos

import numpy as np

F=np.zeros((3,3))
#cuerpo 0 disco inferior
#cuerpo 1 disco superior
#cuerpo 2 medio ambiente

r0=0.02
r1=0.01 
L=0.1 
R0=r0/L
R1=r1/L 
A0=np.pi*r0**2
A1=np.pi*r1**2
g=np.sqrt((r0-r1)**2+L**2)
A2=np.pi*(r0**2+r1**2+(r0+r1)*g)

S=1.0+(1.0+R1**2)/R0**2
F[0,1]=0.5*(S-(S**2-4.0*(r1/r0)**2)**0.5)

S=1.0+(1.0+R0**2)/R1**2
F[1,0]=0.5*(S-(S**2-4.0*(r0/r1)**2)**0.5)

F[0,0]=0.0 
F[1,1]=0.0 

#regla de suma F00+F01+F02=1
F[0,2]=1.0-F[0,0]-F[0,1]
#regla de suma F10+F11+F12=1
F[1,2]=1.0-F[1,0]-F[1,1]

#reciprocidad A0F02=A2F20
F[2,0]=(A0/A2)*F[0,2]

#reciprocidad A0F01=A1F10
F[1,0]=(A0/A1)*F[0,1]

#A2F21=A1F12
F[2,1]=(A1/A2)*F[1,2]

F[2,2]=1.0-F[2,0]-F[2,1]






































