#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 08:24:54 2026

@author: whadymacbook2016
ejemplo radiacion entre cuerpos grises 2026
2 superficies en esquina
"""
import numpy as np

A0=0.2*0.1
A1=0.1*0.1
A2=np.sqrt(0.1**2+0.2**2)*0.1

T0=100.0+273.15
T1=60.0+273.15
T2=25.0+273.15

s=5.67e-8
e0=0.8
e1=0.7

Eb0=s*T0**4
Eb1=s*T1**4
Eb2=s*T2**4

F=np.zeros((3,3))

F[0,0]=0.0 
F[0,1]=0.116 
F[0,2]=1.0-F[0,0]-F[0,1]

#reciprocidad A0F01=A1F10 
F[1,0]=(A0/A1)*F[0,1]
F[1,1]=0.0 
F[1,2]=1.0-F[1,1]-F[1,0]

#reciprocidad 
F[2,2]=0.0 
F[2,0]=A0*F[0,2]/A2
#F[2,1]=A1*F[1,2]/A2
F[2,1]=1.0-F[2,0]-F[2,2]

#ecs de radiacion entre cuerpos grises
#(Eb0-J0)/(1-e0)/e0A0=(J0-J1)/(A0F01)^-1+(J0-Eb2)/(A0F02)^-1
#(Eb1-J1)/(1-e1)/e1A1=(J1-J0)/(A1F10)^-1+(J1-Eb2)/(A1F12)^-1
#J2=Eb2

A=np.zeros((2,2))
B=np.zeros(2)

A[0,0]=A0*F[0,1]+A0*F[0,2]+e0*A0/(1.0-e0)
A[0,1]=-A0*F[0,1]
B[0]=Eb0*e0*A0/(1.0-e0)+Eb2*A0*F[0,2]
A[1,0]=-A1*F[1,0]
A[1,1]=A1*F[1,0]+A1*F[1,2]+e1*A1/(1.0-e1)
B[1]=Eb1*e1*A1/(1.0-e1)+Eb2*A1*F[1,2]

J=np.linalg.solve(A,B)
J2=Eb2

q0=A0*F[0,1]*(J[0]-J[1])+A0*F[0,2]*(J[0]-J2)
q1=A1*F[1,0]*(J[1]-J[0])+A1*F[1,2]*(J[1]-J2)























