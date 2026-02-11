#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 08:25:38 2026

@author: whadymacbook2016
"""
import numpy as np

kb=200
ka=0.5
r1=0.005
r2=0.02
r3=0.03
Tfi=40 
Tfo=25 
hi=100
ho=20
qh_m2=15

Rci=1/(hi*2*np.pi*r1)
RcdB=np.log(r2/r1)/(2*np.pi*kb)
RcdA=np.log(r3/r2)/(2*np.pi*ka)
Rco=1/(ho*2*np.pi*r3)

A=np.zeros((3,3))
RHS=np.zeros(3)

A[0,:]=[1,0,-1/(Rci+RcdB)]
A[1,:]=[0,1,-1/(RcdA+Rco)]
A[2,:]=[1,1,0]


RHS[0]=-Tfi/(Rci+RcdB)
RHS[1]=-Tfo/(Rco+RcdA)
RHS[2]=qh_m2*2*np.pi*r2

sol=np.linalg.solve(A,RHS)


