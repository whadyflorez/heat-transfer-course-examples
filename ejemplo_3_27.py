#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 09:21:30 2026

@author: whadymacbook2016
"""
import numpy as np
qg=30000 
ho=1000
To=20
Rc=1e-4
Lb=5e-3 
kb=1 
hi=40
Ti=20
Rcvo=1/ho
Rcn=Lb/kb
Rcvi=1/hi

A=np.zeros((5,5))
RHS=np.zeros(5)

A[0,:]=[0,0,0,1,1]
A[1,:]=[1/Rcvo,0,0,-1,0]
A[2,:]=[1/(Rc+Rcn+Rcvi),0,0,0,-1]
A[3,:]=[1/Rc,-1/Rc,0,0,-1]
A[4,:]=[0,1/Rcn,-1/Rcn,0,-1]
RHS[0]=qg
RHS[1]=To/Rcvo
RHS[2]=Ti/(Rc+Rcn+Rcvi)
RHS[3]=0
RHS[4]=0

X=np.linalg.solve(A,RHS)














