#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 11:14:35 2021

@author: whadyimac
"""

#EJEMPLO CALENTAMIENTO FOCO
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

h=5.0
Tamb=25.0+273.15
s=5.67e-8
e=0.9
P=100.0
r=4.0e-2
As=4.0*np.pi*r**2
m=107.7e-3
cp=835.0
tfinal=20*60.0

def dTdt(T,t):
    dy=(-h*As*(T[0]-Tamb)-s*e*As*(T[0]**4-Tamb**4)+P)/(m*cp)
    return dy

t=np.linspace(0,tfinal,100)

sol = odeint(dTdt, [Tamb], t)

qperdidas=np.zeros(100)
for i in range(0,100):
    qperdidas[i]=s*e*As*(sol[i]**4-Tamb**4)


plt.plot(t,sol,'o-')
plt.figure()
plt.plot(t,qperdidas,'--')
