#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 11:02:52 2024

@author: whadymacbook2016
"""
import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp as CP

D=76.45e-3  
L=332.56e-3 
g=9.8
Tm=2+273.15
Ts=20+273.15 
Tf=0.5*(Tm+Ts)

rho=PropsSI('D','T',Tf,'P',101325,'air')
mu=PropsSI('V','T',Tf,'P',101325,'air')
Pr=PropsSI('PRANDTL','T',Tf,'P',101325,'air')
k=PropsSI('CONDUCTIVITY','T',Tf,'P',101325,'air')
cp=PropsSI('C','T',Tf,'P',101325,'air')
beta=PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',Tf,'P',101325,'air')

nu=mu/rho
alfa=k/(rho*cp)

RaL=g*beta*(Ts-Tm)*L**3/(nu*alfa)
RaD=g*beta*(Ts-Tm)*D**3/(nu*alfa)
NuL=(0.825+0.387*RaL**(1/6)/(1+(0.492/Pr)**(9/16))**(8/27))**2
hL=NuL*k/L

NuD=(0.6+0.387*RaL**(1/6)/(1+(0.559/Pr)**(9/16))**(8/27))**2
hD=NuD*k/D

















