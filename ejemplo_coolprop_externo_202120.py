#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 08:58:02 2021

@author: whadyimac
Ejemplo coeficientes de conveccion externos sobre placas
"""
from CoolProp.CoolProp import PropsSI
import CoolProp as CP
import numpy as np

#http://www.coolprop.org/coolprop/examples.html
#http://www.coolprop.org/coolprop/HighLevelAPI.html#propssi-function
#http://www.coolprop.org/coolprop/HighLevelAPI.html#fluid-information

Tm=20+273.15
Tw=80+273.15
Tf=0.5*(Tm+Tw)
L=2
v=0.1
dens_air=PropsSI('D','T',Tf,'P',101325,'air'); print(dens_air)
visc_air=PropsSI('V','T',Tf,'P',101325,'air'); print(visc_air)
Pr_air=PropsSI('PRANDTL','T',Tf,'P',101325,'air'); print(Pr_air)
k_air=PropsSI('CONDUCTIVITY','T',Tf,'P',101325,'air'); print(Pr_air)
Re=dens_air*v*L/visc_air
print('Re',Re)
Nu=0.664*Re**0.5*Pr_air**(1/3)
print('Nu',Nu)
#Nu=(1/L)*0.0296*dens_air*v/visc_air*Pr_air**(1/3)*L**2/2
h=Nu*k_air/L
print('h',h)
flux_q=h*(Tw-Tm)
print('q',flux_q)

#caso turbulento
v=10
Re=dens_air*v*L/visc_air
print('Re',Re)
Nu=(0.037*Re**(4/5)-871)*Pr_air**(1/3)
h=Nu*k_air/L
print('h',h)
flux_q_turb=h*(Tw-Tm)