#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 08:33:18 2020

@author: whadymacbook2016
"""
#ejemplo condensacon
import numpy as np
from CoolProp.CoolProp import PropsSI

Ts=-2.0+273.15
Tsat=PropsSI('T','P',101325,'Q',1,'Water')
Tf=(Ts+Tsat)/2.0
g=9.8
L=1
dg=PropsSI('D','T',Tsat,'Q',1,'Water')
dl=PropsSI('D','T',Tf,'Q',0,'Water')
hg=PropsSI('H','T',Tsat,'Q',1,'Water')
hl=PropsSI('H','T',Tsat,'Q',0,'Water')
mul=PropsSI('viscosity','T',Tf,'Q',0,'Water')
kl=PropsSI('conductivity','T',Tf,'Q',0,'Water')
cpl=PropsSI('C','T',Tf,'Q',0,'Water')
hfg=hg-hl
hfgp=hfg+0.68*cpl*(Tsat-Ts)

hprom=0.943*(g*dl*(dl-dg)*kl**3*hfgp/(mul*(Tsat-Ts)*L))**0.25
q=hprom*(Tsat-Ts)
mcond=q/hfgp
print(hprom,q,mcond)