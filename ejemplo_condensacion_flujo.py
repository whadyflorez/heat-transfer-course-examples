#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:11:54 2024

@author: whadymacbook2016
"""
import numpy as np
from CoolProp.CoolProp import PropsSI

Ts=25+273.15
p=0.15*1e5
Tsat=PropsSI('T','P',p,'Q',1,'Water')
Tf=(Ts+Tsat)/2.0
g=9.8
D=6e-3
dg=PropsSI('D','T',Tsat,'Q',1,'Water')
dl=PropsSI('D','T',Tf,'Q',0,'Water')
hg=PropsSI('H','T',Tsat,'Q',1,'Water')
hl=PropsSI('H','T',Tsat,'Q',0,'Water')
hfg=hg-hl
cpl=PropsSI('C','T',Tf,'Q',0,'Water')
hfgp=hfg+0.68*cpl*(Tsat-Ts)
mul=PropsSI('viscosity','T',Tf,'Q',0,'Water')
kl=PropsSI('conductivity','T',Tf,'Q',0,'Water')

Nu=0.729*((dl*g*(dl-dg)*hfgp*D**3)/(mul*kl*(Tsat-Ts)))**0.25
h=Nu*kl/D
Atot=100*np.pi*D
q=h*Atot*(Tsat-Ts)
m=q/hfgp













