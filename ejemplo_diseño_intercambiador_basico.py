#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:21:18 2024

@author: whadymacbook2016
"""
import numpy as np
mh=2
cph=3500
Thi=80
Tho=50
Tci=15
mc=2.5
U=2000
cpc=4180
D=2e-2 

Ch=mh*cph
Cc=mc*cpc

q=Ch*(Thi-Tho)


#flujo pararlelo
Tco=q/Cc+Tci
dT1=Thi-Tci
dT2=Tho-Tco
DTln=(dT1-dT2)/np.log(dT1/dT2)


A=q/(DTln*U)


L=A/(np.pi*D)
pasos=L/7


#contraflujo
dT1cf=Thi-Tco
dT2cf=Tho-Tci
DTlncf=(dT1cf-dT2cf)/np.log(dT1cf/dT2cf)

Acf=q/(DTlncf*U)
Lcf=Acf/(np.pi*D)
pasoscf=Lcf/7


#NTU contraflujo
Cmin=Ch
eff=q/(Cmin*(Thi-Tci))
Cr=Cmin/Cc

NTUcf=(1/(Cr-1))*np.log((eff-1)/(eff*Cr-1))
A_NTUcf=NTUcf*Cmin/U















