#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:24:58 2024

@author: whadymacbook2016

ejemplo flujo internos nusselt 2024
"""
import numpy as np

m=0.5
Ti=25+273.15
Ts=100+273.15
mu=48.6e-2
k=145e-3
cp=1909
D=25e-3
L=100

Re=4*m/(np.pi*D*mu)
Nu=3.66
h=Nu*k/D

To=Ts+(Ti-Ts)*np.exp(-h*np.pi*D*L/(m*cp))

Tprom=0.5*(Ti+To)
cp=0.5*(1909+1951)
mu=0.5*(48.6e-2+25.3e-2)
Re=4*m/(np.pi*D*mu)
h=Nu*k/D

dT1=373-Ti
dT2=373-To
dTln=(dT1-dT2)/np.log(dT1/dT2)
A=np.pi*D*L
Req=1/(h*A)
UA=1/Req
q=UA*dTln















