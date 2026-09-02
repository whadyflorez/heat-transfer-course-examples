#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:22:46 2026

@author: whadymacbook2016
"""
import numpy as np
#horno 50 cm lado cubico. Espesor=5 cm
k=0.72
Ti=300.0  #C
To=25.0 

#perdidas por las paredes
A=0.4*0.4 #m2
qw=k*A/0.05*(Ti-To)
qw_tot=6*qw

#perdidas en bordes
Sb=0.54*0.4
qb=Sb*k*(Ti-To)
qb_tot=12*qb

#perdida por esquinas
Se=0.15*0.05 
qe=Se*k*(Ti-To)
qe_tot=8*qe

