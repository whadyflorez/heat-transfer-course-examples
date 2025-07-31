#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 12:39:26 2025

@author: whadymacbook2016
"""
import numpy as np
d1=5e-3
d2=10e-2 
h1=100
h2=20
Tg=600
Tm=25
kl=0.72
ka=26.3e-3
A=1

Rconv1=1/(h1*A)
Rcondl=d1/(kl*A)
Rconda=d2/(ka*A)
Rconv2=1/(h2*A)

Req=Rconv1+Rcondl+Rconda+Rcondl+Rconv2

q=(Tg-Tm)/Req

Rwlad=(2*d1+d2)/(kl*A)

qlad=(Tg-Tm)/(Rconv1+Rwlad+Rconv2)

Ts1=Tg-q/(h1*A)
Ts2=Ts1-q*Rcondl