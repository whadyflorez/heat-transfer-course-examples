#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 07:28:00 2026

@author: whadymacbook2016
"""
import numpy as np

#problema 3.2

L=200e-3
kc=0.15
kps=0.04
kp=1.4
Ti=20 
Te=0  

delta=(L/kc-L/kp)*kps

Rc=L/kc
q=(Ti-Te)/Rc


#problema 3.3
Tin=20
Tout=-10 
k_vidrio=0.8
d_vidrio=4e-3 
hi=30 
ho=65 
Req=1/hi+d_vidrio/k_vidrio+1/ho 
q_vidrio=(Tin-Tout)/Req


Ts1=Tin-q_vidrio*(1/hi)
Ts2=q_vidrio*(1/ho)+Tout


























