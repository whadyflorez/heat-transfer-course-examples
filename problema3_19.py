#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 08:05:05 2026

@author: whadymacbook2016
"""
import numpy as np

k=0.05
Tfi=300 
Tfo=25 
hi=30 
ho=10
qr=100 
To=40 

Rcv=1/hi 
Ro=1/ho

qp=(To-Tfo)/Ro
Ti=Tfi-(qp-qr)*Rcv
Rcd=(Ti-To)/qp
L=k*Rcd