#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 07:47:24 2026

@author: whadymacbook2016
"""
from sympy import symbols, linsolve,latex
from sympy import init_printing
from IPython.display import display

init_printing(use_latex=True)


# 1. Definir símbolos
c1, c2,Ta,Tb,L,q,k = symbols('c1 c2 Ta Tb L q k')

# 2. Definir el sistema (se asume que cada expresión es igual a 0)
# Ejemplo: x + y = 5  y  x - y = 1
sistema = [Ta - c2, Tb+q*L**2/(2*k)-c1*L-c2]

# 3. Resolver
solucion = linsolve(sistema, (c1, c2))
print(solucion)  # Salida: {(3, 2)}
display(solucion)
solucion_latex=latex(solucion)
