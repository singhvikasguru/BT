# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 17:32:19 2019

@author: King
"""
import numpy as np
import cmath
import matplotlib.pyplot as plt

kappa_0=2.5
tau_0=0.5
sigma=0
t=np.logspace(-2, 2,num=100, base=10.0)
R_sq=[]
R_sq_t=[]
p=0
fp=(p+2)*((p+2+sigma)*(p+2+sigma)+tau_0*tau_0)+kappa_0*kappa_0*(p+2+sigma)
C_inf=2*((2+sigma)**2+tau_0**2)/fp
C_0=-0.5*C_inf*C_inf+2*kappa_0*kappa_0*((2+sigma)**2-tau_0**2)/(fp*fp)
r=[-2, -2-2.54951j, -2+2.54951j]
C_1=2*((2+sigma-r[0])**2+tau_0**2)/(r[0]*r[0]*(r[0]-r[1])*(r[0]-r[2]))
C_2=2*((2+sigma-r[1])**2+tau_0**2)/(r[1]*r[1]*(r[1]-r[2])*(r[1]-r[0]))
C_3=2*((2+sigma-r[2])**2+tau_0**2)/(r[2]*r[2]*(r[2]-r[1])*(r[2]-r[0]))
for x in t:
    R_sq.append((C_inf*x+C_0+ C_1*cmath.exp(r[0]*x)+ C_2*cmath.exp(r[1]*x) + C_3*cmath.exp(r[2]*x)).real)
    R_sq_t.append(R_sq[-1]/x)
plt.semilogx(t, R_sq)
