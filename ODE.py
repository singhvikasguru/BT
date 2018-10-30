# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 17:58:06 2018

@author: PRAVESH
"""

import numpy as np
tanh=np.tanh
import scipy.integrate
#from scipy.integrate import odeint
import matplotlib.pyplot as plt
#from Functions import RK4, fermenter_RK
import math
pow=math.pow
from collections import Counter
exp=np.exp
import csv
import itertools
from collections import defaultdict
import time
from numpy import genfromtxt

c = genfromtxt('Temperature.csv', delimiter=',')
c_x=c[:, 0]
c_y=c[:, 1]

rho_p=1200
f=0.58
Kd_0=6.3*pow(10, 16)
Kpo_0=298.0*pow(10, 2)
Ktdo_0=588.0*pow(10, 4)
Ktm_0=279.6*pow(10, 6)
Ktc=0
Ed=128.54*1000
Ep=18.22*1000
Etd=2.937*1000
MWm=0.10013
MWl=0.0680
Ktc=0
A11=-2.293E-2
A12=7.973E-1
A21=-2.379E-1
A22=8.784
A31=1.385E-1
A32=-3.609
A41=-1.745E-3
A42=2.223
B11=-4.02E-3
B12=1.246E-1
B21=-3.054E-1
B22=9.734
B31=2.706E-1
B32=8.278
B41=2.493E-1
B42=-4.459E-1
Rg=8.314
I0=.258
M0=15000
Rli=0
Rlm=1000
Rvm=0
#%%
def odeint(fun, t, y0, u, integrator=0, method=0, rtol=0, atol=1e-12):
    from scipy.integrate import ode
    s = ode(fun).set_f_params(u)
    integrators = ['vode', 'zvode', 'lsoda', 'dopri5', 'dop853']
    methods = ['adams', 'bdf']
    s.set_integrator(integrators[2],
                     method=methods[1],
#                     order=15,
                     nsteps=30,
                     rtol=0.001)
#                     atol=atol)
#                     with_jacobian=False)
    t0 = t[0]
    dt = t[1]-t0
    y = [y0]
    s.set_initial_value(y0, t0)
    while s.successful() and s.t < t[-1]:
        s.integrate(s.t+dt)
        y.append(s.y)
    y = np.array(y)
    return y

#%%
def PMMA(t, x, T):
    rho_m=966.5-1.1*(T-273.15)
#    x[11]=(x[1]*MWm)/rho_m+(x[9]-x[1])*MWm/rho_p
    shy_m=((x[1]*MWm)/rho_m)/((x[1]*MWm)/rho_m+(x[9]-x[1])*MWm/rho_p)
    Kd=Kd_0*np.exp(-Ed/(Rg*T))
    Kp_0=Kpo_0*np.exp(-Ep/(Rg*T))
    Ktd_0=Ktdo_0*exp(-Etd/(Rg*T))
#    Mw=(MWm)*(x[5]+x[8])/(x[4]+x[7])
    xm=1-x[1]/x[10]
    A1=A11*(T-273.15)+A12
    A2=A21*(T-273.15)+A22
    A3=A31*(T-273.15)+A32
    A4=A41*(T-273.15)+A42
    B1=B11*(T-273.15)+B12
    B2=B21*(T-273.15)+B22
    B3=B31*(T-273.15)+B32
    B4=B41*(T-273.15)+B42
#    if xm>0:
    Kt=Ktd_0*np.exp(A1+A2*xm+A3*xm*xm+A4*xm*xm*xm)
    Kp=Kp_0*np.exp(B1+B2*xm+B3*xm*xm+B4*xm*xm*xm)
    Ki=Kp
    Ktd=Kt
    Ktm=Ktm_0*(Kp/Kp_0)
    xdot=[0 for var in range(0, 16)]
    xdot[0]=-Kd*x[0]+Rli
    xdot[1]=-(x[14]+Ktm)*(x[1]*x[3]/x[11])+Rlm-Rvm-Ki*(x[2]*x[1]/x[11])
    xdot[2]=2*f*Kd*x[0]-Ki*(x[2]*x[1]/x[11])
    xdot[3]=x[14]*(x[2]*x[1]/x[11])-Kt*(x[3]**2/x[11])
    xdot[4]=x[14]*(x[2]*x[1]/x[11])+Kp*x[1]*(x[3]/x[11])-x[13]*(x[3]*x[4]/x[11])+Ktm*x[1]*(x[3]-x[4])/x[11]
    xdot[5]=x[14]*(x[2]*x[1]/x[11])+Kp*x[1]*((x[3]+2*x[4])/x[11])-Kt*(x[3]*x[5]/x[11])+Ktm*x[1]*(x[3]-x[5])/x[11]
    xdot[6]=Ktm*x[1]*x[3]/x[11]+(x[15]+Ktc*0.5)*x[3]*x[3]/x[11]
    xdot[7]=Ktm*x[1]*x[4]/x[11]+x[13]*x[3]*x[4]/x[11]
    xdot[8]=Ktm*x[1]*x[5]/x[11]+x[13]*x[3]*x[5]/x[11]+Ktc*x[4]*x[4]/x[11]
    xdot[9]=Rlm-Rvm
    xdot[10]=Rlm
    xdot[11]=MWm*((1/rho_m)-(1/rho_p))*(-(x[14] + Ktm)*((x[3]*x[1])/(x[11])) + Rlm - Rvm - x[14]*((x[2]*x[1]/x[11])))
    xdot[12]=(-1/M0)*(-(x[14] + Ktm)*(x[3]*x[1]/x[11])+Rlm - Rvm - x[14]*(x[2]*x[1]/x[11]))
    xdot[13] =  -x[13]+(Ktd_0)*exp(A1 + A2*(x[12]) + A3*(x[12]**2) + A4*((x[12]**3)))
    xdot[14]=  -x[14] +(Kp_0)*exp(B1 + B2*(x[12]) + B3*((x[12])**2) + B4*((x[12])**3))
    xdot[15] =  -x[15]+ Ktd_0*exp(A1 + A2*x[12] + A3*(x[12]**2) + A4*(x[12]**3))
    return (xdot)
    

 #%%

temp=[]
sol=[]
Tt=[]
Xm=list()
rho_m = 966.5 - 1.1*(50)
rho_p = 1200;
V=((M0*(MWm))/(rho_m))
kt_0 = Ktdo_0*exp(-Etd/(Rg*(273.15+50)))
init_val=[I0, M0, 0.01, 0.01, 0, 0, 0, 0, 0, M0, M0, V, 0, (kt_0)*exp(A11*(50)+A12), (Kpo_0)*exp(B11*(50)+B12), (kt_0)*exp(A11*(50)+A12)]
I0=init_val
xvals=np.linspace(0, 150, 300)
for iter_n in range((xvals.shape[0])-1):
    #print (xvals[iter_n])
    t=xvals[iter_n]
    if t<120:
        c_p=50
    elif t<127.63:
        c_p=49.574+2.9406*(t-120)+0.5098*(t-120)**2-0.0733*(t-120)**3
    else:
        c_p=70
    c_p=50+273.15
    times = np.linspace(xvals[iter_n],xvals[iter_n+1],10)
    temp = odeint(PMMA, times, I0, c_p, integrator=0, method=1)
    I0 = temp[-1,:]
#    if iter_n==3:
    Xm.append(1-(I0[1]/I0[10]))
#    print(I0)
    sol.append(I0[4])
    Tt.append(times[0:-1])
#sol = np.concatenate(sol)
#Time = np.concatenate(T)
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('time (min)')
ax1.set_ylabel('Conversion', color=color)
#    ax1.plot(time, final_cb, color=color)
plt.plot(xvals[0:len(xvals)-1],Xm, color=color)
plt.show()  
