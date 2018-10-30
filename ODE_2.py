# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 17:58:06 2018

@author: PRAVESH
"""

import numpy as np
tanh=np.tanh
import scipy.integrate
from scipy.integrate import odeint
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
I0=0.00258
M0=1E6
Rli=0
Rlm=0
Rvm=0
#%%
def PMMA(x, t, T):
    rho_m=966.5-1.1*(T-273.15)
    V1=(x[1]*MWm)/rho_m+(x[9]-x[1])*MWm/rho_p
    shy_m=((x[1]*MWm)/rho_m)/((x[1]*MWm)/rho_m+(x[9]-x[1])*MWm/rho_p)
    shy_p=1-shy_m
    Kd=Kd_0*np.exp(-Ed/(Rg*T))
    Kp_0=Kpo_0*np.exp(-Ep/(Rg*T))
    Ktd_0=Ktdo_0*exp(-Etd/(Rg*T))
    Mw=(MWm)*(x[5]+x[8])/(x[4]+x[7])
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
    xdot=[0 for var in range(0, 11)]
    xdot[0]=-Kd*x[0]+Rli
    xdot[1]=-(Kp+Ktm)*(x[1]*x[3]/V1)+Rlm-Rvm-Ki*(x[2]*x[1]/V1)
    xdot[2]=2*f*Kd*x[0]-Ki*(x[2]*x[1]/V1)
    xdot[3]=Ki*(x[2]*x[1]/V1)-Kt*((x[3]**2)/V1)
    xdot[4]=Ki*(x[2]*x[1]/V1)+Kp*x[1]*(x[3]/V1)-Kt*(x[3]*x[4]/V1)+Ktm*x[1]*(x[3]-x[4])/V1
    xdot[5]=Ki*(x[2]*x[1]/V1)+Kp*x[1]*((x[3]+2*x[4])/V1)-Kt*(x[3]*x[5]/V1)+Ktm*x[1]*(x[3]-x[5])/V1
    xdot[6]=Ktm*x[1]*x[3]/V1+(Ktd+Ktc*0.5)*x[3]*x[3]/V1
    xdot[7]=Ktm*x[1]*x[4]/V1+Kt*x[3]*x[4]/V1
    xdot[8]=Ktm*x[1]*x[5]/V1+Kt*x[3]*x[5]/V1+Ktc*x[4]*x[4]/V1
    xdot[9]=Rlm-Rvm
    xdot[10]=Rlm
    return (xdot)
    

#%%

temp=[]
sol=[]
T=[]
Xm=list()
init_val=[I0, M0, 0, 0, 0, 0, 0, 0, 0, M0, M0]
I0=init_val
for nTot in range(200,201,1):
    Xm=[]
    xvals=np.linspace(0, nTot, nTot)
    for iter_n in range((xvals.shape[0])-1):
        #print (xvals[iter_n])
        t=xvals[iter_n]
        if t<120:
            c_p=50
        elif t<127.63:
            c_p=49.574+2.9406*(t-120)+0.5098*(t-120)**2-0.0733*(t-120)**3
        else:
            c_p=70
        times = np.linspace(xvals[iter_n],xvals[iter_n+1],15)
        temp = odeint(PMMA, I0, times, args=(c_p,))
        I0 = temp[-1,:]
        Xm.append(1-I0[1]/I0[10])
        sol.append(I0[4])
        T.append(times[0:-1])
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('time (min)')
    ax1.set_ylabel('Conversion', color=color)
#    ax1.plot(time, final_cb, color=color)
    plt.plot(xvals[0:len(xvals)-1],Xm, color=color)
    plt.show()    
