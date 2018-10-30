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

c = genfromtxt('c_profile.csv', delimiter=',')
c_x=c[:, 0]
c_y=c[:, 1]
m_dot_in_max=0.00756
t_in_M_0=30*60
t_in_M_1=100*60
im=0.9 # impurity factor in the range 0.8 to 1.2, value taken randomly by myself
k0=55
k1=1000
k2=0.4
E=29560.89
R=8.314
a0=555.56
c0=5.2*pow(10, -5)
c1=16.4
c2=2.3
c3=1.563
C_p_M=1.675
C_p_P=3.140
C_p_W=4.187
C_p_C=4.187
T_amb=280.382 # assumed winter, for summer it is 305.382
rho_m=900
rho_p=1040
rho_w=1000 
m_W=205.5
m_C=42.750
UA_loss=0.00567567
d0=0.814
d1=-5.13
hf_inv=0.0 #using the value of 2nd batch
mc_dot=0.9412
theta_1=22.8
theta_2=15
tau_p=40.2
T_inlet=294.26
T_steam=449.82
P=7.693
B1=4.6659
B2=4.7427
delta_H=70152.16
#%%
def PMMA(x, t, c):
    if t>=t_in_M_0 and t<=t_in_M_1:
        m_dot_M=m_dot_in_max
    else:
        m_dot_M=0
    print('----T_wall--------')
    print(t)
    print(x)
    f=x[1]/(x[0]+x[1]+m_C)
    print(f)
    myu=c0*exp(c1*f)*pow(10, c2*(a0/x[2]-c3))
    print(myu)
    k=k0*exp(-E/(R*x[2]))*pow(k1*myu, k2)
    print(k)
    R_p=im*k*x[0]
    A=(x[0]/rho_m+x[1]/rho_p+m_W/rho_w)*(P/B1)+B2
    Tj=(x[3]+x[4])/2
    T_wall=(x[2]+Tj)/2
    myu_wall=c0*exp(c1*f)*pow(10, c2*(a0/T_wall-c3))
    h=d0*exp(d1*myu_wall)
    U=1/(1/h+ hf_inv)
    if c<0.5:
        Kp=0.8*pow(30, -c/50)*(T_inlet-x[4])
    elif c==0.5:
        Kp=0
    else:
       Kp=0.15*pow(30, c/50-2)*(T_steam-x[4])
    mcp=x[0]*C_p_M+x[1]*C_p_P+m_W*C_p_W
    xdot=[0 for var in range(0,5)]
    xdot[0]=m_dot_M-R_p
    xdot[1]=R_p
    xdot[2]=(1/(mcp))*(m_dot_M*C_p_M*(T_amb-x[2])-U*A*(x[2]-Tj)-UA_loss*(x[2]-T_amb)-delta_H*R_p)
    xdot[3]=(1/(m_C*C_p_C))*(mc_dot*C_p_C*(x[4]*(t-theta_1))-x[3])+U*A*(x[2]-Tj)
    xdot[4]=xdot[3]*(t-theta_2)+(x[3]*(t-theta_2)-x[4])/tau_p+Kp/tau_p
    return (xdot)

#%%

temp=[]
sol=[]
T=[]
init_val=[0, 53.98, T_amb, T_amb, T_amb]
I0=init_val
xvals=np.linspace(0, 4000, 201)
for iter_n in range((xvals.shape[0])-1):
    #print (xvals[iter_n])
    t=xvals[iter_n]
    for i in range(0, len(c_x)):
        if c_x[i]>t/60:
            ind=i
            break
    a=c_y[i]
    b=c_y[i-1]
    x2=c_x[i]
    x1=c_x[i-1]
    c_p=a*(t/60-x1)/(x2-x1)+b*(x2-t/60)/(x2-x1)
    times = np.linspace(xvals[iter_n],xvals[iter_n+1],3)
    temp = odeint(PMMA, I0, times, args=(c_p,))
    I0 = temp[-1,:]
#    if I0[4]>368:
#        I0[4]=368
#    if I0[4]<278:
#        I0[4]=278
    print(I0)
    sol.append(I0[4])
    T.append(times[0:-1])
#sol = np.concatenate(sol)
#Time = np.concatenate(T)

plt.plot(Time,(sol[:,4]),'b')
plt.show()    
