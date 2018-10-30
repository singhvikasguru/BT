# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 11:44:15 2018

@author: 1544540
"""
import numpy as np
tanh=np.tanh
import matplotlib.pyplot as plt
import math
pow=math.pow
exp=np.exp
import itertools
from collections import defaultdict
import time
import xlwt
import matlab.engine
eng=matlab.engine.start_matlab()

book = xlwt.Workbook(encoding="utf-8")
sheet1 = book.add_sheet("M0")
sheet2 = book.add_sheet("I0")
sheet3 = book.add_sheet("Rlm")
sheet4 = book.add_sheet("Reward")
sheet5 = book.add_sheet("ODE")


low_I0=0
up_I0=.01
low_M=0
up_M=2000000

T=120
global dt
dt=1
S_I0=np.linspace(low_I0, up_I0, 51)
S_M0=np.linspace(low_M, up_M, 51)
states=list(itertools.product(S_I0,S_M0))
moves=np.linspace(10,10000,51)
#CA=np.logspace(3.301029995, 3.579783596,num=100, base=10.0)
#CA=np.logspace(3.30102999,3.5797358, num=101, base=10.0)

#%%

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
M0=5000000
Rli=0.00000000
Rlm=1000
Rvm=0

#%%
def PMMA(t, x, u, T):
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
    xdot[1]=-(x[14]+Ktm)*(x[1]*x[3]/x[11])+u-Rvm-Ki*(x[2]*x[1]/x[11])
    xdot[2]=2*f*Kd*x[0]-Ki*(x[2]*x[1]/x[11])
    xdot[3]=x[14]*(x[2]*x[1]/x[11])-Kt*(x[3]**2/x[11])
    xdot[4]=x[14]*(x[2]*x[1]/x[11])+Kp*x[1]*(x[3]/x[11])-x[13]*(x[3]*x[4]/x[11])+Ktm*x[1]*(x[3]-x[4])/x[11]
    xdot[5]=x[14]*(x[2]*x[1]/x[11])+Kp*x[1]*((x[3]+2*x[4])/x[11])-Kt*(x[3]*x[5]/x[11])+Ktm*x[1]*(x[3]-x[5])/x[11]
    xdot[6]=Ktm*x[1]*x[3]/x[11]+(x[15]+Ktc*0.5)*x[3]*x[3]/x[11]
    xdot[7]=Ktm*x[1]*x[4]/x[11]+x[13]*x[3]*x[4]/x[11]
    xdot[8]=Ktm*x[1]*x[5]/x[11]+x[13]*x[3]*x[5]/x[11]+Ktc*x[4]*x[4]/x[11]
    xdot[9]=u-Rvm
    xdot[10]=u
    xdot[11]=MWm*((1/rho_m)-(1/rho_p))*(-(x[14] + Ktm)*((x[3]*x[1])/(x[11])) + u - Rvm - x[14]*((x[2]*x[1]/x[11])))
    xdot[12]=(-1/M0)*(-(x[14] + Ktm)*(x[3]*x[1]/x[11])+u - Rvm - x[14]*(x[2]*x[1]/x[11]))
    xdot[13] =  -x[13]+(Ktd_0)*exp(A1 + A2*(x[12]) + A3*(x[12]**2) + A4*((x[12]**3)))
    xdot[14]=  -x[14] +(Kp_0)*exp(B1 + B2*(x[12]) + B3*((x[12])**2) + B4*((x[12])**3))
    xdot[15] =  -x[15]+ Ktd_0*exp(A1 + A2*x[12] + A3*(x[12]**2) + A4*(x[12]**3))
    return (xdot)

#%%
def odeint(fun, t, y0, u, T, integrator=0, method=0, rtol=0, atol=1e-12):
    from scipy.integrate import ode
    s = ode(fun).set_f_params(u, T)
    integrators = ['vode', 'zvode', 'lsoda', 'dopri5', 'dop853']
    methods = ['adams', 'bdf']
    s.set_integrator(integrators[2],
                     method=methods[1],
#                     order=15)
                     nsteps=30,
                     rtol=0.001)
#                     atol=atol,
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

def get_state(state, action_a, time, Io):
    if (Io[4]+Io[7])==0:
        Mw_avg_old=0.10013
    else:
        Mw_avg_old=MWm*(Io[5]+Io[8])/(Io[4]+Io[7])
#    Cb_ref=1
    init_t=time
    u=moves[action_a]
#    times=np.linspace(init_t, init_t+dt, 10)
    n=2
#    T=50+273.5
#    temp = odeint(PMMA, times, Io, u, T, integrator=0, method=1)
    ti, temp2=eng.MMA_Simulation(float(init_t+dt), float(init_t), float(n), float(u), float(Io[0]), float(Io[1]), float(Io[2]), float(Io[3]), float(Io[4]), float(Io[5]), float(Io[6]), float(Io[7]), float(Io[8]), float(Io[9]), float(Io[10]), float(Io[11]), float(Io[12]), float(Io[13]), float(Io[14]), float(Io[15]), 1500000.0, 0.0258, nargout=2)
    All = temp2[-1]
    Mw_avg_new=MWm*(All[5]+All[8])/(All[4]+All[7])
    reward=(Mw_avg_new-Mw_avg_old)*1000000000000

    return All, reward

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def make_policy(Q, epsilon, leng):
    nA=leng
    def policy_fn(observation):
        A = np.ones(nA, dtype=float) * epsilon / nA
        best_action = np.argmax(Q[observation])
        A[best_action] += (1.0 - epsilon)
        return A
    return policy_fn
#%%
def plot_G(cb_ep, ca_ep, act_ep, All):
    final_act = act_ep[-1,:]
    final_ca = ca_ep[-1,:]
    final_cb = cb_ep[-1,:]
#    final_tc = tc_ep[-1,:]
    
    final_act = final_act[0:int(T/dt)]
    final_ca= final_ca[0:int(T/dt)]
    final_cb= final_cb[0:int(T/dt)]
#    final_tc= final_tc[0:80]

#    if len(All)>=80:
#        
#        x1 = ca_ep[-1,:]
#        x1 = x1[0:80]
#        All1 = np.vstack(All)
#        x2 = All1[:,0]
#        x2 = x2[0:80]
#        plt.plot(x2,x1,'ro')
#        c = x2**2 + x1**2
#    
#        fig, ax = plt.subplots()
#        ax.scatter(x2, x1, s=25, c=c, cmap=plt.cm.coolwarm, zorder=10)
#        lims = [
#                np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#                np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
#                ]
#    
#    # now plot both limits against eachother
#        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
#        ax.set_aspect('equal')
#        ax.set_xlim(lims)
#        ax.set_ylim(lims)
#        fig.savefig('so.png', dpi=300)
#        
#        
#        if len(All)>=80:
#            x1 = cb_ep[-1,:]
#            x1 = x1[0:80]
#            All1 = np.vstack(All)
#            x2 = All1[:,0]
#            x2 = x2[0:80]
#            plt.plot(x2,x1,'ro')
#            c = x2**2 + x1**2
#        
#            fig, ax = plt.subplots()
#            ax.scatter(x2, x1, s=25, c=c, cmap=plt.cm.coolwarm, zorder=10)
#            lims = [
#                    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#                    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
#                    ]
#        
#        # now plot both limits against eachother
#            ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
#            ax.set_aspect('equal')
#            ax.set_xlim(lims)
#            ax.set_ylim(lims)
#            fig.savefig('so.png', dpi=300)



    fig, ax1 = plt.subplots()
    time = np.linspace(0,T,int(T/dt))
    color = 'tab:red'
    ax1.set_xlabel('time (min)')
    ax1.set_ylabel('Avg. Molecular weight', color=color)
    ax1.plot(time, final_cb, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'tab:blue'
    ax2.set_ylabel('Flowrate', color=color)  # we already handled the x-label with ax1
    #ax2.plot(time, final_act, color=color)
    ax2.step(time, final_act,  where='post', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.grid(color='g', linestyle='-', linewidth=1)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
#%%
    
u = np.linspace(0,1,10)
sol = list()
rho_m = 966.5 - 1.1*(50)
V=((M0*(MWm))/(rho_m))
kt_0 = Ktdo_0*exp(-Etd/(Rg*(273.15+50)))
times=np.linspace(0, 5, 50000)
init_val=[I0, M0, 0.01, 0.01, 0, 0, 0, 0, 0, M0, M0, V, 0, (kt_0)*exp(A11*(50)+A12), (Kpo_0)*exp(B11*(50)+B12), (kt_0)*exp(A11*(50)+A12)]
Io=init_val
#for i in range(len(u)):
#    temp = odeint(PMMA, times, Io, u[i], integrator=0, method=1)
#    All = temp[:,0]
#    sol.append(All)
#    sol1 = np.vstack(sol)
#fig, ax1 = plt.subplots()
#for i in range(sol1.shape[0]):
#    ax1.plot(times, sol1[i][0:50000], label= i)
#    print( max(sol1[i]))
#plt.legend()
#plt.show()
#%%
Q=np.zeros([len(states), len(moves)])
#for i in range(0, len(states)):
#    for j in range(0, len(moves)):
#        state=i
#        next_action=j
#        All, reward = get_state(state, next_action, 0, Io) 
#        Q[i][j]=reward    
#%%
    
def q_mat(N_episodes, gamma=0.99, alpha=0.45):
    init_val=[I0, M0, 0.01, 0.01, 0, 0, 0, 0, 0, M0, M0, V, 0, (kt_0)*exp(A11*(50)+A12), (Kpo_0)*exp(B11*(50)+B12), (kt_0)*exp(A11*(50)+A12)]
    Io=init_val
    epsilon=0.1
    temp_lst=list()
    M0_ep=np.zeros([N_episodes, int(T/dt)+1])
    I0_ep=np.zeros([N_episodes, int(T/dt)+1])
    Rlm_ep=np.zeros([N_episodes, int(T/dt)+1])
    rew_ep=np.zeros([N_episodes, int(T/dt)+1])
    ind_ep=np.zeros([N_episodes, int(T/dt)+1])
#    Num_ep=np.asarray([1, 10, 30, 50, 100, 200, 400, 700, 900, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000])
    Num_ep=np.linspace(1,150,150)
    Num_ep=[int(x) for x in Num_ep]
    for n in range(N_episodes):
        temp_lst=list()
        temp_lst.append(Io)   
        if n in Num_ep:
            print('No. of episodes = ', n)
            plot_G(M0_ep[0:n], I0_ep[0:n], Rlm_ep[0:n], temp_lst[0:n])
        policy = make_policy(Q, epsilon, len(moves))
        j=0
        t=0
        i0 = find_nearest(S_I0, I0)
        m0 = find_nearest(S_M0, M0)
        state=states.index((i0, m0))
        print(M0_ep[n-1][int(t/dt)-1])
        while t<T:
            action_probs = policy(state)
            next_action = np.random.choice(np.arange(len(action_probs)), p=action_probs)
            Rlm_ep[n][j]=moves[next_action]
            ind_ep[n][j]=state
            i0, m0=states[state]
            I0_ep[n][j]=i0
#            Io=[temp_lst[-1][0], temp_lst[-1][1]]
            All, reward = get_state(state, next_action, t, Io)
            t+=dt
            Io=All
            M0_ep[n][j]=MWm*(Io[5]+Io[8])/(Io[4]+Io[7])
            temp_lst.append(All)
        
            new_I0=All[0]
            new_M0=All[1]

            rew_ep[n][j]=reward
                
            new_I0=find_nearest(S_I0, new_I0)
            new_M0=find_nearest(S_M0, new_M0)
            
            next_state=states.index((new_I0, new_M0))
            Val = np.max(Q[next_state])                     # Maximum Q-value of the states accessible from the next state
            Q[state][next_action] = (1-alpha)*Q[state][next_action] + alpha*(reward + gamma*Val)      # Update Q-values
            state=next_state
            j+=1
    return Q, M0_ep, I0_ep, Rlm_ep, rew_ep, temp_lst, ind_ep


N_epiosdes=600
Q, M0_ep, I0_ep, Rlm_ep, rew_ep, All, ind_ep = q_mat(N_epiosdes)
#All=np.asarray(All)
#for i in range(0, cb_ep.shape[0]):
#    for j in range(0, cb_ep.shape[1]):
#        sheet1.write(i, j, cb_ep[i][j])
#        sheet2.write(i, j, act_ep[i][j])
#        sheet3.write(i, j, ca_ep[i][j])
#        sheet4.write(i, j, rew_ep[i][j])
#for i in range(0, All.shape[0]):
#    sheet5.write(i, 0, All[i][0])
#    sheet5.write(i, 1, All[i][1])
#    
#print(N_epiosdes)
#plot_G(act_ep, ca_ep, cb_ep, tc_ep, All)
#book.save("CSTR_RL3_20_20_20_1_point1_3K_reward_1bydelta_or_-50k_U_Qby10_and_inf_epsilol_point1.xls")
# 