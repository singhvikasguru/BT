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
import torch
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
#d_ca=2.5
#d_cb=2.5
S_I0=np.linspace(low_I0, up_I0, 51)
#CA=np.logspace(3.301029995, 3.579783596,num=100, base=10.0)
#CA=np.logspace(3.30102999,3.5797358, num=101, base=10.0)
S_M0=np.linspace(low_M, up_M, 51)
states=list(itertools.product(S_I0,S_M0))
moves=np.linspace(10,10000,51)

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
def odeint(fun, t, y0, u, integrator=0, method=0, rtol=0, atol=1e-12):
    from scipy.integrate import ode
    s = ode(fun).set_f_params(u)
    integrators = ['vode', 'zvode', 'lsoda', 'dopri5', 'dop853']
    methods = ['adams', 'bdf']
    s.set_integrator(integrators[0],
                     method=methods[0],
                     order=10,
                     rtol=rtol,
                     atol=atol,
                     with_jacobian=False)
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
    Cb_ref=1
    init_t=time
    u=moves[action_a]
    times=np.linspace(init_t, init_t+dt, 2)

    temp = odeint(func, times, Io, u, integrator=0, method=1)
#    temp = odeint(CSTR, Io, times, args=(u,))
    All = temp[-1,:]
    new_Cb=All[1]
#    if time>=4:
#        if new_Cb >= .9 and new_Cb <= 1.1:
#            reward = 100000
#        else:
#            reward = -np.inf
#    else:
#        reward = 0
    del_X = np.abs(new_Cb-Cb_ref)
    if np.abs(del_X) <0.1:
        reward=(np.abs(del_X))
    else:
        reward = -1

    return All, reward

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def make_policy(epsilon, leng):
    nA=leng
    def policy_fn(observation):
        A = np.ones(nA, dtype=float) * epsilon / nA
        best_action = np.argmax(predict(observation))
        A[best_action] += (1.0 - epsilon)
        return A
    return policy_fn
#%%
def plot_G(act_ep, ca_ep, cb_ep, tc_ep, n):
    final_act = act_ep[-1,:]
    final_ca = ca_ep[-1,:]
    final_cb = cb_ep[-1,:]
    final_tc = tc_ep[-1,:]

    final_act = final_act[0:int(T/dt)]
    final_ca= final_ca[0:int(T/dt)]
    final_cb= final_cb[0:int(T/dt)]
    final_tc= final_tc[0:int(T/dt)]



    fig, ax1 = plt.subplots()
    time = np.linspace(0,T,int(T/dt))
    color = 'tab:red'
    ax1.set_xlabel('time (min)')
    ax1.set_ylabel('C_b', color=color)
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
    plt.savefig('n')


#%%
Q=np.zeros([len(states), len(moves)])
for i in range(0, len(states)):
    for j in range(0, len(moves)):
        state=i
        next_action=j
        All, reward = get_state(state, next_action, 0, [2, 1.117])
        Q[i][j]=reward
#%%

device = torch.device('cpu')
# N is batch size; D_in is input dimension;
# H is hidden dimension; D_out is output dimension
N, D_in, H, D_out = 20, 3, 5, 1
model = torch.nn.Sequential(
              torch.nn.Linear(D_in, H),
              torch.nn.ReLU(),
              torch.nn.ReLU(),
              torch.nn.ReLU(),
              torch.nn.Linear(H, D_out),
            ).to(device)
def nn_model(x, y):
    # device = torch.device('cuda') # Uncomment this to run on GPU
  
    # Create random Tensors to hold inputs and outputs
#    x2=np.array([[40, 13.2, 13],[61.1,41,30], [10,5,1.0]], dtype=np.float32)
#    y=np.array([[40],[13],[40.1]], dtype=np.float32)
    x=torch.from_numpy(x)
    y=torch.from_numpy(y)
    x=x.float()
    y=y.float()
    #x = torch.randn(N, D_in, device=device)
    #y = torch.randn(N, D_out, device=device)
#    print('-----------------x,y----------------')
#    print(x)
#    print(y)
    # Use the nn package to define our model as a sequence of layers. nn.Sequential
    # is a Module which contains other Modules, and applies them in sequence to
    # produce its output. Each Linear Module computes output from input using a
    # linear function, and holds internal Tensors for its weight and bias.
    # After constructing the model we use the .to() method to move it to the
    # desired device.
    
    
    # The nn package also contains definitions of popular loss functions; in this
    # case we will use Mean Squared Error (MSE) as our loss function.
    loss_fn = torch.nn.MSELoss(size_average=False)
    
    learning_rate = 1e-4
    for t in range(100):
      # Forward pass: compute predicted y by passing x to the model. Module objects
      # override the __call__ operator so you can call them like functions. When
      # doing so you pass a Tensor of input data to the Module and it produces
      # a Tensor of output data.
      y_pred = model(x)
    
      # Compute and print loss. We pass Tensors containing the predicted and true
      # values of y, and the loss function returns a Tensor containing the loss.
      loss = loss_fn(y_pred, y)
#      print(t, loss.item())
      
      # Zero the gradients before running the backward pass.
      model.zero_grad()
    
      # Backward pass: compute gradient of the loss with respect to all the learnable
      # parameters of the model. Internally, the parameters of each Module are stored
      # in Tensors with requires_grad=True, so this call will compute gradients for
      # all learnable parameters in the model.
      loss.backward()
    
      # Update the weights using gradient descent. Each parameter is a Tensor, so
      # we can access its data and gradients like we did before.
      with torch.no_grad():
        for param in model.parameters():
          param.data -= learning_rate * param.grad
    
    
#    print(model(x))     

#%%

def predict(state_s):
    Mat=np.zeros([1, 3])
    state_s=state_s.tolist()
    state_s2=list()
    for i in range(0, len(moves)):
        state_s2=state_s[:]
        action_a=moves[i]
        state_s2.append(action_a)
        state_s2=np.asarray(state_s2)
        Mat=np.vstack([Mat, state_s2])
    Mat=np.delete(Mat, (0), axis=0)
    Mat=torch.from_numpy(Mat)
    Mat=Mat.float()
    res=model(Mat)
    res=res.data.numpy()
    return(res)
#%%
r_memory=list()
m_limit=20000
discount=0.6
def Save(All):
    if len(r_memory)>m_limit:
        del r_memory[0]
    r_memory.append(All)

def get_data(batches):
    inputs=list()
    target=list()
    len_memory=len(r_memory)
    if len_memory<batches:
        for i in range(0, len_memory):
            [state_s, action_a, reward_r,state_ss]=r_memory[i]
            init_action=action_a
            #final_action=action_a[1]
            state_s2=list()
            state_s2=state_s[:]
            state_s2.append(init_action)
            #state_s2.append(final_action)
            inputs.append(state_s2)
            Q_sa=np.max(predict(state_ss))
#            print('--------------here----------------')
#            print(Q_sa)
#            time.sleep(1)
            target.append([reward_r+discount*Q_sa])
    else:
        for i, value in enumerate(np.random.randint(0, len_memory, size=batches)):
            [state_s, action_a, reward_r,state_ss]=r_memory[i]
            init_action=action_a
            #final_action=action_a[1]
            state_s2=list()
            state_s2=state_s[:]
            state_s2.append(init_action)
            #state_s2.append(final_action)
            inputs.append(state_s2)
            Q_sa=np.max(predict(state_ss))
#            print(Q_sa)
            target.append([reward_r+discount*Q_sa])
    return (np.asarray(inputs),np.asarray(target)) 

#%%

batches=50
def q_mat(N_episodes, gamma=1, alpha=0.4):
    #Q = defaultdict(lambda: np.zeros(len(moves)))

    epsilon=0.1
    temp_lst=list()
    act_ep=np.zeros([N_episodes, int(T/dt)+1])
    ca_ep=np.zeros([N_episodes, int(T/dt)+1])
    cb_ep=np.zeros([N_episodes, int(T/dt)+1])
    tc_ep=np.zeros([N_episodes, int(T/dt)+1])
    rew_ep=np.zeros([N_episodes, int(T/dt)+1])
    ind_ep=np.zeros([N_episodes, int(T/dt)+1])
#    Num_ep=np.asarray([1, 2, 3, 4, 5, 6, 30, 50, 100, 200, 400, 700, 900, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000])
    Num_ep=np.linspace(1,N_episodes, N_episodes)
    for n in range(N_episodes):
        if n in Num_ep:
            print('No. of episodes = ', n)
            plot_G(act_ep[0:n], ca_ep[0:n], cb_ep[0:n], tc_ep[0:n], n)
        policy = make_policy(epsilon, len(moves))
        j=0
        t=0
        u_old=0
        u_new=0
#        ca = find_nearest(CA, 2)
#        cb = find_nearest(CB, 1.117)
##        tc = find_nearest(TC, 90)
#        state=states.index((ca, cb))
        temp_lst=list()
        temp_lst.append([2, 1.117])
        while t<T:
            time.sleep(2)
            Io=[temp_lst[-1][0], temp_lst[-1][1]]
            action_probs = policy(np.asarray(Io))
            next_action = np.random.choice(np.arange(len(action_probs)), p=action_probs)
            act_ep[n][j]=moves[next_action]
            ind_ep[n][j]=state
            ca, cb=states[state]
            ca_ep[n][j]=Io[0]
            cb_ep[n][j]=Io[1]
            new_state, reward = get_state(state, next_action, t, Io)
           
            t+=dt
            temp_lst.append(new_state)
            u_new=moves[next_action]
            del_u=np.abs(u_new-u_old)
            u_old=u_new
#            if del_u>=0.40:
#                if reward>=0:
#                    reward/=10
#                else:
#                    reward=-1

            All=[Io, moves[next_action], reward, new_state]
            Save(All)
            input_t, target=get_data(batches)
            nn_model(input_t, target)
            rew_ep[n][j]=reward



#            if new_ca <= 21000 or new_ca >= 3700:
#                reward = -np.inf
#            elif new_cb <= 1000 or new_cb >= 1100:
#                reward = -np.inf
#            elif new_tc <= 80 or new_tc >= 150:
#                reward = -np.inf
#            else:
#                reward = reward

#            next_state=states.index((new_ca, new_cb))
            Val = np.max(Q[next_state])                     # Maximum Q-value of the states accessible from the next state
            Q[state][next_action] = (1-alpha)*Q[state][next_action] + alpha*(reward + gamma*Val)      # Update Q-values
            state=next_state
            j+=1

    return Q, act_ep, ca_ep, cb_ep, tc_ep, rew_ep, temp_lst, ind_ep


N_episodes=15000
Q, act_ep, ca_ep, cb_ep, tc_ep, rew_ep, All, ind_ep = q_mat(N_episodes)
All=np.asarray(All)
for i in range(0, cb_ep.shape[0]):
    for j in range(0, cb_ep.shape[1]):
        sheet1.write(i, j, cb_ep[i][j])
        sheet2.write(i, j, act_ep[i][j])
        sheet3.write(i, j, ca_ep[i][j])
        sheet4.write(i, j, rew_ep[i][j])
for i in range(0, All.shape[0]):
    sheet5.write(i, 0, All[i][0])
    sheet5.write(i, 1, All[i][1])

print(N_episodes)
plot_G(act_ep, ca_ep, cb_ep, tc_ep, N_episodes)
book.save("CSTR_DRL_withoutCAS_point1_5K_reward_1bydelta_or_-50k_U_Qby10_and_inf_epsilol_point1.xls")
   
