#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.optimize import curve_fit


# In[ ]:





# In[2]:




def SAIR(y,t,beta,q, gamma,kappa, rho,delta):
    S, A, I, R = y
    dS = -beta * S * (A+I) +delta*R
    dA = q*beta * S * (A+I) - rho*A-kappa*A
    dI = (1-q)*beta * S * (A+I) +rho*A-gamma*I
    dR = kappa*A+gamma * I-delta*R
    return dS, dA, dI, dR


# In[ ]:





# In[3]:


data=[137,134,117,114,110,114,117,99,100,94,98]
l=[]
ntested=31600
for i in range(len(data)):
    l.append(ntested)
    ntested+=900
print(l)
for i in range(len(data)):
    data[i]=data[i]/l[i]
print(data)    

beta = 0.05  # infected person infects 1 other person per day
q=0.67
gamma = 0.1
kappa=0.08
rho=0.17
delta = 0.0001 # incubation period of three days

A0=0.0020
I0=0.0031
R0=0.0231
S0=1-A0-I0-R0
y0=S0,A0,I0,R0
t = np.linspace(0, 11, 11)
sol = odeint(SAIR, y0, t, args=(beta,q, gamma,kappa, rho,delta))
print(sol.T)

fig = plt.figure(figsize=(12,4))
plt.plot(t,sol.T[1])
plt.plot(t,sol.T[2])
plt.plot(t,sol.T[3])
plt.plot(np.arange(11),data,"k*:")
plt.grid("True")
plt.legend(["Asymptonic","Infected","Recovered","Original Data"])


error=[]
for i in range (len(data)):
    error.append((sol.T[2][i]-data[i])/data[i])
print(error)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




