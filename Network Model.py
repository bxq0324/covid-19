#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


# In[2]:


df=pd.read_csv("Champaign.csv")
E=df.values[:,1:]
E


# In[3]:


X=pd.read_csv("subnet.csv").values
X[:,0][::-1]


# In[4]:


E=np.array([[0.5,0.5,0],[1/3,1/3,1/3],[0,1/2,1/2]])
A0=[15,6,5]
I0=[16,8,6]
R0=[126,108,47]
D0=[0,0,0]
N=[6145,4245,3844]
S0=[]
for i in range(3):
    S0.append(N[i]-A0[i]-I0[i]-R0[i])  
para=[0.05,0.5,0.1,0.08,0.17,0.0001]
T=16
S,A,I,R,D=S0,A0,I0,R0,D0
beta,q, gamma,kappa, rho,delta=para
print(S0)


# In[5]:


def SAIR_simulation(S0,A0,I0,R0,N,para,E,Data,district):
    S=S0
    A=A0
    I=I0
    R=R0
    data=[]  
    for t in range(T):
        b=[]
        a=0
        for i in range(3):
            for j in range(3):
                e=E[i][j]
                #print(A[j]+I[j])
                a+=e*(A[j]+I[j])/N[j]
            b.append(a)
        for i in range(3):
            #print(A)
            AA=A[i]
            II=I[i]
            SS=S[i]
            RR=R[i]
            #print(AA)
            #print(b)
            S[i]+=-beta * SS * b[i] +delta*RR
            #print(AA)
            A[i]+= q*beta * SS * b[i]- rho*AA-kappa*AA
            #print(A)
            #print(AA)
            I[i]+= (1-q)*beta * SS* b[i] +rho*AA-gamma*II
            R[i]+= kappa*AA+gamma * II-delta*RR
        d={"Susceptible":S[district],"Asymtomatic":A[district],"Infected":I[district],"Recovered":R[district]}
        data.append(d)
        N=[s+180 for s in N]
    df=pd.DataFrame(data)
    t = np.linspace(0, 16, 16)
    fig = plt.figure(figsize=(12,4))
    #plt.plot(t,df["Susceptible"])
    plt.plot(t,Data[:,district][::-1])
    plt.plot(t,df["Asymtomatic"])
    plt.plot(t,df["Infected"])
    plt.plot(t,df["Recovered"])
    plt.grid("True")
    plt.legend(["Original Data","Asymptonic","Infected","Recovered"])
    return df
            
SAIR_simulation(S0,A0,I0,R0,N,para,E,X,2)       


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




