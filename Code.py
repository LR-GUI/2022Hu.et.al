"""
Created on Tue Feb 28 16:40:48 2023

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint




def LotkaVolterra(y0,a,D,time):
# Function defining the model. Inputs initial conditions, parameters and maxtime, outputs time and trajectories arrays
    
    def System(y,t,a,D): # Function defining the system of differential equations
        
        
        dy = []
        for i in range(len(y)):
            A = 0
            diff = 0

            for j in range(len(y)):
                A += y[j]*a[i][j]
                
            diff = y[i]*(1-A)+D

            dy.append(diff)
        

        return dy
    
    

    T = np.linspace(0,time,100*time+1) # Time array. We integrate the function for these time values
    sol = odeint(System, y0, T, args=(a,D,)) # Call odeint to integrate the system

    return sol,T






# For plots 1.A and 1.B

S = 50
D = 0.000001
am = 0.64
maxtime = 10000
y0 = [1/S]*S
a = np.zeros((S,S))
for i in range(S):
    for j in range(S):
        if(i==j):
            a[i][j] = 1
        else:
            a[i][j] = 2*am*np.random.rand()




sol,solT = LotkaVolterra(y0,a,D,maxtime)

ss = np.zeros(S)
ss = sol.T

for n in ss:    
    plt.plot(solT,n)

plt.xscale('log')
plt.yscale('log')
plt.ylim(0.000001,10)
plt.xlim(1,10000)

plt.ylabel('Densities')
plt.xlabel('Time')
plt.title('S='+str(S)+' and '+'<a>='+str(am))

plt.savefig('S'+str(S)+'a'+str(am)+'.png', dpi=300,format='png')




# For plot 1.C

S = 50
D = 0.000001
am = np.linspace(0.01,1,50)
y0 = [1/S]*S
maxtime = 200
mean_survival = []
for am_idx in range(len(am)):
    
    survival = []
    for iter in range(100):
        a = np.zeros((S,S))
        for i in range(S):
            for j in range(S):
                if(i==j):
                    a[i][j] = 1
                else:
                    a[i][j] = 2*am[am_idx]*np.random.rand()
    
    
    
    

        sol,solT = LotkaVolterra(y0,a,D,maxtime)
    
        ss = np.zeros(S)
        ss = sol.T
    
        survival_species = 0
        for n in ss:
            
            amount = np.mean(n[-100:-1])
            if amount>0.001:
                survival_species += 1
        survival_species = survival_species/S
        survival.append(survival_species)
    mean_survival.append(np.mean(survival))
plt.plot(am,mean_survival)
plt.ylabel('Survival fraction')
plt.xlabel('Interaction strength')
plt.title('Survival. S='+str(S))
plt.savefig('Survival.png', dpi=300,format='png')



# For plot 1.D

S = 50
D = 0.000001
am = np.linspace(0.01,1,25)
y0 = [1/S]*S
maxtime = 1000
n_iter=100
mean_communities = []
for am_idx in range(len(am)):
    
    fluctuating_communities = 0
    for iter in range(n_iter):
        a = np.zeros((S,S))
        for i in range(S):
            for j in range(S):
                if(i==j):
                    a[i][j] = 1
                else:
                    a[i][j] = 2*am[am_idx]*np.random.rand()
    
        
    
    
        sol,solT = LotkaVolterra(y0,a,D,maxtime)
    
        ss = np.zeros(S)
        ss = sol.T
    
        fluctuating_species = []
        for n in ss:

            fluctuating_species.append(np.std(n[-500:-1])/np.mean(n[-500:-1]))
            
        amount = np.mean(fluctuating_species)
        if amount>0.001:
            fluctuating_communities += 1
    mean_communities.append(fluctuating_communities/n_iter)
plt.plot(am,mean_communities)
plt.ylabel('Fluctuation fraction')
plt.xlabel('Interaction strength')
plt.title('Fluctuation. S='+str(S))
plt.savefig('Fluctuation.png', dpi=300,format='png')



# For plot in 1.E

D = 0.000001
S = np.arange(0,110,10)
S[0] = 2
am = np.arange(0,1.1,0.1)
am[0] = 0.01

maxtime = 200

solSurvival = [[0]*len(S) for i in range(len(am))]

for i in range(len(S)):
  for j in range(len(am)):
      
    y0 = [1/S[i]]*S[i]
    survival = []
    for iter in range(100):
        a = np.zeros((S[i],S[i]))
        for ii in range(S[i]):
            for jj in range(S[i]):
                if(ii==jj):
                    a[ii][jj] = 1
                else:
                    a[ii][jj] = 2*am[j]*np.random.rand()
    
    
    
    

        sol,solT = LotkaVolterra(y0,a,D,maxtime)
    
        ss = np.zeros(S[i])
        ss = sol.T
    
        survival_species = 0
        for n in ss:
            
            amount = np.mean(n[-100:-1])
            if amount>0.001:
                survival_species += 1
        survival_species = survival_species/S[i]
        survival.append(survival_species)
    solSurvival[j][i] = np.mean(survival)
      
    
    
plt.contourf(S,am,solSurvival,50,vmin=0,vmax=1,cmap=plt.cm.viridis)
plt.colorbar()
plt.ylabel('Interaction strength, <a>')
plt.xlabel('Size of species pool, S')
plt.title('Fraction of surviving species')
plt.savefig('SurvivalSpace.png', dpi=300,format='png')




# For plot in 1.F

D = 0.000001
S = np.arange(0,110,20)
S[0] = 2
am = np.arange(0,1.1,0.2)
am[0] = 0.01
n_iter = 100
maxtime = 200

solFluctuation = [[0]*len(S) for i in range(len(am))]

for i in range(len(S)):
  for j in range(len(am)):
      
    y0 = [1/S[i]]*S[i]
    fluctuating_communities = 0
    for iter in range(n_iter):
        a = np.zeros((S[i],S[i]))
        for ii in range(S[i]):
            for jj in range(S[i]):
                if(ii==jj):
                    a[ii][jj] = 1
                else:
                    a[ii][jj] = 2*am[j]*np.random.rand()
    
    
    
    

        sol,solT = LotkaVolterra(y0,a,D,maxtime)
    
        ss = np.zeros(S[i])
        ss = sol.T
    
        fluctuating_species = []
        for n in ss:

            fluctuating_species.append(np.std(n[-500:-1])/np.mean(n[-500:-1]))
            
        amount = np.mean(fluctuating_species)
        if amount>0.001:
            fluctuating_communities += 1
    solFluctuation[j][i] = fluctuating_communities/n_iter
      
    
    
plt.contourf(S,am,solFluctuation,50,vmin=0,vmax=1,cmap=plt.cm.viridis)
plt.colorbar()
plt.ylabel('Interaction strength, <a>')
plt.xlabel('Size of species pool, S')
plt.title('Fraction of fluctuating communities')
plt.savefig('FluctuationSpace.png', dpi=300,format='png')