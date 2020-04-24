# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a Heat Transfer script file.
"""

import numpy as np
from matplotlib import pyplot

''' Geometry '''
s_B  = 0.20

''' Spacial grid '''
dx   =  0.01
nx   = int(s_B/dx) + 1
x    = np.linspace(0,s_B,nx)

'''Time grid'''
t_sim = 14400
dt   =  0.5
nt   = int(t_sim/dt)

'''Material Property'''
rho   = 2400
cp    = 1000
lamda = 1.13

'''Initial Condition'''
T_A   = 303.15

'''Define Results Parameter'''
T    = np.ones(nx)*T_A

'''simulation standard fire curve'''
'''Initial Condition'''
T_A   = 303.15

'''Define Results Parameter'''
T    = np.ones(nx)*T_A

for n in range(1,nt):
        Tn = T.copy()
        T[1:-1] = Tn[1:-1] + dt * lamda/(rho*cp) * (Tn[2:] - 2*Tn[1:-1] +Tn[0:-2])/dx**2
        T[0] = 293.15 + 345*np.log10(8*n*dt/60+1)
        T[-1] = T[-2]
print (T-273.15)
plotT = pyplot.figure()
pyplot.plot(x,T-273.15,'b')

'''simulation HC fire curve'''
'''Initial Condition'''
T_A   = 303.15

'''Define Results Parameter'''
T    = np.ones(nx)*T_A

for n in range(1,nt):
        Tn = T.copy()
        T[1:-1] = Tn[1:-1] + dt * lamda/(rho*cp) * (Tn[2:] - 2*Tn[1:-1] +Tn[0:-2])/dx**2
        T[0] = 293.15 + 1060*(1 - 0.325*np.exp(-0.167*n*dt/60) - 0.675*np.exp(-2.5*n*dt/60))
        T[-1] = T[-2]
print (T-273.15)
pyplot.plot(x,T-273.15,'g')

'''simulation HCM fire curve'''
'''Initial Condition'''
T_A   = 303.15

'''Define Results Parameter'''
T    = np.ones(nx)*T_A

for n in range(1,nt):
        Tn = T.copy()
        T[1:-1] = Tn[1:-1] + dt * lamda/(rho*cp) * (Tn[2:] - 2*Tn[1:-1] +Tn[0:-2])/dx**2
        T[0] = 293.15 + 1280*(1 - 0.325*np.exp(-0.167*n*dt/60) - 0.675*np.exp(-2.5*n*dt/60))
        T[-1] = T[-2]
        if n%1440 == 0 :
                pyplot.plot(x,T-273.15,'r')
print (T-273.15)
pyplot.plot(x,T-273.15,'r')
pyplot.show(plotT)

''' Fire curve '''
T    = np.ones(nt)
for n in range(1,nt):
        T[n] = 293.15 + 345*np.log10(8*n*dt/60+1)
print (T-273.15)
pyplot.plot(T-273.15,'b')
T    = np.ones(nt)
for n in range(1,nt):
        T[n] = 293.15 + 1060*(1 - 0.325*np.exp(-0.167*n*dt/60) - 0.675*np.exp(-2.5*n*dt/60))
pyplot.plot(T-273.15,'g')
T    = np.ones(nt)
for n in range(1,nt):
        T[n] = 293.15 + 1280*(1 - 0.325*np.exp(-0.167*n*dt/60) - 0.675*np.exp(-2.5*n*dt/60))
print (T-273.15)
pyplot.plot(T-273.15,'r')
pyplot.show(plotT)
