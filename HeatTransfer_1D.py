# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a Heat Transfer script file.
"""

import numpy
from matplotlib import pyplot

''' Geometry '''
s_B  = 0.2

''' Spacial grid '''
dx   =  0.01
nx   = int(s_B/dx) + 1
x    = numpy.linspace(0,s_B,nx)

'''Time grid'''
t_sim = 7200
dt   =  0.5
nt   = int(t_sim/dt)

'''Material Property'''
rho   = 2400
cp    = 1000
lamda = 1.13

'''Initial Condition'''
T_A   = 303.15

'''Define Results Parameter'''
T    = numpy.ones(nx)*T_A

'''simulation'''
for n in range(1,nt):
        Tn = T.copy()
        T[1:-1] = Tn[1:-1] + dt * lamda/(rho*cp) * (Tn[2:] - 2*Tn[1:-1] +Tn[0:-2])/dx**2
        T[0] = 1573.
        T[-1] = T[-2]

#plotT = pyplot.figure(figsize=(6,7),dpi=100)
plotT = pyplot.figure()
pyplot.plot(x,T-273.15,'k')
pyplot.show(plotT)
