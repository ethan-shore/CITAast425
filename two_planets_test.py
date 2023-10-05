#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 15:51:01 2023

@author: ethanshore

Python script to test two planet system simulations for the cluster
"""

import rebound
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import gridspec

m1 = 5e-5 # Neptune mass
m2 = 5e-5
delta = 3*(m1 + m2) ** (1/3) # equation 10

def setupSimulation(mass1, a1, f1, r1, mass2, a2, f2, r2):
    sim = rebound.Simulation()
    sim.units = ("Msun","AU","yr")
    sim.add(m=1.)                      # Solar Mass Star
    sim.add(m=mass1,a=a1, f=f1, r=r1)  # Neptune 1
    sim.add(m=mass2,a=a2, f=f2, r=r2)  # Neptune 2
    sim.move_to_com()
    return sim

def runSimulation(simulation, min_dist, max_dist, iterations):
    sim.exit_min_distance = min_dist
    sim.exit_max_distance = max_dist
    Noutputs = 1000
    times = np.linspace(0,iterations*2.*np.pi,Noutputs) 
    distances = np.zeros(Noutputs)
    a1 = np.zeros(Noutputs) # initializing arrays for a    
    a2 = np.zeros(Noutputs)
    ps = simulation.particles   # ps is now an array of pointers. It will update as the simulation runs.
    try:
        for i,time in enumerate(times):
            simulation.integrate(time)
            dp = ps[1] - ps[2]   # Calculates the coponentwise difference between particles 
            distances[i] = np.sqrt(dp.x*dp.x+dp.y*dp.y+dp.z*dp.z)
            a1[i] = ps[1].a
            a2[i] = ps[2].a
    except rebound.Escape as error:
        print(error, distances[0])
        # plot orbit with escaped particle
        '''
        %matplotlib inline
        op = rebound.OrbitPlot(sim,Narc=300)

        # plot distance between planets
        plt.figure()
        plt.plot(times/(2.*np.pi), distances, lw=0.5)
        plt.xlabel('Time [Orbits]')
        plt.ylabel('Distance [AU]')
        '''
        
    except rebound.Encounter as error:
        print(error)
        
deltaFactor = float(sys.argv[1])


sim = setupSimulation(m1, 30.05, 0, 1.65e-4, m2, 30.05*(1. + deltaFactor * delta), np.pi, 1.65e-4)
sim.ri_ias15.adaptive_mode = 2
Tfin = float(sys.argv[2])
print("Argument List:", str(sys.argv))
Nout  = 512*4
interval = Tfin / Nout
sasave = '/Users/ethanshore/Downloads/AST425_SH/run_' + sys.argv[1] + '_'+ sys.argv[2]+'.sa'
sim.automateSimulationArchive(sasave,interval=interval,deletefile=True)
sim.integrate(Tfin)

sa = rebound.SimulationArchive(sasave)
N = len(sa)

times = np.zeros(N)
a1,a2 = np.zeros((2,N))
e1,e2 = np.zeros((2,N))
E = np.zeros(N)
for i,sim in enumerate(sa):
    ps = sim.particles
    times[i] = sim.t
    a1[i] = ps[1].a
    a2[i] = ps[2].a
    e1[i] = ps[1].e
    e2[i] = ps[2].e
    E[i] = sim.energy()

fig = plt.figure()

# set up the figure for two plots
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 

# add the semi major axes and apogees/perigees
ax0 = plt.subplot(gs[0])
line00, = ax0.plot(times, a1, 'k')
line01, = ax0.plot(times,a1*(1-e1),'k--')
line02, = ax0.plot(times,a1*(1+e1),'k--')

line03, = ax0.plot(times, a2, 'r')
line04, = ax0.plot(times,a2*(1-e2),'r--')
line05, = ax0.plot(times,a2*(1+e2),'r--')

ax0.set_ylabel('Semimajor axis [AU]')
plt.title(r'Separation $=${}$\Delta $'.format(deltaFactor))

# add the energy
ax1 = plt.subplot(gs[1], sharex = ax0)
line10, = ax1.plot(times, E/np.mean(E)-1, 'blue')

ax1.set_ylabel('Normalized Energy')
ax1.set_xlabel('Time [yr]')

plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
ax0.legend((line00, line01, line03,
            line04, line10), ('Planet 1', 'Apogee/perigee 1',
                              'Planet 2', 'Apogee/perigee 2', 'Energy'),
                              prop = { "size": 6 })

plt.subplots_adjust(hspace=.0)

filename = '/Users/ethanshore/Downloads/AST425_SH/testdelta_' + sys.argv[1] + '_'+ sys.argv[2]+'.jpg'

plt.savefig(filename)
#plt.show()
'''
for i in $(seq 0.4 0.1 0.4); do python /Users/ethanshore/Downloads/AST425_SH/two_planets_test.py $i 1e7; done
'''