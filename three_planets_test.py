#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:27:04 2023

@author: ethanshore

Python script to test three planet system simulations
"""

import rebound
import numpy as np
#import matplotlib.pyplot as plt
import sys
#from matplotlib import gridspec
savedir = './'

m1 = 5e-5 # Neptune mass
m2 = 5e-5
m3 = 5e-5

# from p. 13 of Obertas et al., checking for Mutual hill radii between 6 and 8
delta = 5.9 + float(sys.argv[1]) * 0.1 

def setupSimulation(mass1, a1, f1, r1, mass2, a2, f2, r2, mass3, a3, f3, r3):
    sim = rebound.Simulation()
    sim.units = ("Msun","AU","yr")
    sim.add(m=1.)                      # Solar Mass Star
    sim.add(m=mass1,a=a1, f=f1, r=r1)  # Neptune 1
    sim.add(m=mass2,a=a2, f=f2, r=r2)  # Neptune 2
    sim.add(m=mass3, a=a3, f=f3, r=r3)
    sim.move_to_com()
    return sim

def runSimulation(simulation, min_dist, max_dist, iterations):
    sim.exit_min_distance = min_dist
    sim.exit_max_distance = max_dist
    Noutputs = 1000
    times = np.linspace(0,iterations*2.*np.pi,Noutputs) 
    distances12 = np.zeros(Noutputs)
    distances13 = np.zeros(Noutputs)
    distances23 = np.zeros(Noutputs)
    a1 = np.zeros(Noutputs) # initializing arrays for a    
    a2 = np.zeros(Noutputs)
    a3 = np.zeros(Noutputs)
    ps = simulation.particles   # ps is now an array of pointers. It will update as the simulation runs.
    try:
        for i,time in enumerate(times):
            simulation.integrate(time)
            dp12 = ps[1] - ps[2]   # Calculates the coponentwise difference between particles 
            dp13 = ps[1] - ps[3]
            dp23 = ps[2] - ps[3]
            distances12[i] = np.sqrt(dp12.x*dp12.x+dp12.y*dp12.y+dp12.z*dp12.z)
            distances13[i] = np.sqrt(dp13.x*dp13.x+dp13.y*dp13.y+dp13.z*dp13.z)
            distances23[i] = np.sqrt(dp23.x*dp23.x+dp23.y*dp23.y+dp23.z*dp23.z)
            a1[i] = ps[1].a
            a2[i] = ps[2].a
            a3[i] = ps[3].a
    except rebound.Escape as error:
        print(error, distances12[0], distances13[0], distances23[0])
        
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
        
#deltaFactor = float(sys.argv[1])


sim = setupSimulation(m1, 30.05, 0, 1.65e-4, m2, 30.05*(1. + delta), (2/3) * np.pi, 1.65e-4, m3, 30.05*(1. + 2 * delta), (4/3)*np.pi, 1.65e-4)
sim.ri_ias15.adaptive_mode = 2
Tfin = float(sys.argv[2])
print("Argument List:", str(sys.argv))
Nout  = 512*4
interval = Tfin / Nout
sasave = savedir + '/run_' + delta + '_'+ sys.argv[2]+'.sa'
sim.automateSimulationArchive(sasave,interval=interval,deletefile=True)
sim.integrate(Tfin)

'''

sa = rebound.SimulationArchive(sasave)
N = len(sa)

times = np.zeros(N)
a1,a2,a3 = np.zeros((3,N))
e1,e2,e3 = np.zeros((3,N))
E = np.zeros(N)
for i,sim in enumerate(sa):
    ps = sim.particles
    times[i] = sim.t
    a1[i] = ps[1].a
    a2[i] = ps[2].a
    a3[i] = ps[3].a
    e1[i] = ps[1].e
    e2[i] = ps[2].e
    e3[i] = ps[3].e
    E[i] = sim.energy()

fig = plt.figure()

# set up the figure for two plots
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 

# add the semi major axes and apogees/perigees
ax0 = plt.subplot(gs[0])
# planet 1
line00, = ax0.plot(times, a1, 'k')
line01, = ax0.plot(times,a1*(1-e1),'k--', alpha=0.5)
line02, = ax0.plot(times,a1*(1+e1),'k--', alpha=0.5)
# planet 2
line03, = ax0.plot(times, a2, 'r')
line04, = ax0.plot(times,a2*(1-e2),'r--', alpha=0.5)
line05, = ax0.plot(times,a2*(1+e2),'r--', alpha=0.5)
# planet 3
line06, = ax0.plot(times, a3, 'g')
line07, = ax0.plot(times,a3*(1-e3),'g--', alpha=0.5)
line08, = ax0.plot(times,a3*(1+e3),'g--', alpha=0.5)

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
            line04, line06, line07, line10), ('Planet 1', 'Apogee/perigee 1',
                              'Planet 2', 'Apogee/perigee 2', 'Planet 3',
                              'Apogee/perigee 3', 'Energy'),
                              loc='upper left', prop = { "size": 6 })

plt.subplots_adjust(hspace=.0)

filename = '/Users/ethanshore/Downloads/AST425_SH/Three_Planet_sims/threeplanet_testdelta_' + delta + '_'+ sys.argv[2]+'.jpg'

plt.savefig(filename)
#plt.show()
'''
'''
for i in $(seq 0.7 0.1 0.7); do python /Users/ethanshore/Downloads/AST425_SH/Three_Planet_sims/three_planets_test.py $i 1e7; done
'''