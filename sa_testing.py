#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:22:48 2023

@author: ethanshore
"""

import rebound
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec

sa = rebound.SimulationArchive('/Users/ethanshore/Documents/AST425_SH/run_0.5_1e6.sa')
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
    E[i] = sim.calculate_energy()

fig = plt.figure()

gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 

ax0 = plt.subplot(gs[0])
line00, = ax0.plot(times, a1, 'k')
line01, = ax0.plot(times,a1*(1-e1),'k--', lw=0.5)
line02, = ax0.plot(times,a1*(1+e1),'k--', lw=0.5)

line03, = ax0.plot(times, a2, 'r')
line04, = ax0.plot(times,a2*(1-e2),'r--', lw=0.5)
line05, = ax0.plot(times,a2*(1+e2),'r--', lw=0.5)

ax0.set_ylabel('Semimajor axis [AU]')
plt.title(r'Separation $=${}$\Delta $'.format(1.0))


ax1 = plt.subplot(gs[1], sharex = ax0)
line10, = ax1.plot(times, E/np.mean(E)-1, 'blue')

ax1.set_ylabel('Normalized Energy')
ax1.set_xlabel('Time [yr]')

plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

plt.subplots_adjust(hspace=.0)
#filename = '/Users/ethanshore/Documents/AST425_SH/testdelta' + sys.argv[1] + '.jpg'
#plt.savefig(filename)
#plt.show()