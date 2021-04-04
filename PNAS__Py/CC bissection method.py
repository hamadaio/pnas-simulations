#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 10:41:26 2016

@author: sarahgoethals
"""

from brian2 import *

matplotlib.rcParams['figure.figsize'] = (10, 8)

defaultclock.dt = 0.015*ms

# Passive parameters
EL = -80.*mV
Cm = 0.9*uF/cm**2
gL = 1.*(siemens/meter**2)
Ri = 150.*ohm*cm

# Na channels parameters
ENa = 60.*mV
ka = 6.*mV
va = -30.*mV
gNa_surf = 3000.*(siemens/meter**2)
taum = 0.05*ms
Rm = Cm*taum

# Morphology: somatodendritic cylinder and axon
# Rmq: the "soma" is located at the 499th compartment.
dend_length = 500.*um
axon_diam = 1.5*um
axon_length = 300.*um

# Na channels distribution
AIS_length = 45.*um


# Bissection method
dv_max_required = 350.
Na_starts = linspace(5.,40.,36)*um 
d_dend_values = []

for Na_start in Na_starts:
    d_dend_min = 2.*um
    d_dend_max = 8.*um
    d_dend_current = 5.*um
    while True:
        
        # Create SpatialNeuron with x_ais and d_dend_current
        # simulate it, recording V
        
        morpho = Cylinder(diameter=d_dend_current, length=dend_length, n=500)
        axon = Cylinder(diameter=axon_diam, length=axon_length, n=300)
        morpho.axon = axon
    
        AIS_position = dend_length + Na_start
        Na_end = AIS_position + AIS_length
    
        # Equations
        eqs='''
        Im = gL*(EL - v) + gNa*m*(ENa - v) : amp/meter**2
        dm/dt = (minf - m) / taum: 1  # simplified Na channel
        minf = 1 / (1 + exp((va - v) / ka)) : 1
        gNa : siemens/meter**2
        I : amp (point current)
        '''
        neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                       method="exponential_euler")

        # Initial segment
        initial_segment = morpho.axon[Na_start:Na_end-dend_length]
        neuron.gNa[initial_segment] = gNa_surf
    
        # Initialisation
        neuron.I = 0*amp
        neuron.v= EL

        # Recording the voltage across compartments and time 
        M=StateMonitor(neuron,['v','m'],record=True)

        # Current injected at the junction between the somatodendritic compartment and the axon
        run(10*ms)
        neuron.I[499] = 0.2*nA
        run(9*ms)   
        
        # Max of dv/dt at the soma
        dv_max = max(diff(M[499].v)/defaultclock.dt)/(mV/ms)
        
        #print Na_start, d_dend_current, dv_max
        
        # Bisection method
        if abs(dv_max - dv_max_required) < 0.05: 
            break
        if dv_max < dv_max_required:
            d_dend_max = d_dend_current
        else:
            d_dend_min = d_dend_current
        
        d_dend_current = (d_dend_max - d_dend_min)/2 + d_dend_min
    
    d_dend_values.append(d_dend_current)