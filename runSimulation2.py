# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 11:57:16 2015

@author: leonardo
"""

import matplotlib.pyplot as plt
import numpy
import pyNN.neuron
import time
from neuron import *
from pyNN.neuron.cells import NativeCellType
from pyNN.random import RandomDistribution
from nrnutils import Mechanism, Section

#load_mechanisms("/home/leonardo/workspace/network/")

class SpinalNeuron(object):
    
    def __init__(self, **parameters):
        # define mechanisms, parameters and topology
        
        # for soma
        traub = Mechanism('traub', gl=parameters['gl_soma'], el=0.0,
                    gnabar=parameters['gnabar'], gkfbar=parameters['gkfbar'],
                    gksbar=parameters['gksbar'], vtraub=parameters['vtraub'])
        
        self.soma = Section(L=parameters['Ls'], diam=parameters['Ds'],
                    Ra=parameters['Ra_s'], mechanisms=[traub])                
        
        self.soma.ena = parameters['ena']
        self.soma.ek = parameters['ek']
                    
        # for dendrite
        self.is_motoneuron = parameters['is_motoneuron']
        self.is_active = parameters['is_active']
        
        if parameters['is_active'] == 1:
            active = Mechanism('active', gl=parameters['gl_dend'], el=0.0,
                            gcabar=parameters['gcabar'], theta=parameters['theta'])
            self.dend = Section(L=parameters['Ld'], diam=parameters['Dd'],
                            Ra=parameters['Ra_d'], mechanisms=[active])
            self.dend.eca = parameters['eca']
        else:
            pas = Mechanism('pas', g=parameters['gl_dend'], e=0.0)
            self.dend = Section(L=parameters['Ld'], diam=parameters['Dd'],
                            Ra=parameters['Ra_d'], mechanisms=[pas])

        if parameters['is_motoneuron'] == 1:
            self.dend.connect(self.soma, 0, 1)
            self.dend.add_synapse('ampa', 'Exp2Syn', e=0.0, tau1=0.1, tau2=5.0)
        
        #needed for PyNN
        self.source_section = self.soma
        self.source = self.soma(0.5)._ref_v
        self.parameter_names = ('gl_soma', 'gnabar', 'gkfbar', 'gksbar', 'ena',
                                'ek', 'vtraub', 'is_motoneuron', 'Ls', 'Ds',
                                'Ra_s', 'is_active', 'gl_dend', 'gcabar', 'eca',
                                'theta', 'Ld', 'Dd', 'Ra_d')
        self.traces = {}
        self.recording_time = False
    
    def _set_gl_soma(self, value):
        for seg in self.soma:
            seg.traub.gl = value
    def _get_gl_soma(self):
            return self.soma(0.5).traub.gl
    gl_soma = property(fget=_get_gl_soma, fset=_set_gl_soma)
        
    def _set_gnabar(self, value):
        for seg in self.soma:
            seg.traub.gnabar = value
    def _get_gnabar(self):
            return self.soma(0.5).traub.gnabar
    gnabar = property(fget=_get_gnabar, fset=_set_gnabar)
    
    def _set_gkfbar(self, value):
        for seg in self.soma:
            seg.traub.gkfbar = value
    def _get_gkfbar(self):
            return self.soma(0.5).traub.gkfbar
    gkfbar = property(fget=_get_gkfbar, fset=_set_gkfbar)

    def _set_gksbar(self, value):
        for seg in self.soma:
            seg.traub.gksbar = value
    def _get_gksbar(self):
            return self.soma(0.5).traub.gksbar
    gksbar = property(fget=_get_gksbar, fset=_set_gksbar)

    def _set_ena(self, value):
        self.soma.ena = value
    def _get_ena(self):
        return self.soma.ena
    ena = property(fget=_get_ena, fset=_set_ena)

    def _set_ek(self, value):
        self.soma.ek = value
    def _get_ek(self):
        return self.soma.ek
    ek = property(fget=_get_ek, fset=_set_ek)

    def _set_vtraub(self, value):
        for seg in self.soma:
            seg.traub.vtraub = value
    def _get_vtraub(self):
            return self.soma(0.5).traub.vtraub
    vtraub = property(fget=_get_vtraub, fset=_set_vtraub)

    def _set_Ls(self, value):
        self.soma.L = value
    def _get_Ls(self):
        return self.soma.L
    Ls = property(fget=_get_Ls, fset=_set_Ls)

    def _set_Ds(self, value):
        self.soma.diam = value
    def _get_Ds(self):
        return self.soma.diam
    Ds = property(fget=_get_Ds, fset=_set_Ds)

    def _set_Ra_s(self, value):
        self.soma.Ra = value
    def _get_Ra_s(self):
        return self.soma.Ra
    Ra_s = property(fget=_get_Ra_s, fset=_set_Ra_s)

    def _set_gl_dend(self, value):
        for seg in self.dend:
            if parameters['is_active'] == 1:
                seg.active.gl = value
            else:
                seg.pas.g= value
    def _get_gl_dend(self):
        if parameters['is_active'] == 1:
            return self.dend(0.5).active.gl
        else:
            return self.dend(0.5).pas.g                
    gl_dend = property(fget=_get_gl_dend, fset=_set_gl_dend)

    def _set_gcabar(self, value):
        for seg in self.dend:
            if parameters['is_active'] == 1:
                seg.active.gcabar = value
            else:
                return 0
    def _get_gcabar(self):
        if parameters['is_active'] == 1:
            return self.dend(0.5).active.gcabar
        else:
            return 0
    gcabar = property(fget=_get_gcabar, fset=_set_gcabar)

    def _set_eca(self, value):
        self.dend.eca = value
    def _get_eca(self):
        return self.dend.eca
    eca = property(fget=_get_eca, fset=_set_eca)
        
    def _set_theta(self, value):
        for seg in self.dend:
            if parameters['is_active'] == 1:
                seg.active.theta = value
            else:
                return 0
    def _get_theta(self):
        if parameters['is_active'] == 1:
            return self.dend(0.5).active.theta
        else:
            return 0
    theta = property(fget=_get_theta, fset=_set_theta)

    def _set_Ld(self, value):
        self.dend.L = value
    def _get_Ld(self):
        return self.dend.L
    Ld = property(fget=_get_Ld, fset=_set_Ld)

    def _set_Dd(self, value):
        self.dend.diam = value
    def _get_Dd(self):
        return self.dend.diam
    Dd = property(fget=_get_Dd, fset=_set_Dd)

    def _set_Ra_d(self, value):
        self.dend.Ra = value
    def _get_Ra_d(self):
        return self.dend.Ra
    Ra_d = property(fget=_get_Ra_d, fset=_set_Ra_d)        
    
    def memb_init(self):
        for sec in (self.soma, self.dend):
            for seg in sec:
                seg.v = self.v_init

class MotoneuronType(NativeCellType):
    default_parameters = {'gl_soma': 1/1150.0, 'gnabar': 0.030, 'gkfbar': 0.004,
                          'gksbar': 0.016, 'ena': 120.0, 'ek': -10.0,
                          'vtraub': 0.0, 'is_motoneuron': 1, 'Ls': 77.5,
                          'Ds': 77.5, 'Ra_s': 70.0, 'is_active': 0,
                          'gl_dend': 1/14400.0, 'gcabar': 0.038, 'eca': 140.0,
                          'theta': 25.0, 'Ld': 5500.0, 'Dd': 41.5, 'Ra_d': 70.0}
    default_initial_values = {'v': 0.0}
    recordable = ['spikes', 'soma(0.5).v', 'dend(0.5).v']
    model = SpinalNeuron

nrn = pyNN.neuron

dt = 0.025
nrn.setup()
nrn.setup(timestep=dt)

# creating a population of neurons

num_MN = 120

# constant parameters
parameters = {'gnabar': 0.030, 'gkfbar': 0.004, 'ena': 120.0, 'ek': -10.0,
              'vtraub': 0.0, 'is_motoneuron': 1, 'Ra_s': 70.0, 'is_active': 0,
              'eca': 140.0, 'theta': -40.0, 'Ra_d': 70.0}
p1 = nrn.Population(num_MN, MotoneuronType, parameters)

# range of parameters
gl_soma_range = numpy.linspace(1/1150.0, 1/1050.0, num=num_MN)
gksbar_range = numpy.linspace(0.016, 0.025, num=num_MN)
Ls_range = numpy.linspace(77.5, 82.5, num=num_MN)
Ds_range = numpy.linspace(77.5, 82.5, num=num_MN)
gl_dend_range = numpy.linspace(1/14400.0, 1/10700.0, num=num_MN)
gcabar_range = numpy.linspace(0.00005, 0.00005, num=num_MN)
Ld_range = numpy.linspace(5500.0, 6800.0, num=num_MN)
Dd_range = numpy.linspace(41.5, 62.5, num=num_MN)

p1.tset('gl_soma', gl_soma_range)
p1.tset('gksbar', gksbar_range)
p1.tset('Ls', Ls_range)
p1.tset('Ds', Ds_range)
p1.tset('gl_dend', gl_dend_range)
p1.tset('gcabar', gcabar_range)
p1.tset('Ld', Ld_range)
p1.tset('Dd', Dd_range)

p1.initialize('v', 0.0)

pulse = nrn.DCSource(amplitude=3.0, start=100.0, stop=5000.0)
pulse.inject_into(p1)
p1.inject(pulse)

#times = numpy.arange(0.0,5000.0,0.025)
#amplitudes = 3*((1.0+0.5*numpy.sin(2.0*numpy.pi*1.0*times/1000.0))*numpy.sin(2.0*numpy.pi*20.0*times/1000.0)) + 2.0
#amplitudes = 1.5*numpy.sin(2.0*numpy.pi*2.0*times/1000.0) + 3.0
#plt.plot(times, amplitudes)
#plt.show()
#sine_wave = nrn.StepCurrentSource(times,amplitudes)
#p1.inject(sine_wave)
#sine_wave.inject_into(p1)

p1._record('soma(0.5).v')

#prj = nrn.Projection(p1, p1, nrn.AllToAllConnector(), target='soma.ampa')
#t0 = time.time()
nrn.run(5000.0)

#t1 = time.time()
#print (t1-t0)/num_MN
id, t, v = p1.recorders['soma(0.5).v'].get().T

#plot the results
plt.plot(t[0:200000.0-1],v[0:200000.0-1])
plt.show()
plt.plot(t[1800000.0:2000000.0-1],v[1800000.0:2000000.0-1])
plt.show()

import scipy.io as io
io.savemat('Membrane1',{'v': v})