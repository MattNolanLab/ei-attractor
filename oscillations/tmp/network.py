from brian import *
from brian.library.IF import *
from brian.library.synapses import *
from brian.membrane_equations import *
from brian.library.ionic_currents import *

from scipy import linspace
#from scipy.io import loadmat
from optparse import OptionParser
from datetime import datetime

#import time
#import math
#import random

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

options = Bunch(
        stellate = Bunch(
            Ei   = -10 * mV,
            taui = 10 * msecond,
            C    = 1 * uF,
            gl   = 0.3 * msiemens,
            El   = 10.6 * mV,
            gK   = 36 * msiemens,
            EK   = -12 * mV,
            gNa  = 120 * msiemens,
            ENa  = 120 * mV,
            N    = 100),
        interneurons = Bunch(
            Ee   = 10 * mV,
            taue = 2 * msecond,
            C    = 1 * uF,
            gl   = 0.3 * msiemens,
            El   = 10.6 * mV,
            gK   = 36 * msiemens,
            EK   = -12 * mV,
            gNa  = 120 * msiemens,
            ENa  = 120 * mV,
            N    = 100),
        exc_weight = 0,
        inh_weight = 0)


options.exc_weight = 1/(options.stellate.N + options.stellate.N)
options.inh_weight = 1/(options.stellate.N + options.stellate.N)

def createStellateCells(options):
    st = options.stellate

    eqs = MembraneEquation(st.C) + leak_current(st.gl, st.El)
    eqs += K_current_HH(st.gK, st.EK) + Na_current_HH(st.gNa, st.ENa)
    eqs += Current('gi*(Ei - vm) : amp', Ei=st.Ei)
    eqs += Equations('dgi/dt = -gi/taui : siemens', taui=st.taui)
    eqs += Current('I:amp')

    return NeuronGroup(st.N, eqs, implicit=True, freeze=True)


def createInterneurons(options):
    int = options.interneurons

    Ee = int.Ee
    taue = int.taue

    eqs = MembraneEquation(int.C) + leak_current(int.gl, int.El)
    eqs += K_current_HH(int.gK, int.EK) + Na_current_HH(int.gNa, int.ENa)
    eqs += Current('ge*(Ee - vm): amp')
    eqs += 'dge/dt = -ge/taue : siemens'
    eqs += Current('I:amp')

    return NeuronGroup(int.N, eqs, implicit=True, freeze=True)


def createNetwork(options):

    stellates = createStellateCells(options)
    interneurons = createInterneurons(options)

    # Connect the populations together, stellates exciting interneurons,
    # interneuons exciting stellates
    C_exc = Connection(stellates, interneurons,
            state='ge',
            delay=options.delay,
            spareness=options.sparseness,
            weight=options.exc_weight)

    C_inh = Connection(interneurons, stellates,
            state='gi',
            delay=options.delay,
            spareness=options.sparseness,
            weight=options.inh_weight)



###############################################################################

net = createNetwork(options)


