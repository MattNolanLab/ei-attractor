from brian import *
from brian.library.IF import *
from brian.library.synapses import *
from brian.membrane_equations import *

from scipy import linspace
#from scipy.io import loadmat
from optparse import OptionParser
from datetime import datetime

#import time
#import math
#import random


def createStellateCells(options):
    st = options.stellate

    Ei = st.Ei
    taui = st.taui

    eqs = MembraneEquation(st.C) + leak_current(st.Gl, st.El)
    eqs += K_current_HH(st.GK, st.EK) + Na_current_HH(st.GNa, st.ENa)
    eqs += Current('gi*(Ei - vm): amp')
    eqs += 'dgi/dt = -gi/taui : siemens')
    eqs += Current('I:amp')

    return NeuronGroup(st.N, eqs, implicit=True, freeze=True)


def createInterneurons(options):
    int = options.interneurons

    Ee = int.Ee
    taue = int.taue

    eqs = MembraneEquation(int.C) + leak_current(int.Gl, int.El)
    eqs += K_current_HH(int.GK, int.EK) + Na_current_HH(int.GNa, int.ENa)
    eqs += Current('ge*(Ee - vm): amp')
    eqs += 'dge/dt = -ge/taue : siemens')
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
