from brian import *
from brian.library.IF import *
from brian.library.synapses import *
from brian.membrane_equations import *

from scipy import linspace
from scipy.io import loadmat
from optparse import OptionParser
from datetime import datetime

import time
import math
import random


# Get a preferred direction for a neuron
def getPreferredDirection1D(pos_x):
# pos_x/y - position of neuron in 2d sheet
    pos2_x = pos_x % 2
    if pos2_x == 0:
        return 1  # right
    else:
        return -1 # left

#def getPreferredDirectionRandom(pos_x, pos_y):
#    # return random preferred direction on the sheet
#    return random.choice([[0, 1], [0, -1], [-1, 0], [1, 0]]);


def createNetwork1D(options, clock, W = None):

    C=200*pF
    refractory = 20*ms;
    
    # Synapse parameters
    Ee=0*mvolt
    Ei=-80*mvolt

    sheet_size = options.sheet_size
    l = options.l
    a = options.a
    connMult = options.connMult

    taum = options.taum*msecond
    taui = options.taui*msecond
    threshold = options.threshold*mvolt
    noise_sigma = options.noise_sigma*mvolt

    gL=C/taum
    EL=options.EL*mV
    VT=-55*mV
    DeltaT=3*mV

    beta = 3.0 / options.lambda_net**2
    gamma = 1.05 * beta

    # Setup neuron equations
    # Using exponential integrate and fire model
    eqs = '''
        dvm/dt = 1/C*Im + (noise_sigma*xi/taum**.5): volt
        Im = gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT) + gi*(Ei - vm) + B  : amp
        dgi/dt = -gi/taui : siemens
        B : amp
        '''

    sheetGroup = NeuronGroup(
            options.sheet_size,
            model=eqs,
            threshold=threshold,
            reset=EL,
            refractory=refractory,
            clock=clock)

    if (W == None):
        inhibConn = Connection(sheetGroup, sheetGroup, 'gi', structure='dense');
        inh_matrix = asarray(inhibConn.W)
    
        # Create toroidal connections matrix on the 2d sheet
        for j in xrange(len(sheetGroup)):
            j_x = j;
            prefDir = getPreferredDirection1D(j_x)
            #prefDir = getPreferredDirectionRandom(j_x, j_y)
            for i in xrange(len(sheetGroup)):
                i_x = i 
                if abs(j_x - i_x) > sheet_size/2:
                    if i_x > sheet_size/2:
                        i_x = i_x - sheet_size
                    else:
                        i_x = i_x + sheet_size
        
                abs_x_sq = (i_x - j_x -l*prefDir)**2
                w = a*math.e**(-gamma*(abs_x_sq)) - math.e**(-beta*(abs_x_sq));
                inh_matrix[j, i] = connMult*abs(w)*1e-9
    else:
        inhibConn = Connection(sheetGroup, sheetGroup, 'gi', structure='dense')
        print 'Initializing connections from file...'
        inhibConn.connect(sheetGroup, sheetGroup, W)


    # Initialize membrane potential randomly
    if (noise_sigma == 0):
        sheetGroup.vm = EL + (VT-EL) * rand(len(sheetGroup))
    else:
        sheetGroup.vm = EL + zeros(len(sheetGroup))

    return [sheetGroup, inhibConn]


#def printConn(sheet_size, conn, write_data, print_only_conn):
#    # Plot connection matrix for neuron at position defined by row_i
#    rows_i = [0.25*sheet_size**2 + 0.25*sheet_size, 0.25*sheet_size**2 +
#            0.75*sheet_size, 0.75*sheet_size**2 + 0.25*sheet_size,
#            0.75*sheet_size**2 + 0.75*sheet_size]
#    plot_i = 1
#    for row_i in rows_i:
#        #col_i = row_i
#        connRow = conn[row_i, :]
#        #connCol = conn[:, col_i]
#        
#        x = arange(0, sheet_size)
#        y = arange(0, sheet_size)
#        X, Y = meshgrid(x, y)
#        subplot(2,2, plot_i)
#        contour(X, Y, reshape(connRow, (sheet_size, sheet_size)))
#        #figure()
#        #contour(X, Y, reshape(connCol, (sheet_size, sheet_size)))
#        plot_i += 1
#    xlabel("Neuron number")
#    ylabel("Neuron number")
#    if write_data:
#        savefig(conn_fname, format="eps")
#    elif print_only_conn == True:
#        show()
#        exit()
#
