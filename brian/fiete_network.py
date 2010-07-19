import brian_no_units
from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from scipy import linspace
from scipy.io import loadmat
from optparse import OptionParser
from datetime import datetime

import time
import math

# define provisional model parameters - these might be changed in the future

threshold = -20*mvolt;
refractory = 20*ms;

# Synapse parameters
Ee=0*mvolt
Ei=-80*mvolt
#taue=5*msecond



def get_exp_IF(C, gL, EL, VT, DeltaT, Ei, taui):
    eqs=exp_IF(C,gL,EL,VT,DeltaT)
    #eqs=leaky_IF(taum, EL)
    # Use only inhibitory connections from Burak&Fiete, 2009. Should work if the
    # velocity input is non-zero even when speed is zero.
    eqs+=exp_conductance('gi',Ei,taui) # from library.synapses
    eqs+=Current('''B : amp''')
    return eqs


# Get a preferred direction for a neuron
def getPreferredDirection(pos_x, pos_y):
# pos_x/y - position of neuron in 2d sheet
    pos4_x = pos_x % 4
    pos2_y = pos_y % 2
    if pos4_x == 0:
        if pos2_y == 0:
            return [0, 1] # North
        else:
            return [1, 0] # East
    elif pos4_x == 1:
        if pos2_y == 0:
            return [0, -1] # South
        else:
            return [-1, 0] # West
    elif pos4_x == 2:
        if pos2_y == 0:
            return [1, 0]
        else:
            return [0, 1]
    else:
        if pos2_y == 0:
            return [-1, 0]
        else:
            return [0, -1]

def createNetwork(sheet_size, lambda_net, l, a, connMult, clock, taum_ms, taui_ms):
    C=200*pF
    taum=taum_ms*msecond
    taui=taui_ms*msecond
    gL=C/taum
    EL=-70*mV
    VT=-55*mV
    DeltaT=3*mV

    beta = 3.0 / lambda_net**2
    gamma = 1.05 * beta

    sheetGroup=NeuronGroup(sheet_size**2,model=get_exp_IF(C, gL, EL, VT, DeltaT,
        Ei, taui),threshold=threshold,reset=EL,refractory=refractory,
        clock=clock)
    inhibConn = Connection(sheetGroup, sheetGroup, 'gi', structure='dense');
    inh_matrix = asarray(inhibConn.W)
    
    for j in xrange(len(sheetGroup)):
        j_x = j % sheet_size
        j_y = j // sheet_size
        prefDir = getPreferredDirection(j_x, j_y)
        for i in xrange(len(sheetGroup)):
            i_x = i % sheet_size
            i_y = i // sheet_size
            
            if abs(j_x - i_x) > sheet_size/2:
                if i_x > sheet_size/2:
                    i_x = i_x - sheet_size
                else:
                    i_x = i_x + sheet_size
            if abs(j_y - i_y) > sheet_size/2:
                if i_y > sheet_size/2:
                    i_y = i_y - sheet_size
                else:
                    i_y = i_y + sheet_size
    
            abs_x_sq = (i_x - j_x -l*prefDir[0])**2 + (i_y - j_y - l*prefDir[1])**2
            w = a*math.e**(-gamma*(abs_x_sq)) - math.e**(-beta*(abs_x_sq));
            inh_matrix[j, i] = connMult*abs(w)*nS

    sheetGroup.vm = EL + (VT-EL) * rand(len(sheetGroup))

    return [sheetGroup, inhibConn]


def printConn(sheet_size, conn, write_data, print_only_conn):
    # Plot connection matrix for neuron at position defined by row_i
    rows_i = [0.25*sheet_size**2 + 0.25*sheet_size, 0.25*sheet_size**2 +
            0.75*sheet_size, 0.75*sheet_size**2 + 0.25*sheet_size,
            0.75*sheet_size**2 + 0.75*sheet_size]
    plot_i = 1
    for row_i in rows_i:
        #col_i = row_i
        connRow = conn[row_i, :]
        #connCol = conn[:, col_i]
        
        x = arange(0, sheet_size)
        y = arange(0, sheet_size)
        X, Y = meshgrid(x, y)
        subplot(2,2, plot_i)
        contour(X, Y, reshape(connRow, (sheet_size, sheet_size)))
        #figure()
        #contour(X, Y, reshape(connCol, (sheet_size, sheet_size)))
        plot_i += 1
    xlabel("Neuron number")
    ylabel("Neuron number")
    if write_data:
        savefig(conn_fname, format="eps")
    elif print_only_conn == True:
        show()
        exit()

