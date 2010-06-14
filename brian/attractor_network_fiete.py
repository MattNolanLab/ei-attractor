from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from scipy import linspace


import time
import math

# define provisional model parameters - these might be changed in the future
C=200*pF
taum=10*msecond
gL=C/taum
EL=-70*mV
VT=-55*mV
DeltaT=3*mV

threshold = -20*mvolt;
refractory = 2*ms;

# Synapse parameters
Ee=0*mvolt
Ei=-80*mvolt
taue=5*msecond
taui=10*msecond

eqs=exp_IF(C,gL,EL,VT,DeltaT)
# Use only inhibitory connections from Burak&Fiete, 2009. Should work if the
# velocity input is non-zero even when speed is zero.
eqs+=exp_conductance('gi',Ei,taui) # from library.synapses
eqs+=Current('''B : amp''')

# Definition of 2d topography on a sheet and connections
sheet_size = 40;  # Total no. of neurons will be sheet_size^2

# Definition of connection parameters and setting the connection iteratively:
# Mexican hat connectivity
lambda_net = 13;
beta = 3.0 / lambda_net**2;
gamma = 10 * beta;
l = 2;
a = 1;

def translateIndexTo2d(index, size):
# index - linear index of the neuron
# size - size of the sheet (the 2d sheet is always square)
    return complex(index % size - size/2, index // size - size/2)

def getPreferredDirection(pos):
# pos - complex number defining position of neuron in 2d sheet
    pos4_x = pos.real % 4
    pos2_y = pos.imag % 2
    if pos4_x == 0:
        if pos2_y == 0:
            return complex(0, 1) # North
        else:
            return complex(1, 0) # East
    elif pos4_x == 1:
        if pos2_y == 0:
            return complex(0, -1) # South
        else:
            return complex(-1, 0) # West
    elif pos4_x == 2:
        if pos2_y == 0:
            return complex(1, 0)
        else:
            return complex(0, 1)
    else:
        if pos2_y == 0:
            return complex(-1, 0)
        else:
            return complex(0, -1)



sheetGroup=NeuronGroup(sheet_size**2,model=eqs,threshold=threshold,reset=EL,refractory=refractory)

inhibConn = Connection(sheetGroup, sheetGroup, 'gi', structure='dense');
for i in range(len(sheetGroup)):
    for j in range(len(sheetGroup)):
        pos_i = translateIndexTo2d(i, sheet_size);
        pos_j = translateIndexTo2d(j, sheet_size);
        prefDir = getPreferredDirection(pos_j); # all as complex numbers
        x = pos_i - pos_j
        abs_x = abs(x)
        if abs_x > sheet_size/2:
            abs_x = sheet_size - abs_x
        abs_x -= l*prefDir;
        #print math.e**(-gamma*abs(x)**2)
        w = a*math.e**(-gamma*(abs(x)**2)) - math.e**(-beta*(abs(x)**2));
        inhibConn[j, i] = abs(w)*nS


# Velocity inputs - for now zero velocity
input = 0.3*namp
sheetGroup.B = linspace(input, input, sheet_size**2)


# Record the number of spikes
M=SpikeMonitor(sheetGroup)
MV = StateMonitor(sheetGroup, 'vm', record = [480, 494])

sheetGroup.vm = EL + (VT-EL) * rand(len(sheetGroup))


print inhibConn.W[:][16]

print "Simulation running..."
start_time=time.time()
run(5000*ms)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"
print M.nspikes/sheet_size**2.,"spikes per neuron"
raster_plot(M)
figure();
plot(MV.times/ms, MV[494]/mV)
figure();
plot(MV.times/ms, MV[480]/mV)
show()
