import brian_no_units
from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from scipy import linspace
from optparse import OptionParser
from datetime import datetime


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

# Network parameters definition
optParser = OptionParser()
optParser.add_option("--lambda-net", type="float", default=13,
        dest="lambda_net")
optParser.add_option("-s", "--sheet-size", type="int", default=40, dest="sheet_size")
optParser.add_option("-l", type="float", default=2.0, dest="l")
optParser.add_option("-a", type="float", default=1.0, dest="a");
optParser.add_option("-c", "--conn-mult", type="float", default=10.0,
        dest="connMult")
optParser.add_option("-w", "--write-data", action="store_true",
        dest="write_data", default=False)
optParser.add_option("-i", "--input", type="float", default=0.3, dest="input",
        help="Feedforward constant input to all neurons")
optParser.add_option("-p", "--print-only-conn", action="store_true",
        dest="print_only_conn", default=False)

(options, args) = optParser.parse_args()
print "Options:"
print options


# Definition of 2d topography on a sheet and connections
sheet_size = options.sheet_size;  # Total no. of neurons will be sheet_size^2

# Definition of connection parameters and setting the connection iteratively:
# Mexican hat connectivity
lambda_net = options.lambda_net
beta = 3.0 / lambda_net**2
gamma = 1.05 * beta
l = options.l
a = options.a

connMult = options.connMult


#def translateIndexTo2d(index, size):
# index - linear index of the neuron
# size - size of the sheet (the 2d sheet is always square)
#    return complex(index % size - size/2, index // size - size/2)

def getPreferredDirection(pos_x, pos_y):
# pos - complex number defining position of neuron in 2d sheet
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



start_time=time.time()

sheetGroup=NeuronGroup(sheet_size**2,model=eqs,threshold=threshold,reset=EL,refractory=refractory)
inhibConn = Connection(sheetGroup, sheetGroup, 'gi', structure='dense');

for j in range(len(sheetGroup)):
    j_x = j % sheet_size
    j_y = j // sheet_size
    prefDir = getPreferredDirection(j_x, j_y) # all as complex numbers
    for i in range(len(sheetGroup)):
        i_x = i % sheet_size
        i_y = i // sheet_size
        
        #if abs(j_x - i_x) > sheet_size/2:
        #    if i_x > sheet_size/2:
        #        i_x = i_x - sheet_size
        #    else:
        #        i_x = i_x + sheet_size
        #if abs(j_y - i_y) > sheet_size/2:
        #    if i_y > sheet_size/2:
        #        i_y = i_y - sheet_size
        #    else:
        #        i_y = i_y + sheet_size

        abs_x = math.sqrt((i_x - j_x -l*prefDir[0])**2 + (i_y - j_y -
            l*prefDir[1])**2)
        if (abs_x > sheet_size/2):
            abs_x = sheet_size - abs_x

        w = a*math.e**(-gamma*(abs_x**2)) - math.e**(-beta*(abs_x**2));
        inhibConn[j, i] = connMult*abs(w)*nS

duration=time.time()-start_time
print "Connection setup time:",duration,"seconds"


# Plot connection matrix for neuron at position position defined by row_i
if options.print_only_conn == True:
    rows_i = [325, 355, 1405, 1435]

    for row_i in rows_i:
        #col_i = row_i
        connRow = inhibConn[row_i, :]
        #connCol = inhibConn[:, col_i]
        
        x = arange(0, sheet_size)
        y = arange(0, sheet_size)
        X, Y = meshgrid(x, y)
        figure()
        contour(X, Y, reshape(connRow, (sheet_size, sheet_size)))
        #figure()
        #contour(X, Y, reshape(connCol, (sheet_size, sheet_size)))
    show()
    exit()

# Velocity inputs - for now zero velocity
input = options.input*namp
sheetGroup.B = linspace(input, input, sheet_size**2)


# Record the number of spikes
Monitor = SpikeCounter(sheetGroup)
#MV = StateMonitor(sheetGroup, 'vm', record = [480, 494])

sheetGroup.vm = EL + (VT-EL) * rand(len(sheetGroup))


print "Simulation running..."
start_time=time.time()
run(5*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"

print Monitor.count

timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
dirName = 'results/'
population_fname = dirName + timeSnapshot + '_spacePlot.eps'
options_fname = dirName + timeSnapshot + '.params'
count_fname = dirName + timeSnapshot + '_linePlot.eps'

if (options.write_data):
    f = open(options_fname, 'w')
    f.write(str(options))
    f.close()

# print contour plot of firing
x = arange(0, sheet_size);
y = arange(0, sheet_size);
X, Y = meshgrid(x, y)
contour(X, Y, reshape(Monitor.count, (sheet_size, sheet_size)))
if options.write_data:
    savefig(population_fname, format="eps")

figure()
plot(Monitor.count)
if options.write_data:
    savefig(count_fname, format="eps")

#figure();
#plot(MV.times/ms, MV[494]/mV)
#figure();
#plot(MV.times/ms, MV[480]/mV)

if options.write_data == False:
    show()
