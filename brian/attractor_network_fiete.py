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

from fiete_network import *


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



# Directory and filenames constants
timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
dirName = 'results/'
population_fname = dirName + timeSnapshot + '_spacePlot.eps'
options_fname = dirName + timeSnapshot + '.params'
count_fname = dirName + timeSnapshot + '_linePlot.eps'
conn_fname = dirName + timeSnapshot + "_conn.eps"



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
sheet_size = options.sheet_size  # Total no. of neurons will be sheet_size^2

# Definition of connection parameters and setting the connection iteratively:
# Mexican hat connectivity

start_time=time.time()

sheetGroup=NeuronGroup(sheet_size**2,model=get_exp_IF(C, gL, EL, VT, DeltaT, Ei,
    taui),threshold=threshold,reset=EL,refractory=refractory)
inhibConn = createConn(sheetGroup, options.sheet_size, options.lambda_net, options.l, options.a, options.connMult)


duration=time.time()-start_time
print "Connection setup time:",duration,"seconds"



# Velocity inputs - for now zero velocity
input = options.input*namp
sheetGroup.B = linspace(input, input, sheet_size**2)


# Record the number of spikes
Monitor = SpikeCounter(sheetGroup)
#MV = StateMonitor(sheetGroup, 'vm', record = [480, 494])

sheetGroup.vm = EL + (VT-EL) * rand(len(sheetGroup))

printConn(sheet_size, inhibConn, options.write_data, options.print_only_conn)

print "Simulation running..."
start_time=time.time()
net = Network(sheetGroup, inhibConn, Monitor)
net.run(5*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"

print Monitor.count

if (options.write_data):
    f = open(options_fname, 'w')
    f.write(str(options))
    f.close()


# print contour plot of firing
x = arange(0, sheet_size);
y = arange(0, sheet_size);
X, Y = meshgrid(x, y)
figure()
contour(X, Y, reshape(Monitor.count, (sheet_size, sheet_size)))
xlabel("Neuron number")
ylabel("Neuron number")
if options.write_data:
    savefig(population_fname, format="eps")

figure()
plot(Monitor.count)
xlabel("Neuron number")
ylabel("Spiking activity")
if options.write_data:
    savefig(count_fname, format="eps")

#figure();
#plot(MV.times/ms, MV[494]/mV)
#figure();
#plot(MV.times/ms, MV[480]/mV)

if options.write_data == False:
    show()
