#
#   simulation_bump_fitting.py
#
#   Main simulation run: Fitting a Gaussian to the bump and frequency analysis.
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from numpy.random       import choice

from models.parameters  import getOptParser
from models.gc_net_nest import BasicGridCellNetwork


parser          = getOptParser()
parser.add_option("--gammaNSample",   type="float",   help="Fraction of neurons in the network to sample from, for the frequency analysis.")

(options, args) = parser.parse_args()


################################################################################

ei_net = BasicGridCellNetwork(options, simulationOpts=None)

const_v = [00.0, 0.0]
ei_net.setConstantVelocityCurrent_e(const_v)


NSample = int(options.gammaNSample * ei_net.net_Ne)
stateRecF_e = choice(ei_net.E_pop, NSample, replace=False)

stateMonF_params = {
        'withtime' : True,
        'interval' : options.sim_dt*10,
        'record_from' : ['I_clamp_GABA_A']
}
stateMonF_e = ei_net.getGenericStateMonitor(stateRecF_e, stateMonF_params,
        'stateMonF_e')


################################################################################

ei_net.simulate(options.time, printTime=True)
ei_net.endSimulation()


output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
        options.fileNamePrefix, options.job_num)
ei_net.saveAll(output_fname)

################################################################################

ei_net.printTimes()
