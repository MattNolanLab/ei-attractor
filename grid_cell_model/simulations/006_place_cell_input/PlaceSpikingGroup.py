#
#   PlaceSpikingGroup.py
#
#   A Brian NeuronGroup that produces Poisson place cell spikers.
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
import numpy as np

import brian
from brian import *

from place_cells import *


class PlaceSpikingGroup(brian.PoissonGroup):
    '''
    A class that implements Poisson spiking group, where each neuron is a place
    cell with rate dependent on the position of an animal in an environment.
    '''

    def __init__(self, PlaceCellModel, trajectory, traj_dt, clock):
        '''
        PlaceCellModel      Reference to PlaceCells object
        trajectory          A tuple containing evenly sampled animal 2D trajectory
        traj_dt             A dt for trajectory (sec)
        '''
        self._PC        = PlaceCellModel
        self.trajectory = trajectory
        self.traj_dt    = traj_dt
        self.clock      = clock


        ratesFunc = None
        PoissonGroup.__init__(self, self._PC.N, ratesFunc, self.clock) 


    def update(self):
        '''
        Override the update method so that we can simply update the rates our
        own way, i.e. using place cell generator
        '''
        traj_it = self.clock.t // self.traj_dt
        self._S[0, :] = self._PC.getFiringRates(self.trajectory[traj_it])
        NeuronGroup.update(self)


    def getPlaceCells(self):
        '''
        Return an object of a PlaceCells class
        '''
        return self._PC


# Tests
if __name__ == '__main__':
    boxSize = (121, 200)
    N = (10, 10)
    totalSz = N[0]*N[1]
    maxRates = 15
    widths = 40

    T = 10*second
    traj_dt = 20*ms
    traj = np.zeros((T/traj_dt + 1, 2))
    
    PC = UniformBoxPlaceCells(boxSize, N, maxRates, widths, random=False)
    SG = PlaceSpikingGroup(PC, traj, traj_dt)

    net = Network(SG)

    spikeMon = SpikeMonitor(SG)
    net.add(spikeMon)

    net.run(T)

    raster_plot(spikeMon)
    show()

