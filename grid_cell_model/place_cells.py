#
#   place_cells.py
#
#   Place cell simulation class.
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
from matplotlib.pyplot import *



class PlaceCells(object):
    '''
    A class to simulate place cells. This base class simply takes as an argument
    the number of place cells and the distribution of preferred positions and
    widths of the fields. The fields are then simulated as Gaussian functions of
    a position of the animal. The output is an "instantaneous" firing rate of
    each place cell, given a position.
    '''

    def __init__(self, N, maxRates, centers, widths):
        '''
        '''
        self.N          = N
        self.maxRates   = np.array(maxRates)
        self.centers    = np.array(centers)
        self.widths     = np.array(widths)

        if (self.maxRates.shape != (self.N, 2)):
            raise Exception('maxRates must be an array with dimensions (N, 2)')

        if (self.centers.shape != (self.N, 2)):
            raise Exception('centers must be an array with dimensions (N, 2)')

        if (self.widths.shape != (self.N, 2)):
            raise Exception('widths must be an array with dimensions (N, 2)')

    
    def getFiringRates(self, pos):
        '''
        Return a vector of firing rates for each place cell, given position
        'pos'
        '''
        return maxRates * np.exp(- ((pos[0] - centers[:, 0])**2 + (pos[1] -
            centers[:, 1])**2) / (2.0*widths**2)))


    def remap(self, envN):
        '''
        Remap the centers and possibly firing rates to the environment number
        envN
        '''
        pass



class UniformBoxPlaceCells(PlaceCells):
    '''
    Generate a uniform distribution of place cells in a 2D environment of a
    rectangular shape.
    '''

    def __init__(self, boxSize, N, maxRates, widths):
        '''
        N        A tuple containing number of place cells in each dimension (X, Y)
        maxRates An array of size Nx*Ny specifying max. firing rate in the place
                    field.
        widths   An array of widths of place fields (they are circular)
        '''
        self.boxSize = boxSize

        if (boxSize[0] <= 0 or boxSize[1] <= 0):
            raise Exception('boxSize dimenstions must be positive!')

        cx           = np.linspace(0, self.boxSize[0], N[0])
        cy           = np.linspace(0, self.boxSize[1], N[1])
        ctr_x, ctr_y = np.meshgrid(px, py)
        centers = np.hstack((cx.flatten(), cy.flatten())).T
        PlaceCells.__init__(self, N[0]*N[1], maxRates, centers, widths)


    def remap(self, envN):
        pass



if __name__ == '__main__':
    boxSize = (1.0, 1.0)
    N = (10, 10)
    maxRates = 15 # Hz
    widths = 0.1

    PC = UniformBoxPlaceCells(boxSize, N, maxRates, widths)
    print PC.centers

