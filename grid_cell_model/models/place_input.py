'''Place cell input simulation class.

.. currentmodule:: grid_cell_model.models.place_cells

.. note::

    This is an older code which was used in [PASTOLL2013]_. It is only used to
    define connection weights from place cells to E cells in the current
    version of the model. Otherwise use the NEST place cell model which is
    defined in
    :meth:`~grid_cell_model.models.gc_net_nest.NestGridCellNetwork.setPlaceCells`

'''
import numpy as np

from scipy.io import loadmat


class PlaceCellInput(object):
    '''
    Place cell input class. One can define an arena size (circular) and then
    probe the template for an input to a population of neurons that should
    center the attractor bump of activity to a desired location. This class can
    also be used to set up connection weights from place cells to E cells.

    This can be used to reset the bump to a correct position
    '''

    def __init__(self, Ne_x, Ne_y, arenaSize, gridsep, gridCenter, fieldSigma=7):
        '''
        Create a place cell input template. Here we work with a twisted torus
        topology with dimensions X x Y = 1 x sqrt(3)/2.
        Ne_x/y      size of the neural sheet
        arenaSize   Arena radius (cm)
        gridSep     Distance between grid field peaks (cm)
        gridCenter  Offset of the center of grid fields (cm)
        fieldSigma  Sigma of the gaussians used to produce the template (cm)
        '''

        self.sigma = fieldSigma


        self.Ne_x = Ne_x
        self.Ne_y = Ne_y
        self.arenaSize = arenaSize

        self.gridsep_x = gridsep
        self.gridCenter = gridCenter
        self.gridsep_y = self.gridsep_x * np.sqrt(3) # Assuming twisted torus
        self.dx = self.gridsep_x/Ne_x;

        #import pdb; pdb.set_trace()

        X, Y = np.meshgrid(np.arange(arenaSize*2/self.dx) * self.dx, np.arange(arenaSize*2/self.dx) * self.dx);

        X = 1. * (X - arenaSize)
        Y = 1. * (Y - arenaSize)

        X_mod = np.abs(np.mod(X - self.gridsep_x/2 - gridCenter[0], self.gridsep_x) - self.gridsep_x/2)
        Y_mod = np.abs(np.mod(Y - self.gridsep_y/2 - gridCenter[1], self.gridsep_y) - self.gridsep_y/2)
        arena1 = np.exp(-(X_mod**2 + Y_mod**2)/2/self.sigma**2)

        shift_x = self.gridsep_x*np.cos(np.pi/3)
        shift_y = self.gridsep_x*np.sin(np.pi/3)
        arena2 = np.exp(-( (X_mod - shift_x)**2 + (Y_mod-shift_y)**2)/2/self.sigma**2)

        self.arena = arena1 + arena2
        self.arena = self.arena / np.max(self.arena)
        self.X= X
        self.Y= Y


    def getSheetInput(self, pos_x, pos_y):
        '''
        Return the sheet input for a given position of the animal.
        Rows are the Y dimension.
        '''
        x = (pos_x - self.arenaSize)/self.dx
        y = (pos_y - self.arenaSize)/self.dx

        return self.arena[y-self.Ne_y/2:y+self.Ne_y/2, x-self.Ne_x/2:x+self.Ne_x/2]



if __name__=="__main__":
    from matplotlib.pyplot import *
    Ne_x = 34
    Ne_y = 30
    arenaSize = 180.
    gridsep = 70            # cm
    gridCenter = [0, 0]

    pc = PlaceCellInput(Ne_x, Ne_y, arenaSize, gridsep, gridCenter)
    pcolormesh(pc.X, pc.Y, pc.arena); axis('equal'); show()
    pcolormesh(pc.arena); axis('equal'); show()
    pcolormesh(pc.getSheetInput(.0, .0)); show()


    vel_fname = '../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat'
    ratData = loadmat(vel_fname)
    pos_x = ratData['pos_x'].flatten()
    pos_y = ratData['pos_y'].flatten()

    print "Testing dimensions of input..."
    for it in xrange(len(pos_x)):
        tmp = pc.getSheetInput(pos_x[it], pos_y[it])
        if len(tmp) != Ne_y or len(tmp[0]) != Ne_x:
            raise Exception()

        #print len(tmp[0]), len(tmp)
        #print it
    print "OK\n"


    print "Testing grid cell heat map"
    neuronNums = [10, 20, 10*Ne_x+10, 20*Ne_x+10]

    dx = 3      # cm
    xedges = np.linspace(-arenaSize/2, arenaSize/2, arenaSize/dx + 1)
    yedges = np.linspace(-arenaSize/2, arenaSize/2, arenaSize/dx + 1)
    for neuronNum in neuronNums:
        print "  neuron no. " + str(neuronNum)
        hist = np.ndarray((len(xedges), len(yedges)), dtype=object)
        avgHist = np.ndarray((len(xedges), len(yedges)))
        nInput = []
        for it in xrange(len(pos_x)):
            x = pos_x[it]
            y = pos_y[it]

            x_bin = np.floor(x/dx) + np.floor(len(xedges)/2)
            y_bin = np.floor(y/dx) + np.floor(len(yedges)/2)
            if (x_bin < 0 or y_bin < 0):
                raise Exception()

            if hist[y_bin, x_bin] == None:
                hist[y_bin, x_bin] = []

            hist[y_bin, x_bin].append(pc.getSheetInput(x, y).ravel()[neuronNum])


        for x_it in range(len(xedges)):
            for y_it in range(len(yedges)):
                if hist[y_it, x_it] == None:
                    avgHist[y_it, x_it] = 0
                else:
                    avgHist[y_it, x_it] = np.mean(hist[y_it, x_it])


        X, Y = np.meshgrid(xedges, yedges)
        pcolormesh(X, Y, avgHist)
        xlabel('X position (cm)')
        ylabel('Y position (cm)')
        title('Neuron no. ' + str(neuronNum))
        savefig('grid_field_heat_map_{0:04}.png'.format(neuronNum))


