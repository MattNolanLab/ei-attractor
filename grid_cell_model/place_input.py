#
#   place_input.py
#
#   Place cell input simulation class.
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

from scipy.io import loadmat


class PlaceCellInput(object):

    def __init__(self, Ne_x, Ne_y, arenaSize, gridsep, gridCenter):

        self.sigma = 7;


        self.Ne_x = Ne_x
        self.Ne_y = Ne_y
        self.arenaSize = arenaSize

        self.gridsep_x = gridsep
        self.gridCenter = gridCenter
        self.gridsep_y = self.gridsep_x * np.sqrt(3);
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
        self.X= X
        self.Y= Y

    def getSheetInput(self, pos_x, pos_y):
        '''
        Return the sheet input for a given position of the animal.
        Rows are the X dimension.
        '''
        x = (pos_x - self.arenaSize)/self.dx
        y = (pos_y - self.arenaSize)/self.dx

        return self.arena[y-self.Ne_y/2:y+self.Ne_y/2, x-self.Ne_x/2:x+self.Ne_x/2]



if __name__=="__main__":
    Ne_x = 68
    Ne_y = 58
    arenaSize = 180.
    gridsep = arenaSize/3.
    gridCenter = [0, gridsep/2]
    
    pc = PlaceCellInput(Ne_x, Ne_y, arenaSize, gridsep, gridCenter)
    #pcolormesh(pc.X, pc.Y, pc.arena); axis('equal'); show()
    #pcolormesh(pc.arena); axis('equal'); show()
    
    
    vel_fname = '../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat'
    ratData = loadmat(vel_fname)
    pos_x = ratData['pos_x'].flatten()
    pos_y = ratData['pos_y'].flatten()
    
    for it in xrange(len(pos_x)):
        tmp = pc.getSheetInput(pos_x[it], pos_y[it])
        if len(tmp) != Ne_y or len(tmp[0]) != Ne_x:
            raise Exception()
    
        print len(tmp[0]), len(tmp)
        print it
    
