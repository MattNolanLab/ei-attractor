#
#   colors.py
#
#   Color manipulations.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
import matplotlib.colors

class Colormap2D(matplotlib.colors.ListedColormap):

    def __init__(self, size, XColor, YColor, baseColor=[0, 0, 0],
            name='from_list'):
        self.XSize, self.YSize = size[0], size[1]
        self.baseColor = baseColor
        self.XColor    = np.asarray(XColor)
        self.YColor    = np.asarray(YColor)

        totalN = self.XSize * self.YSize
        colors = np.ndarray((totalN, 3))
        X, Y = np.meshgrid(np.linspace(0, 1, self.XSize), np.linspace(0, 1, self.YSize))
        it = 0
        for y in xrange(self.YSize):
            for x in xrange(self.XSize):
                colors[it, :] = baseColor + X[y, x]*self.XColor + \
                        Y[y,x]*self.YColor
                it += 1
        matplotlib.colors.ListedColormap.__init__(self, colors, name, totalN)




if __name__ == "__main__":
    import matplotlib.pyplot as plt
    xSize, ySize = (100, 100)
    size = xSize, ySize

    baseColor = (0, 0, 0)
    XColor    = (1, 0, 0)
    YColor    = (0, 1, 0)

    cmap2D = Colormap2D(size, baseColor, XColor, YColor, name='2dcmap')
    X, Y = np.meshgrid(range(size[0]), range(size[1]))
    plt.pcolor(X + Y*xSize, cmap=cmap2D)
    plt.show()
