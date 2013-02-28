#
#   stat.py
#
#   Statistical analysis packages for GridCells
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


__all__ = ['circularVariance']


def circularvariance(x, range):
    ''' Returns circular moments of a vector of circular variable x, defined in
    the 'range' (this will be mapped to a circle)'''
    c = np.exp(1j*2*np.pi*x/range)
    avg = np.mean(c)
    theta_avg = np.angle(avg)
    theta_var = 1 - np.abs(avg)
    return (theta_avg, theta_var)
    


