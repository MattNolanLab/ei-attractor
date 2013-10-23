#
#   settings.py
#
#   Global settings for the thesis figures.
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
from matplotlib import rc

# Global plotting settings
fontsize = 11
letterSize = 14

# Matplotlib
rc('pdf', fonttype=42)
rc('mathtext', default='regular')
rc('font', size=fontsize)


def setFName(fname, outputRoot='output'):
    return "{0}/{1}".format(outputRoot, fname)

