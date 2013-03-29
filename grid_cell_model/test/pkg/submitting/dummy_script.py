#!/usr/bin/env python
#
#   dummy_script.py
#
#   Dummy simulation for testing submitters.
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
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--job_num",     type="int")
parser.add_option("--output_data", type="int")
parser.add_option("--output_dir",  type="str")
o, args = parser.parse_args()


f = open('{0}/dummy_output_{1}.txt'.format(o.output_dir, o.output_data), 'w')
f.write(str(o.output_data))
f.close()
