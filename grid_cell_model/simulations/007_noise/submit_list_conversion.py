#! /usr/bin/env python
#
#   submit_list_conversion.py
#
#   Submit a list index format conversion script
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
import sys
import os

scriptPath = '../../data_storage/list_conversion.py'
numFiles = 100

if (len(sys.argv) != 2):
    print("Usage:\n{0} h5ListFile".format(sys.argv[0]))
fileList = list(open(sys.argv[1], 'r'))

def concatFileNames(l):
    ret = ''
    for name in l:
        if (name[-1] == '\n'):
            name = name[0:-1]
        ret += name + ' '
    return ret

it = 0
for fileIdx in xrange(0, len(fileList), numFiles):
    print(it)
    cmd = "qsub -N conv{0} -P inf_ndtc -o list_log -l h_rt=0:30:00 ./cluster_submit.sh {1} {2}"
    #cmd = "python {1} {2} > list_log/conv{0}"
    fileNames = concatFileNames(fileList[fileIdx:fileIdx+numFiles])
    cmd = cmd.format(it, scriptPath, fileNames)
    print cmd
    os.system(cmd)
    it += 1
