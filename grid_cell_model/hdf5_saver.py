#
#   hdf5_saver.py
#
#   HDF5 wrapper for saving simulation data to HDF format files 
#   Requires the h5py module
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
import h5py

__all__ = ['OutputDump']


class OutputDump(object):
    '''
    Creates a HDF5 output dump and allows to write datasets into this file.
    '''

    def __init__(self, output_dir, fileNamePrefix, job_num, trial_it, dump_it):
        self.output_dir = output_dir
        self.filenamePrefix = fileNamePrefix
        self.job_num = job_num
        self.trial_it = trial_it
        self.dump_it = dump_it

        self.filenameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"
        self.filename = self.filenameTemp.format(output_dir, fileNamePrefix,
                job_num, trial_it, dump_it)

        self.f = h5py.File(self.filename + '.h5', 'w')

    def filename(self):
        '''Output file name'''
        return self.f.filename


    def saveDictionary(self, d, dict_name):
        '''
        Serialize a dictionary into the output file under a group name 'dict_name'.
        Each item in the dicionary must be an array type or string
        '''
        dict_grp = self.f.create_group(dict_name)
        for key, val in d.iteritems():
            if val is not None:
                tmp = np.array(val)
                dict_grp.create_dataset(str(key), tmp.shape, tmp.dtype)
                dict_grp[str(key)][...] = tmp

    def saveArray(self, name, arr):
        '''Save an array to the output file, under name "name"'''
        self.f.create_dataset(name, data=arr)


    def close(self, cleanPreviousDump=False):
        self.f.close()
        if cleanPreviousDump:
            self.cleanPrevDump()

    def cleanPrevDump(self):
        rmFilename = self.filenameTemp.format(output_dir, fileNamePrefix,
                job_num, trial_it, dump_it - 1)
        try:
            os.remove(rmFileName)
        except:
            print "Couldn't remove previous dump file: " + rmFileName


