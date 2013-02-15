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
import numpy as np

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

    def hdf_filename(self):
        '''Output file name'''
        return self.f.filename


    def saveDictionary(self, dict_name, d):
        '''
        Serialize a dictionary into the output file under a group name 'dict_name'.
        Each item in the dicionary must be an array type or string
        '''
        dict_grp = self.f.create_group(dict_name)
        for key, val in d.iteritems():
            try:
                tmp = np.array(val)
                dict_grp.create_dataset(str(key), tmp.shape, tmp.dtype)
                dict_grp[str(key)][...] = tmp
            except Exception as e:
                print "Could not save a dictionary item '" + str(key) + "' to HDF"
                #print "Exception message" + str(e)

    def saveArray(self, name, arr):
        '''Save an array to the output file, under name "name"'''
        self.f.create_dataset(name, data=arr)


    def saveSpikes(self, grp_name, obj_arr):
        '''
        Save spike times for each neuron from an object array, into a separate
        group grp_name. Each neuron will be labeled by its index conveted into
        a string identifier
        '''
        spike_grp = self.f.create_group(grp_name)
        for n_it in xrange(len(obj_arr)):
            spike_grp.create_dataset(str(n_it), data=np.array(obj_arr[n_it]))


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


# Test this module
if __name__ == "__main__":
    output_dir = "."
    filenamePrefix = ''
    job_num = 1111
    trial_it = 10
    dump_it = 23
    dict = {"Iext": 123, 'Ivel': 10e-3, 'none': None}
    array1 = np.arange(1000)
    array2 = np.arange(100) * 0.01

    od = OutputDump(output_dir, filenamePrefix, job_num, trial_it, dump_it)
    output_fname = od.filename

    od.saveDictionary('options', dict)
    od.saveArray('array1', array1)
    od.saveArray('array2', array2)


    od.close()

