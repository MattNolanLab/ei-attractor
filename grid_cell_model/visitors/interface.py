'''Interface definitions for the visitors package.

.. currentmodule:: grid_cell_model.visitors.interface

Classes
-------

.. autosummary::

    Visitor
    DictDSVisitor
    EmptyDSVisitor
'''
from __future__ import absolute_import, print_function

from abc import ABCMeta, abstractmethod
import os

from ..data_storage.sim_models.ei import extractSpikes


class Visitor(object):
    '''
    An abstract visitor class.

    Normally, the base class of the Visitor design pattern must contain all the
    methods implemented in the derived classes. Due to duck typing, one does
    not need to declare the specific implementation methods of the visitor.
    '''
    __meta__ = ABCMeta

    @abstractmethod
    def _dummy(self):
        raise NotImplementedError()


    def getFigPath(self, dataFilePath, rootDir, r, c, trialNum, ext="pdf"):
        '''Return name of figure related to data being processed.

        Extract the data file name from the path given and return the name
        of the figure related to the data file name. Also, try to create the
        underlying directories (if the resulting path does not exist) so the
        caller can simply output the file with that name into the specified
        path.
        '''
        dataRootDir, dataFileName = os.path.split(dataFilePath)
        figDir = dataRootDir + "/"
        if rootDir is not None:
            figDir += rootDir
            if not os.path.exists(figDir):
                os.makedirs(figDir)
        figFileName = os.path.splitext(dataFileName)[0]
        figPath = "{figDir}/{figFileName}_r{r}_c{c}_tr{tr}.{ext}".format(
                figDir=figDir, figFileName=figFileName, r=r, c=c, tr=trialNum,
                ext=ext)
        return figPath


class DictDSVisitor(Visitor):
    '''
    Dictionary data set visitor.

    A visitor that takes a dictionary data set method of any kind and processes
    it. All the keys in the dictionary must be strings.
    '''

    @abstractmethod
    def visitDictDataSet(self, ds, **kw):
        '''
        Visit the dictionary data set, 'ds', and perform the specific operations
        (defined by the derived classes) on this data set.
        '''
        raise NotImplementedError()

    def getOption(self, data, optStr):
        '''Extract an option from a data dictionary'''
        return data['options'][optStr]

    def getNetParam(self, data, p):
        '''Extract a network parameter (p) from the data dictionary'''
        return data['net_attr'][p]


    def folderExists(self, d, nameList):
        '''
        Check if the folder at the end of the 'nameList' exists in dictionary
        'd'
        '''
        for name in nameList:
            if (name in d.keys()):
                d = d[name]
            else:
                return False
        return True

    def _checkAttrIsNone(self, attr, pName, data):
        '''
        Check if 'attr' is None. If yes, extract it from 'data' under the name
        'pName', otherwise return attr back.
        '''
        if (attr is None):
            return self.getOption(data, pName)
        else:
            return attr


    def _getSpikeTrain(self, data, monName, dimList):
        '''
        Return the senders and spike times from a monitor in the data
        dictionary

        Parameters
        ----------
        data : dict
            A data dictionary containing monName
        monName : str
            The name of the monitor
        dimList : list of string
            A  list of dimensions. This will be used to extract the total
            number of neurons from network parameters.
        output: tuple
            A tuple containing (senders, times, N). 'N' is the total number of
            neurons in the population.
        '''
        N = 1
        for dim in dimList:
            dimN = self.getNetParam(data, dim)
            N *= dimN
        mon = data[monName]
        senders, times = extractSpikes(mon)
        return senders, times, N


    @staticmethod
    def _getNeuronSpikeTrain(data, monName, n):
        '''
        Return a spike train for the particular neuron indexed by `n`.  If
        there are no spikes from the particular neuron, return an empty array.

        Parameters
        ----------
        data : dict-like object
            Data dictionary
        monName : string
            Name of the monitor
        n : int
            Neuron number.
        output : numpy array
            Spikes of the selected neuron
        '''
        mon = data[monName]
        senders, times = extractSpikes(mon)
        idx = (senders == n)
        return times[idx]


class EmptyDSVisitor(DictDSVisitor):
    '''
    A Dictionary date set visitor that will do nothing.
    '''
    def __init__(self):
        pass

