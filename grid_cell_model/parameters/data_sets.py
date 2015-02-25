'''Data sets base classes module.'''


class DataSet(object):
    '''
    The data set object interface
    '''

    def __init__(self):
        raise NotImplementedError()

    @property
    def data(self):
        '''Container that contains data.'''
        raise NotImplementedError()

    @property
    def parameters(self):
        raise NotImplementedError()

    def visit(self, visitor, **kw):
        '''Apply the visitor to the data set.'''
        raise NotImplementedError()


class DictDataSet(DataSet):
    '''A data set that holds a dictionary structure.'''
    def __init__(self, dataDict):
        self._d = dataDict

    @property
    def data(self):
        return self._d

    def visit(self, v, **kw):
        v.visitDictDataSet(self, **kw)
