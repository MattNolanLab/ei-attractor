'''Define basic data structures for pipeline processing during the simulations.'''

class PipelineStage(object):
    '''A base class for callable pipeline stage functors.

    This class can be configured and instances can be called inside a pipeline.
    The class must define/override the __call__ method.
    '''
    def __call__(self, data_in):
        raise NotImplementedError


class PipelineData(object):
    '''A pipeline data object.

    An object to manage data passed through the pipeline transparently. It
    contains one attribute: ``items``, which contains the actual data. ``items``
    has a dictionary interface for now.
    '''
    def __init__(self):
        self.items = {}


class Pipeline(object):
    '''A pipeline that passes a data object between different consecutive
    processing stages.
    '''
    def __init__(self):
        self._data = None
        self._stages = []

    def append(self, stage):
        '''Append a pipeline stage.

        Parameters
        ----------
        stage : PipelineStage
            The processing stage to append. In principle this could be any
            callable object that accepts the :class:`~PipelineData`` instance.
        '''
        self._stages.append(stage)

    @property
    def num_stages(self):
        '''Number of stages in the pipeline.'''
        return len(self._stages)

    def run(self, data):
        '''Run ``data`` through the pipeline.

        Note that this process works on ``data`` **in situ**. If anything gets
        changed during the pipeline, it will be overwritten in the data object.

        Parameters
        ----------
        data : PipeLineData
            Data object that will be passed on to the first stage of the
            pipeline.

        Returns
        -------
        data : PipelineData
            Data after the last stage of the pipeline has been applied to it.
        '''
        for stage in self._stages:
            data = stage(data)
            if not isinstance(data, PipelineData):
                raise TypeError('All pipeline stages must strictly return a '
                                'value of type '
                                'simtools.pipelines.PipelineData.')
        return data

