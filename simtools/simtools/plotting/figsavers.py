'''Automatic figure savers.'''
from __future__ import absolute_import, print_function, division

from matplotlib.backends.backend_pdf import PdfPages


class MultiFigureSaver(object):
    '''Matplotlib figure saver that processes several images.

    This is an abstract class that defines basic functionality. Use some of the
    derived class ideally combined with dependency injection or some kind of
    factory.

    The concept is that the client gets this interface and simply calls savefig
    multiple times. Depending on the concrete implementation, either the
    figures will be saved as `one` combined PDF, or separate PDFs.
    '''
    def __init__(self, file_name, ext='pdf', start_cnt=0):
        self.file_name = file_name
        self.ext = ext
        self._cnt = 0

    def set_file_name(self, file_name):
        self.file_name = file_name
        self.reset()

    def reset(self, cnt_val=None):
        if cnt_val is not None:
            self._cnt = cnt_val

    def set_backend_params(self, **kwargs):
        self._backend_params = kwargs

    def savefig(self, fig):
        raise NotImplementedError()

    def close(self):
        pass


class PdfOutputSaver(MultiFigureSaver):
    '''Create a single-PDF multipage document, by calling the savefig()
    method.
    '''
    def __init__(self, file_name, ext='pdf'):
        super(PdfOutputSaver, self).__init__(file_name, ext, start_cnt=0)
        self._saver = None
        self._reset()

    def reset(self, cnt_val=None):
        self._reset()

    def _reset(self):
        self.close()
        if self.file_name is not None:
            fname = "%s.%s" % (self.file_name, self.ext)
            self._saver = PdfPages(fname)

    def savefig(self, fig):
        self._saver.savefig(fig, **self._backend_params)

    def close(self):
        if self._saver is not None:
            self._saver.close()


class SeparateMultipageSaver(MultiFigureSaver):
    '''Create PDF documents, that are numbered by a counter. This counter is
    appended to the file name when it is being saved.
    '''
    def __init__(self, file_name, ext='pdf'):
        super(SeparateMultipageSaver, self).__init__(file_name, ext,
                                                     start_cnt=0)

    def savefig(self, fig):
        fname = "%s_%d.%s" % (self.file_name, self._cnt, self.ext)
        fig.savefig(fname, **self._backend_params)
        self._cnt += 1


