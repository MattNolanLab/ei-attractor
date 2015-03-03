'''Package initialisation.'''
__all__ = ['DimensionException', 'SubmitError']

class DimensionException(Exception):
    def __init__(self, method, msg=None):
        self.method = method
        self.msg = msg

    def __str__(self):
        retstr = "Wrong dimensions in {0}".format(self.method)
        if self.msg is not None:
            retstr += "Additional message: " + self.msg
        return retstr


class SubmitError(Exception):
    pass



