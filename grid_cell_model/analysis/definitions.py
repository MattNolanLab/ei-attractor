'''
Type definitions common to the whole package.
'''

class Range(object):
    '''A range object.

    Defines low bound (start) and high bound (stop). Pretty much useless
    because python has a native range object.
    '''
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


    def __str__(self):
        return "[start: " + self.start + ', stop: ' + self.stop
