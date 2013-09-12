#! /usr/bin/env python
import numpy as np
from data_storage import DataStorage


d = DataStorage.open('list.h5', 'w')
d1 = {"hola" : [10, 20, 30], "blabla": np.arange(20)}
d['list'] = [1, 2, d1]
d.close()

print("\nReopening the file\n")
d = DataStorage.open('list.h5', 'r+')
l = d['list']
print l
l.append(10)
print l
d['list_append'] = l
d['list'] = l

