#
#   test_data_storage.py
#
#   Unit tests for the data_storage package
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
'''
.. currentmodule: test.test_data_storage

The :mod:`~test.test_data_storage` module defines a set of classes for unit
testing of the most important functions and classes of this package.  The
module is currently based on the unittest module.

A list of currently supported tests:

.. autosummary::
'''

import unittest
import numpy as np
import collections

from data_storage import DataStorage

notImplMsg = "Not implemented"

##############################################################################
class TestHDF5Storage(unittest.TestCase):

    def getItem(self, d, key):
        return d[key]

    def test_basic_types(self):
        ds = DataStorage.open('test_basic_types.h5', 'w')

        test_int      = 124
        test_float    = 0.1
        test_1d_array = np.arange(1000000)
        test_2d_array = np.random.rand(100, 100)
        test_list     = [1, 2, 3, 4]
        test_dict     = dict(
                int=123,
                float=111.1,
                list=test_list)
                
        ds['int']      = test_int
        ds['float']    = test_float
        ds['1d_array'] = test_1d_array
        ds['2d_array'] = test_2d_array
        ds['list']     = test_list
        ds['dict']     = test_dict

        # Test return types
        self.assertTrue(isinstance(ds['int'], int))
        self.assertTrue(isinstance(ds['float'], float))
        self.assertTrue(isinstance(ds['1d_array'], np.ndarray))
        self.assertTrue(isinstance(ds['2d_array'], np.ndarray))
        self.assertTrue(isinstance(ds['list'], collections.MutableSequence))
        self.assertTrue(isinstance(ds['dict'], collections.MutableMapping))

        # Test values of basic types
        self.assertEqual(ds['int'], test_int)
        self.assertEqual(ds['float'], test_float)
        self.assertTrue(np.all(ds['1d_array'] == test_1d_array))
        self.assertTrue(np.all(ds['2d_array'] == test_2d_array))
        self.assertEqual(ds['list'], test_list)
        self.assertEqual(ds['dict'], test_dict)

        ds.flush()

        # Test delete
        del(ds['int'])
        self.assertRaises(KeyError, self.getItem, ds, 'int')
        del(ds['float'])
        self.assertRaises(KeyError, self.getItem, ds, 'float')
        del(ds['1d_array'])
        self.assertRaises(KeyError, self.getItem, ds, '1d_array')
        del(ds['2d_array'])
        self.assertRaises(KeyError, self.getItem, ds, '2d_array')
        del(ds['list'])
        self.assertRaises(KeyError, self.getItem, ds, 'list')
        del(ds['dict'])
        self.assertRaises(KeyError, self.getItem, ds, 'dict')

        ds.close()

        # ds should be empty now
        ds = DataStorage.open('test_basic_types.h5', 'r')
        self.assertEqual(len(ds), 0)

        ds.close()


    def test_lists(self):
        def appendListAndTest(ds, key, test_l, item):
            ds[key].append(item)
            test_l.append(item)
            self.assertEqual(test_l, ds[key])

        ds = DataStorage.open('test_lists.h5', 'w')

        d1 = {
                "hola" : [10, 20, 30],
                "str" : "This is a test string"}

        test_list = [1, 2, d1]
        ds['list'] = test_list
        ds.close()

        ds = DataStorage.open('test_lists.h5', 'r+')
        appendListAndTest(ds, 'list', test_list, 10)
        appendListAndTest(ds, 'list', test_list, 23.5)
        appendListAndTest(ds, 'list', test_list, [1, 2, 3])
        appendListAndTest(ds, 'list', test_list, dict(a=10, b=[1, 2, 3]))
        ds.close()



    def test_iterator(self):
        ds = DataStorage.open('test_iterator.h5', 'w')

        test_list = list(np.arange(100))
        ds['list'] = test_list 

        it1 = 0
        for val1 in ds['list']:
            self.assertEqual(val1, test_list[it1])  
            it1 += 1

            it2 = 0
            for val2 in ds['list']:
                self.assertEqual(val2, test_list[it2])
                it2 += 1

        ds.close()



    def test_empty_arr(self):
        arr = np.array([])

        ds = DataStorage.open('test_empty_arr.h5', 'w')
        ds['empty'] = arr
        ds.close()

        ds = DataStorage.open('test_empty_arr.h5', 'r')
        self.assertTrue(np.all(arr == ds['empty']))
        self.assertEqual(len(ds['empty']), 0)


    def test_chained_getter(self):
        test_dict     = dict(
                int=123,
                float=111.1,
                list=[1, 2, 3, dict(
                    a='blabla',
                    b=10,
                    c=np.random.rand(10))])

        ds = DataStorage.open('test_chained_getter.h5', 'w')
        ds['nested'] = test_dict

        ds_nested = ds['nested']
        self.assertEqual(ds_nested.getItemChained(('int',)), test_dict['int'])
        self.assertEqual(ds_nested.getItemChained(('list', 0)),
                test_dict['list'][0])
        self.assertEqual(ds_nested.getItemChained(('list', 3, 'a')),
                test_dict['list'][3]['a'])
        
        # The same test but with a list as an index
        self.assertEqual(ds_nested.getItemChained(['list', 3, 'a']),
                test_dict['list'][3]['a'])

        ds.close()


    def test_chained_setter(self):
        def testChain(ds, keyList, testValue):
            val = ds
            for key in keyList[0:-1]:
                val = val[key]
                self.assertTrue(isinstance(val, collections.MutableMapping))
            self.assertEqual(val[keyList[-1]], testValue)

        keyList = ['a', 'b', 'c', 'd']
        testValue = [1, 2, 3, 4]
        ds = DataStorage.open('test_chained_setter.h5', 'w')
        

        # Initial test
        ds.setItemChained(keyList, testValue)
        testChain(ds, keyList, testValue)

        # Over-write test
        newTestValue = 'Over-write test string'
        ds.setItemChained(keyList, newTestValue)
        testChain(ds, keyList, newTestValue)

        # Do not overwrite if overwriteLast is False
        noOverWriteTestValue = 'This should not be written'
        ds.setItemChained(keyList, noOverWriteTestValue, overwriteLast=False)
        testChain(ds, keyList, newTestValue)

        # ['a', 'b', 'another', 'xxx']
        otherKeyList = ['a', 'b', 'another', 'xxx']
        otherTestValue = [1, 2, 3, 'ahoy']
        ds.setItemChained(otherKeyList, otherTestValue)
        testChain(ds, otherKeyList, otherTestValue)
        testChain(ds, keyList, newTestValue) # Should not be overwritten

        # Single item in keyList
        singleList = ['single']
        ds.setItemChained(singleList, otherTestValue)
        testChain(ds, singleList, otherTestValue)

        # Test using getItemChained()
        self.assertEqual(ds.getItemChained(keyList), newTestValue)
        self.assertEqual(ds.getItemChained(otherKeyList), otherTestValue)
        self.assertEqual(ds.getItemChained(singleList), otherTestValue)
