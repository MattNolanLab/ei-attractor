/*
 *   _spikesmodule.cpp
 *
 *   Advanced spike analysis
 *
 *       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
 *       
 *       This program is free software: you can redistribute it and/or modify
 *       it under the terms of the GNU General Public License as published by
 *       the Free Software Foundation, either version 3 of the License, or
 *       (at your option) any later version.
 *       
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *       GNU General Public License for more details.
 *       
 *       You should have received a copy of the GNU General Public License
 *       along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <exception>

#include <Python.h>
#include <numpy/arrayobject.h>
#include <blitz/array.h>

#include "definitions.h"
#include "spikes.h"
#include "python_converter.h"

extern "C" {


static PyObject *
_spikes_spike_time_diff(PyObject *self, PyObject *args)
{
    PyObject *in1 = NULL;
    PyObject *in2 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *arr2 = NULL;


    if (!PyArg_ParseTuple(args, "OO", &in1, &in2)) {
        return NULL;
    }


    try {
        arr1 = convertPyToNumpy<double, 1>(in1);
        arr2 = convertPyToNumpy<double, 1>(in2);
        DblVector train1 = convertPyToBlitz<double, 1>(arr1);
        DblVector train2 = convertPyToBlitz<double, 1>(arr2);
        DblVector td = spikes::spike_time_diff(train1, train2);

        PyObject *ret = convertBlitzToPy<double, 1>(td);

        Py_XDECREF(arr1);
        Py_XDECREF(arr2);
        return ret;
    } catch (python_exception) {
        Py_XDECREF(arr1);
        Py_XDECREF(arr2);
        return NULL;
    }
    Py_RETURN_NONE;
}


static PyMethodDef SpikesMethods[] = {

    {"spike_time_diff", _spikes_spike_time_diff, METH_VARARGS, "Distribution of"
       " spike time differences"},

    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
init_spikes(void)
{
    (void) Py_InitModule("_spikes", SpikesMethods);

    import_array();
}

} // extern "C"
