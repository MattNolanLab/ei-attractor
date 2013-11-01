/*
 *   _signalmodule.cpp
 *
 *   Advanced signal analysis
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

#include "signal.h"
#include "python_converter.h"
#include "definitions.h"

extern "C" {


static PyObject *
_signal_correlation_function(PyObject *self, PyObject *args)
{
    PyObject *in1 = NULL;
    PyObject *in2 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *arr2 = NULL;
    int lag_start;
    int lag_end;


    if (!PyArg_ParseTuple(args, "OOii", &in1, &in2, &lag_start, &lag_end)) {
        return NULL;
    }


    try {
        arr1 = convertPyToNumpy<double, 1>(in1);
        arr2 = convertPyToNumpy<double, 1>(in2);
        DblVector v1 = convertPyToBlitz<double, 1>(arr1);
        DblVector v2 = convertPyToBlitz<double, 1>(arr2);
        DblVector c = sig::correlation_function(v1, v2, lag_start, lag_end);

        PyObject *ret = convertBlitzToPy<double, 1>(c);

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



static PyMethodDef SignalMethods[] = {

    {"correlation_function", _signal_correlation_function, METH_VARARGS, "Correlation function"},

    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
init_signal(void)
{
    (void) Py_InitModule("_signal", SignalMethods);

    import_array();
}

} // extern "C"
