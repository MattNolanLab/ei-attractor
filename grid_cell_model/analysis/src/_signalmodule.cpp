/*
 *   _signalmodule.c
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
#include <blitz/array.h>

#include "signal.h"
#include "python_converter.h"

extern "C" {

typedef blitz::Array<double, 1> DblVector;

static PyObject *
_signal_correlation_function(PyObject *self, PyObject *args)
{
    PyObject *in1;
    PyObject *in2;
    int lag_start;
    int lag_end;


    if (!PyArg_ParseTuple(args, "OOii", &in1, &in2, &lag_start, &lag_end)) {
        return NULL;
    }

    try {
        DblVector v1 = convertPyToBlitz<double, 1>(in1);
        DblVector v2 = convertPyToBlitz<double, 1>(in2);
        //std::cout << v1 << v2;

        DblVector c = sig::correlation_function(v1, v2, lag_start, lag_end);
        return convertBlitzToPy<double, 1>(c);

    } catch (python_exception) {
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
