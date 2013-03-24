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

//#include "signal.h"
#include "python_converter.h"


static PyObject *
_signal_dot(PyObject *self, PyObject *args)
{
    PyObject *in1;
    PyObject *in2;

    if (!PyArg_ParseTuple(args, "OO", &in1, &in2)) {
        return NULL;
    }


    try {
        blitz::Array<double, 1> A1 = convertPyToBlitz<double, 1>(in1);
        blitz::Array<double, 1> A2 = convertPyToBlitz<double, 1>(in2);
        std::cout << A1 << A2;
        std::cout << dot(A1, A2) << std::endl;
    } catch (python_exception) {
        return NULL;
    }

    Py_RETURN_NONE;
}




static PyMethodDef SignalMethods[] = {

    {"dot", _signal_dot, METH_VARARGS, "Dot product."},

    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
init_signal(void)
{
    (void) Py_InitModule("_signal", SignalMethods);

    import_array();
}
