/*
 *   testmodule.c
 *
 *   Testing (and learning in fact) interfacing C with Python.
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
//#include <boost/python.hpp>
#include "numpy_boost.hpp"

typedef numpy_boost<double, 1> DoubleVec;

template <typename T>
T &sum(const T &p1, const T &p2) {
    size_t len = p1.shape()[0];
    T *res = new T(p1.shape());
    for (size_t i = 0; i < len; i++) {
        (*res)[i] = p1[i] + p2[i];
    }
    return *res;
}


static PyObject *
_test_sum(PyObject *self, PyObject *args)
{
    PyObject *in1;
    PyObject *in2;

    if (!PyArg_ParseTuple(args, "OO", &in1, &in2)) {
        return NULL;
    }


    try {
        DoubleVec v1(in1);
        DoubleVec v2(in2);

        if (v1.size() != v2.size()) {
            PyErr_SetString(PyExc_TypeError, "Vectors must be of the same length!");
            return NULL;
        }

        DoubleVec &ret = sum<DoubleVec>(v1, v2);

        Py_INCREF(ret.py_ptr());
        return ret.py_ptr();
    } catch (...) {
        std::cout << "test" << std::endl;
        return NULL;
    }
}




static PyMethodDef TestMethods[] = {

    {"sum", _test_sum, METH_VARARGS, "Sum two integers"},

    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
init_test(void)
{
    (void) Py_InitModule("_test", TestMethods);

    import_array();
}
