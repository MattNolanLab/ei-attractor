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

#include <Python.h>

static PyObject *
test_sum(PyObject *self, PyObject *args)
{
    int a1, a2;

    if (!PyArg_ParseTuple(args, "ii", &a1, &a2)) {
        return NULL;
    }

    return Py_BuildValue("i", a1 + a2);
}


static PyMethodDef TestMethods[] = {

    {"sum", test_sum, METH_VARARGS, "Sum two integers"},

    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
inittest(void)
{
    (void) Py_InitModule("test", TestMethods);
}
