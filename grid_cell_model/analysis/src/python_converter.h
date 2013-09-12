/*
 *   python_converter.h
 *
 *   Conversions Python <--> Blitz++.
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

#ifndef _PYTHON_CONVERTER_H
#define _PYTHON_CONVERTER_H

#include <Python.h>
#include <numpy/arrayobject.h>

/**
 * Provides a mapping from C++ datatypes to Numpy type
 * numbers.
 */
namespace detail {
    template<class T>
    class numpy_type_map {
        public:
            static const int typenum;
    };

    template<>
    const int numpy_type_map<float>::typenum = NPY_FLOAT;

    template<>
    const int numpy_type_map<std::complex<float> >::typenum = NPY_CFLOAT;

    template<>
    const int numpy_type_map<double>::typenum = NPY_DOUBLE;

    template<>
    const int numpy_type_map<int>::typenum = NPY_INT;

    template<>
    const int numpy_type_map<std::complex<double> >::typenum = NPY_CDOUBLE;

    template<>
    const int numpy_type_map<long double>::typenum = NPY_LONGDOUBLE;

    template<>
    const int numpy_type_map<std::complex<long double> >::typenum = NPY_CLONGDOUBLE;

    //template<>
    //const int numpy_type_map<boost::int8_t>::typenum = NPY_INT8;

    //template<>
    //const int numpy_type_map<boost::uint8_t>::typenum = NPY_UINT8;

    //template<>
    //const int numpy_type_map<boost::int16_t>::typenum = NPY_INT16;

    //template<>
    //const int numpy_type_map<boost::uint16_t>::typenum = NPY_UINT16;

    //template<>
    //const int numpy_type_map<boost::int32_t>::typenum = NPY_INT32;

    //template<>
    //const int numpy_type_map<boost::uint32_t>::typenum = NPY_UINT32;

    //template<>
    //const int numpy_type_map<boost::int64_t>::typenum = NPY_INT64;

    //template<>
    //const int numpy_type_map<boost::uint64_t>::typenum = NPY_UINT64;
}


class python_exception : public std::exception {
};


/**
 * Conver PyObject* to a PyArrayObject*.
 *
 * @param obj Python object. Any array-like object.
 * @return Numpy PyArrayObject*. The documentation for PyArrayFromObject is
 *      little cryptic, but this function should always increase the reference
 *      count of the **object that is returned**
 * @throw python_exception when the conversion fails
 */
template<class T, int N>
PyArrayObject* convertPyToNumpy(PyObject* obj)
{
    PyArrayObject* arr_obj;

    arr_obj = (PyArrayObject *) PyArray_FromObject(obj,
            detail::numpy_type_map<T>::typenum, N, N);

    if (arr_obj == NULL) {
        throw python_exception();
    }

    return arr_obj;
}


template<class T, int N>
blitz::Array<T,N> convertPyToBlitz(PyArrayObject* arr_obj)
{
    blitz::TinyVector<int,N> shape(0);
    blitz::TinyVector<int,N> strides(0);
    for (int i = 0; i < N; i++)
    {
        shape[i] = arr_obj->dimensions[i];
        strides[i] = arr_obj->strides[i]/sizeof(T);
    }
    return blitz::Array<T,N>((T*) arr_obj->data,shape,strides,
                             blitz::neverDeleteData);
}



/**
 * @file python/blitz_numpy.cc
 * @date Mon Sep 26 11:47:30 2011 +0200
 * @author Andre Anjos <andre.anjos@idiap.ch>
 *
 * @brief Automatic converters to-from python for blitz::Array's
 *
 * Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * Avoids the big number of warnings...
 */
static PyArrayObject* make_pyarray(int nd, npy_intp* dims, int type) {
  return (PyArrayObject*)PyArray_SimpleNew(nd, dims, type);
}

/**
 * Objects of this type bind blitz::Array<T,N> to numpy arrays. Your method
 * generates as output an object of this type and the object will be
 * automatically converted into a Numpy array.
 */
template <typename T, int N>
static PyObject* convertBlitzToPy(const blitz::Array<T, N>& tv)
{
    typedef typename blitz::Array<T,N> array_type;
    typedef typename blitz::TinyVector<int,N> shape_type;

    npy_intp dims[N];
    for (int i=0; i<N; ++i)
        dims[i] = tv.extent(i);

    PyArrayObject* retval = make_pyarray(N, dims,
            detail::numpy_type_map<T>::typenum);

    //wrap new PyArray in a blitz layer and then copy the data
    shape_type shape=0;
    for (int k=0; k<retval->nd; ++k)
        shape[k] = retval->dimensions[k];

    shape_type stride=0;
    for (int k=0; k<retval->nd; ++k)
        stride[k] = (retval->strides[k]/sizeof(T));

    array_type bzdest((T*)retval->data, shape, stride, blitz::neverDeleteData);
    bzdest = tv;

    return reinterpret_cast<PyObject*>(retval);
};

#endif // _PYTHON_CONVERTER_H
