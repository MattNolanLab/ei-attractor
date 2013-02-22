/*
 *   gridcells_exceptions.h
 *
 *   Grid cells model exceptions.
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

#ifndef GRIDCELLS_EXCEPTIONS_H
#define GRIDCELLS_EXCEPTIONS_H

#include "exceptions.h"


class SizeException : public nest::KernelException
{
    std::string msg;

    public:
    //! @param detailed error message
    SizeException() :
        KernelException("SizeException"),
        msg()
    {}


    SizeException(std::string msg) :
        KernelException("SizeException"),
        msg(msg)
    {}


    ~SizeException() throw () {}

    std::string message();
};


#endif // GRIDCELLS_EXCEPTIONS_H
