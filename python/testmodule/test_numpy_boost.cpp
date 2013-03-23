/*
 *   test_numpy_boost.cpp
 *
 *   Test numpy_boost.
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


#include "numpy_boost.hpp"

typedef numpy_boost<int, 1> intArr;

void func()
{

    int dims[] = {10};

    numpy_boost<double, 1> ret(dims);
    ret[0] = 1;
}

int main(void)
{
    func();

    return 0;
}
