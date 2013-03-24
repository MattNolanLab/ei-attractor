/*
 *   signal.h
 *
 *   Signal analysis python independent functions.
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

#ifndef _SIGNAL_H
#define _SIGNAL_H

#include <cmath>
#include <blitz/array.h>

namespace sig {

/**
 * 
 *
 * 
 */
template<typename T>
T *correlation_function(const T &v1, const T &v2, unsigned lag_start, unsigned lag_end) {
    size_t sz = min(v1.size(), v2.size());
    size_t szRes = lag_end - lag_start + 1;
    T *res = new T(sz);

    int i = 0;
    for (int lag = lag_start; i <= lag_end; i++) {
        int s = max(0, -lag);
        int e = min(sz - lag - 1, sz - 1);
        res(i) = dot(v1(Range(s, e)), v2(Range(s + lag, e+lag)));
        i++;
    }

    return *res;
}

} // namespace sig

#endif // _SIGNAL_H
