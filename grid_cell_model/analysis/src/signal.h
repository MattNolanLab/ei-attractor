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
 * Lag-restricted correlation function. 
 *
 * Compute the correlation function (CF) of two blitz vectors (1D Arrays), with
 * lags that are user-defined.
 *
 * @param v1        First vector
 * @param v2        Second vector
 * @param lag_start Starting lag value
 * @param lag_en    End lag value (included in the CF)
 * @return A new blitz 1D array with size lag_end - lag_start + 1
 *
 * @note It is the responsibility of the user to pass correct values of
 *       lag_start and lag_end as these are not tested for boundaries.
 * @note The result for both arrays being empty is undefined.
 */
template<typename T>
blitz::Array<T, 1>
correlation_function(const blitz::Array<T, 1> &v1, const blitz::Array<T, 1>
        &v2, int lag_start, int lag_end) {
    size_t sz = std::min(v1.size(), v2.size());
    size_t szRes = lag_end - lag_start + 1;
    blitz::Array<T, 1> res(szRes);

    int i = 0;
    for (int lag = lag_start; lag <= lag_end; lag++) {
        int s = std::max(0, -lag);
        int e = std::min(sz - lag - 1, sz - 1);
        res(i) = dot(v1(blitz::Range(s, e)), v2(blitz::Range(s + lag, e+lag)));
        i++;
    }

    return res;
}

} // namespace sig

#endif // _SIGNAL_H
