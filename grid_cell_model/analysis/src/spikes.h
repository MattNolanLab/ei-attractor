/*
 *   spikes.h
 *
 *   Spike analysis python independent functions.
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

#ifndef SPIKES_H
#define SPIKES_H

#include <cmath>
#include <blitz/array.h>

namespace spikes {

/**
 * Compute the distribution of the spike time differences of all pairs of
 * spikes of two spike trains (1D arrays) (t2 - t1). The order thus matters: if
 * t1 precedes t2, then the result will be positive.
 *
 * @param train1 First spike train.
 * @param train2 Second spike train.
 * @return An array of spike time differences for each spike pair.
 */
template<typename T>
blitz::Array<T, 1>
spike_time_diff(const blitz::Array<T, 1> &train1, const blitz::Array<T, 1>
        &train2) {
    int sz1 = train1.size();
    int sz2 = train2.size();

    int szRes = sz1 * sz2;
    blitz::Array<T, 1> res(szRes);
    int resIdx = 0;
    for (int tIdx1 = 0; tIdx1 < sz1; tIdx1++) {
        for (int tIdx2 = 0; tIdx2 < sz2; tIdx2++) {
            res(resIdx) = train2(tIdx2) - train1(tIdx1);
            resIdx++;
        }
    }

    return res;
}



} // namespace spikes

#endif // SPIKES_H
