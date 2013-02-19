/*
 *   PulsatingCurrentGenerator.h
 *
 *   Pulsating urrent generator
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


#ifndef PULSATINGCURRENTGENERATOR_H
#define PULSATINGCURRENTGENERATOR_H

#include <cmath>

#include "CurrentGenerator.h"

/**
 * Generate a current of type 0.5 * A * (1.0 + sin(2*pi*f*t + phi)).
 */
class PulsatingCurrentGenerator : public CurrentGenerator {

    public:

    /**
     * Construct the object.
     *
     * @param dt        Time step (ms)
     * @param startT    Time corresponding to the start of the simulation (ms)
     * @param startWhen When to start emitting pulsating current (ms)
     * @param A         Amplitude of the current (max)
     * @param f         Frequency (Hz)
     * @param phi       Phase (rad)
     */
    PulsatingCurrentGenerator(double dt, double startT, double startWhen, double A, double f, double phi) :
        CurrentGenerator(dt, startT), startWhen(startWhen), A(A), f(f*1e-3), phi(phi)
    {
        currTime = startT;
        currSteps = 0;
        I = getCurrentAtTime(startT);
    }


    double getCurrent() {
        return I;
    }


    /** Advance one time step **/
    void advance() {
        currSteps++;
        currTime = startT + currSteps*dt;
        if (currTime < startWhen) {
            I = .0;
        } else {
            I =  getCurrentAtTime(currTime);
        }
    }


    void reset() {
        currTime = startT;
        currSteps = 0;
        I = getCurrentAtTime(startT);
    }



    private:


    /** Return current for the specified time in ms.
     * 
     * @param t Specified time.
     */
    double getCurrentAtTime(double t) {
        return .5 * A * (1. + sin(2. * M_PI * f * currTime + phi));
    }

    double A;           //!< Amplitude (max) (pA)
    double f;           //!< Frequency (1/msec) Different to what is in the constructor
    double phi;         //!< Phase (rad)

    double currTime;    //!< Current simulation time (ms)
    double currSteps;   //!< Number of steps from startT
    double I;           //!< Current current (pA)

    double startWhen;   //!< When to start emitting pulsating current (ms)
};


#endif // PULSATINGCURRENTGENERATOR_H
