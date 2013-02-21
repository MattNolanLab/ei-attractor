/*
 *   CurrentGenerator.h
 *
 *   Generic interface for a current generator
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


#ifndef CURRENTGENERATOR_H
#define CURRENTGENERATOR_H

/**
 * Generate current of different properties.
 *
 * This is an abstract class that defines how a model requests a current value
 * during the simulation.
 * It allows to advance the current simulation time by one or more steps of dt
 * and then request the value. The underlying implementation may or may not
 * perform integration.
 *
 * NOTES: integration, but also function
 *
 */
class CurrentGenerator {
    public:

    /**
     * Construct the current generator.
     *
     * @param dt     Simulation dt
     * @param startT Time corresponding to the start of the simulation.
     */
    CurrentGenerator(double dt, double startT) :
        dt(dt), startT(startT)
    {}


    /**
     * Get current at the current time position.
     * @return Current at the current time position
     */
    virtual double getCurrent() = 0;


    /**
     * Advance one time step.
     */
    virtual void advance() = 0;


    /**
     * Advance several time steps.
     *
     * @param n Number of time steps to advance.
     */
    virtual void advance(int n) {
        for (int i = 0; i < n; i++)
            advance();
    }


    /** Reset the state to when t==startT. **/
    virtual void reset() = 0;

    /** Get start time **/
    double getStartTime() const { return startT; }
    /** Set start time **/
    void   setStartTime(double t) { startT = t; }

    protected:

    double dt;      //!< Simulation time step
    double startT;  //!< Start time in ms
    
};


#endif // CURRENTGENERATOR_H

