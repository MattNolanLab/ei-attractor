/*
 *   VelocityInputGenerator.h
 *
 *   Velocity input current generator.
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

#ifndef VELOCITYINPUTGENERATOR_H
#define VELOCITYINPUTGENERATOR_H

#include <vector>

#include "CurrentGenerator.h"
#include "gridcells_exceptions.h"


/**
 * Velocity inputs structure.
 *
 * Provides several data members for animal movement simulation.
 */
class VelocityInputs {


    public:
    VelocityInputs() {
        dt = .0;
    }


    VelocityInputs(const std::vector<double>& x, const std::vector<double>& y, double dt) :
        dt(dt)
    {
        setPos(x, y);
    }

    void setPos(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size())
            throw nest::DimensionMismatch(x.size(), y.size());

        ratPosX = x;
        ratPosY = y;
        calculateVelocity();
    }

    const std::vector<double>& getVelX() const { return ratVelX; }
    const std::vector<double>& getVelY() const { return ratVelY; }
    double get_dt() const { return dt; }

    int getXSize() const { return ratVelX.size(); }
    int getYSize() const { return ratVelY.size(); }

    const std::vector<double>& getPosX() const { return ratPosX; }
    const std::vector<double>& getPosY() const { return ratPosY; }


    private:
    void calculateVelocity() {
        for (int i = 0; i < ratPosX.size() - 1; i++) {
            ratVelX[i] = (ratPosX[i+1] - ratPosX[i-1]) / (dt * 1e-3);
            ratVelY[i] = (ratPosY[i+1] - ratPosY[i-1]) / (dt * 1e-3);
        }
    }

    std::vector<double> ratPosX; //!< Animal positions X (cm)
    std::vector<double> ratPosY; //!< Animal positions Y (cm)
    std::vector<double> ratVelX; //!< Animal velocity X (cm/s)
    std::vector<double> ratVelY; //!< Animal velocity Y (cm/s)
    double dt;                   //!< Time step for the positions (ms)
};


/**
 * Preferred direction vector
 */
struct PrefDirs {
    double x;   //!< X component
    double y;   //!< Y component

    /** Initialises the preferred direction to [x, y] == [.0, .0] **/
    PrefDirs() {
        x = .0;
        y = .0;
    }

    PrefDirs(double x, double y) :
        x(x), y(y)
    {}
};



/**
 * Generates velocity inputs for each neuron, based on the provided animal
 * trajectory data.
 *
 * The parameters genStart defines when to start reading the velocity input data
 * and generate currents when asked. So the first value will be generated at
 * genStart time, while zeros will be generated before that time.
 *
 * The velocity generator generates current that is a dot product of the
 * specified preferred direction and the current velocity vector (derived from
 * the positional data).
 */
class VelocityInputGenerator : public CurrentGenerator {

    public:


    /**
     * Create the object.
     *
     * @param dt       Simulation time step (ms)
     * @param startT   Simulation start time (ms)
     * @param genStart Start of input generation (ms)
     *
     * @throw BadProperty
     */
    VelocityInputGenerator(double dt, double startT, double genStart, const PrefDirs pd) :
        CurrentGenerator(dt, startT), genStart(genStart), prefDirs(pd)
    {
        veldtSteps = (int) (velInputs.get_dt() / dt);
        if (veldtSteps != veldtSteps)
            throw nest::BadProperty("The velocity input time step must be a multiple of the simulation time step!");

        if (veldtSteps == 0)
            throw nest::BadProperty("Velocity input time step must be >= the simulation time step!");

        stepNow = 0;
        nextChange = stepNow + veldtSteps;
        velInputIt = 0;
    }


    /** Get input current for the current step. But only when we are active **/
    double getCurrent() {
        double currT = stepNow*dt + startT;
        if (currT >= genStart) {
            double x = velInputs.getVelX()[velInputIt] * prefDirs.x;
            double y = velInputs.getVelY()[velInputIt] * prefDirs.y;
            return x + y;
        } else {
            return .0;
        }
    }

    /** Advance one simulation time step **/
    void advance()
    {
        stepNow++;
        if (stepNow >= nextChange) {
            nextChange += veldtSteps;
            velInputIt++;
        }
    }


    /** Reset the state to time startT **/
    void reset()
    {
        stepNow = 0;
        nextChange = stepNow + veldtSteps;
        velInputIt = 0;
    }


    /**
     * Setup velocity inputs for all objects
     * 
     * All object share the same set of velocity inputs
     *
     * @param vi Velocity inputs
     */
    void setVelocityInputs(const VelocityInputs& vi) {
        if (vi.getXSize() == 0)
            throw SizeException("Size of the X part of the velocity vector must be > 0");
        if (vi.getYSize() == 0)
            throw SizeException("Size of the Y part of the velocity vector must be > 0");

        velInputs = vi;
    }


    /** Get velocity inputs (of all objects) **/
    const VelocityInputs& getVelocityInputs() const {
        return velInputs;
    }


    /** Get preferred directions **/
    const PrefDirs& getPrefDirs() const {
        return prefDirs;
    }

    private:

    double genStart;   //!< When to start generating an input
    PrefDirs prefDirs; //!< Preferred direction of this generator

    long   stepNow;
    long   velInputIt;
    int    veldtSteps;
    double nextChange;


    static VelocityInputs velInputs; //!< Velocity inputs, shared by all objects

};


#endif // VELOCITYINPUTGENERATOR_H
