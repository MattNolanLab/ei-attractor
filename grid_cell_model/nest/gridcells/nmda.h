#ifndef NMDA_H
#define NMDA_H

#include <cmath>

/**
 * Computes the magnesium conductance multiplier. See Jahr & Stevens (1990).
 *
 * @param V Membrane voltage (mV)
 * @param C_Mg Mg2+ concentration (mM)
 */
inline
double nmda_mg_multiplier(double V, double C_Mg)
{
    return 1. / (1. + C_Mg * exp(-0.062*V) / 3.57);
}


#endif // NMDA_H
