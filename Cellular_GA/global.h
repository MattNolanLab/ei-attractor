/**
 * Name:        global.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Global declarations
 */

#ifndef GLOBAL_H
#define GLOBAL_H

#include <cstdlib>
#include <cmath>
#include <cassert>

#include "debug.h"

/* MPI root node (gathers statistics) */
#define ROOT_NODE 0

/** Crossover types **/
#define C_RANDOMWALK 0
#define C_AREA       1


/**
 * Compute random number with exponential distribution on both sides.
 *
 * @param lambda Lambda parameter.
 */
inline double expRandom(double lambda)
{
    assert(lambda > 0);

    double x = rand() / (RAND_MAX + 1.0);
    double rnd = - 1/lambda * log(1 - x);

    assert(rnd != NAN);

    return rnd;
}

inline double uniformRandom()
{
    return rand() / (RAND_MAX + 1.0);
}

#endif /* GLOBAL_H */

/* End of file global.h */
