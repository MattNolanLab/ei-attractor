/**
 * Name:        OneMaxChromosome.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: OneMax chromosome class
 */

#ifndef ONEMAXCHROMOSOME_H
#define ONEMAXCHROMOSOME_H

#include <cstdlib>
#include <ostream>

#include "SimpleBitArray.h"
#include "Chromosome.h"


/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
#   include "mpi.h"
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */



/**
 * Chromosome class for OneMax.
 */
class OneMaxChromosome : public Chromosome
{

  public:
    /**
     * Constructor that creates a random chromosome with specified length.
     */
    OneMaxChromosome(unsigned int length);

    /**
     * Create a chromosome from an array.
     *
     * @param d Data
     * @param s Size of the chromosome
     */
    OneMaxChromosome(unsigned char *d, int s);

    /**
     * Copy constructor
     */
    OneMaxChromosome(const OneMaxChromosome& c);

    ~OneMaxChromosome();

    /**
     * Evaluate the chromosome's fitness value and return it.
     */
    double getFitness();

    /**
     * Make a copy of myself.
     */
    Chromosome *copy() const;

    /**
     * Do crossover with Chromosome c
     *
     * @param c Chromosome to do crossover with
     * @return New offspring pointer
     */
    Chromosome *crossover(const Chromosome *c) const;

    /**
     * Mutate myself.
     *
     * @param mRate mutation rate in the <0, 1> range
     */
    void mutate(double mRate);

    /**
     * Print chromosome into a stream
     */
    void print(std::ostream &s) const;

    /**
     * Print as debug info int specified stream.
     */
    void printDebug(std::ostream &s)
    {
        array->printDebug(s);
    }

    static void printNumberOfInstances(ostream& s)
    {
        s << "OneMaxChromosome: number of instances: " << instances << endl;
    }


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED

    /**
     * See Chromosome description for this function.
     */
    void getMPISendData(const void **buf, int *count, MPI_Datatype *type) const ;


    /**
     * See Chromosome description for this function.
     */
    bool replaceWithMPIChromosome(void *buf, int size);


    /**
     * See Chromosome description for this function.
     */
    void getMPIProps(void **buf, int *count, MPI_Datatype *type);

    /**
     * See Chromosome description for this function.
     */
    Chromosome *getNewFromMPI(void *buf, int siz) const;
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */

  private:

    /** Actual bitfield **/
    SimpleBitArray *array;

    bool fitnessValid;
    double fitness;

    /** Number of instances **/
    static int instances;

};

#endif /* ONEMAXCHROMOSOME_H */


/* end of file OneMaxChromosome.h */
