/**
 * Name:        Chromosome.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Chromosome abstract class
 */

#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <cstdlib>
#include <ostream>


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
#   include "mpi.h"
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */


/**
 * Chromosome class.
 * This is a pure abstract class that provides interface for GA chromosomes.
 */
class Chromosome
{
  public:

    Chromosome()
    {
        local = true;
    }

    virtual ~Chromosome() = 0;

    /**
     * Evaluate chromosome and return fittness
     */
    virtual double getFitness() = 0;

    /**
     * Make a copy of myself.
     */
    virtual Chromosome *copy() const = 0;

    /**
     * Do crossover with Chromosome c
     *
     * @param c Chromosome to do crossover with
     * @return New offspring pointer
     */
    virtual Chromosome *crossover(const Chromosome *c) const = 0;

    /**
     * Mutate myself.
     *
     * @param mRate mutation rate in the <0, 1> range
     */
    virtual void mutate(double mRate) = 0;

    /**
     * Print chromosome into a stream
     */
    virtual void print(std::ostream &s) const = 0;

    /**
     * Print as debug info int specified stream.
     */
    virtual void printDebug(std::ostream &s) = 0;

    /**
     * Defines if the Chromosome is in the local grid.
     */
    bool isLocal() const { return local; };


    /**
     * Set if the Chromosome is local or global.
     */
    void setLocal(bool l) { local = l; };


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    /**
     * Get read-only data for this chromosome for use with MPI_Send().
     * The data are handled by output function parameters.
     *
     * @param buf   Buffer data pointer. It will be allocated read-only, don't
     *              delete it.
     * @param count Number of data types to send (for use with MPI_Send)
     * @param type  MPI data type
     */
    virtual void
    getMPISendData(const void **buf, int *count, MPI_Datatype *type) const = 0;


    /**
     * Replace my genetic code with that received by MPI functions. 
     * This function doesn't check if the proper chromosome has been received,
     * it's the user's responsibility.
     *
     * @param buf  Received data buffer. If the return value is true, this buffer 
     *             will be handled properly by the chromosome and should not be
     *             further manipulated by the user.
     * @param size Number of MPI data elements.
     * @return True if the operation succeedd. If this function returns false,
     *         it indicates some fatal error, because the data received doesn't
     *         correspond to the chromosome's properties.
     */
    virtual bool replaceWithMPIChromosome(void *buf, int size) = 0;

    /**
     * Create new chromosome from the MPI data.
     * @param buf MPI buffer.
     * @param siz Number of MPI data elements.
     * @return New Chromosome.
     */
    virtual Chromosome *getNewFromMPI(void *buf, int siz) const = 0;


    /**
     * Get MPI communication properties for this chromosome (for use with
     * MPI_Send() and MPI_Recv().
     *
     * @param buf   Receive buffer. Allocated by the Chromosome itself. Must be
     *              properly handled by the user or returned back in the
     *              ReplaceWithMPIChromosome() function.
     * @param count Number of MPI_Datatype units.
     * @param type  MPI_Datatype for use.
     */
    virtual void getMPIProps(void **buf, int *count, MPI_Datatype *type) = 0;

#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */


  protected:

    /**
     * Is the fitness valid? This flag is used for quick fitness evaluation.
     */
    bool fitnessValid;

    /** Get random double number in the range <0, 1> **/
    double getRandomDouble() const 
    {
        return rand() / (RAND_MAX + 1.0);
    }

  private:
    bool local;
};


#endif /* CHROMOSOME_H */


/* end of file Chromosome.h */
