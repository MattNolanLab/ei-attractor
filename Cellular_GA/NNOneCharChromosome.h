/**
 * Name:        NNOneCharChromosome.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural network chromosome - one character recognition
 */

#ifndef NNONECHARCHROMOSOME_H
#define NNONECHARCHROMOSOME_H

#include <cstdlib>
#include <ostream>

#include "Chromosome.h"
#include "NNChromosome.h"
#include "NeuralNetwork.h"
#include "NNTrainData.h"


/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
#   include "mpi.h"
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */



/**
 * Chromosome class for a neural network.
 */
class NNOneCharChromosome : public NNChromosome
{

  public:
    /**
     * Constructor that creates a random chromosome with specified length.
     * 
     * @param inputSize Number of input neurons
     * @param NNNodes Number of neurons in the hidden layer.
     * @param symbol Training symbol
     * @param td Train data.
     */
    NNOneCharChromosome(unsigned inputSize, unsigned NNNodes, const NNTrainData* td,
            unsigned char symbol);



    /**
     * Evaluate the chromosome's fitness value and return it.
     */
    double getFitness();

    /**
     * Crossover with another chromosome (must be NNOneCharChromosome) and
     * produce offspring.
     */
    Chromosome* crossover(const Chromosome *c) const;

    /**
     * Make a copy of myself.
     */
    Chromosome *copy() const;


    /**
     * Print chromosome into a stream
     */
    void print(std::ostream &s) const;

    /**
     * Print chromosome into stream using special train data.
     */
    void print(std::ostream &s, const NNTrainData* td = NULL) const;

    /**
     * Print as debug info into specified stream.
     */
    void printDebug(std::ostream &s)
    {
        s << "No debug information about NNOneCharChromosome" << std::endl;
    }

    void setNNTrainData(const NNTrainData* td)
    {
        assert(td != NULL);
        trainData = td;
    }

  private:
    /**
     * Copy constructor
     */
    NNOneCharChromosome(const NNOneCharChromosome& c);

    /**
     * Create new chromosome using existing layers (they are copied)
     *
     * @param l  New layers.
     * @param td Train data.
     * @param s  Training symbol.
     */
    NNOneCharChromosome(const NNLayer** l, const NNTrainData* td,
            unsigned char s);

    /**
     * Get number of nodes in the neural network - obscure way.
     */
    unsigned* getLSizes(unsigned nInputs, unsigned nNodes)
    {
        static unsigned sizes[3] = {nInputs, nNodes, 1};
        return sizes;
    }

    /** Training symbol - ASCII **/
    const unsigned char symbol;

    /** Debug **/
    unsigned lastBadBit;
};

#endif /* NNONECHARCHROMOSOME_H */


/* end of file NNOneCharChromosome.h */
