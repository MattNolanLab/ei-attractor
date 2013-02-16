/**
 * Name:        NNChromosome.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural network chromosome.
 */

#ifndef NNCHROMOSOME_H
#define NNCHROMOSOME_H

#include <cstdlib>
#include <ostream>

#include "Chromosome.h"
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
class NNChromosome : public Chromosome
{

  public:
    /**
     * Constructor that creates a random chromosome with specified length.
     */
    NNChromosome(unsigned numLayers, const unsigned* lSizes,
            const NNTrainData* td);


    /**
     * Make a new chromosome, copying existing layers.
     *
     * @param numLayers Number of layers.
     * @param layers Layers array of pointers. Copies of these layers will
     *   be made
     * @param td Train data.
     */
    NNChromosome(unsigned numLayers, const NNLayer** layers,
            const NNTrainData* td);

    ~NNChromosome();

    /**
     * Evaluate the chromosome's fitness value and return it.
     */
    double getFitness();


    /**
     * Get neural network
     */
    NeuralNetwork* getNN() const { return nn; }

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
     * Print chromosome into stream using special train data.
     */
    void print(std::ostream &s, const NNTrainData* td = NULL) const;

    /**
     * Print as debug info into specified stream.
     */
    void printDebug(std::ostream &s)
    {
        s << "No debug information about NNChromosome" << std::endl;
    }

    static void printNumberOfInstances(std::ostream& s)
    {
        s << "NNChromosome: number of instances: " << instances << std::endl;
    }

    void setNNTrainData(const NNTrainData* td)
    {
        assert(td != NULL);
        trainData = td;
    }

  protected:
    /**
     * Copy constructor
     */
    NNChromosome(const NNChromosome& c);

    /** Check if the topologies of the two specified chromosomes match **/
    bool checkTopologyMatch(const NNChromosome* c1, const NNChromosome* c2) const;
    
    /** Combine layers in a neural network **/
    const NNLayer** combineLayers(const NNChromosome* chr) const;

    /** Neural network for this chromosome **/
    NeuralNetwork* nn;

    /** A trainset for this chromosome **/
    const NNTrainData* trainData;

    bool fitnessValid;
    double fitness;

    /** Number of instances **/
    static int instances;

    static const double lambda;

    /** 
     * Treshold for the error sum to be considered OK
     */
    static const double errSumTreshold;

    /** Treshold for the bit error to be considered the right bit **/
    static const double errBitTreshold;


};

#endif /* NNCHROMOSOME_H */


/* end of file NNChromosome.h */
