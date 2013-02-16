/**
 * Name:        ChromosomeFactory.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Abstract class for chromosome creation.
 */

#ifndef CHROMOSOMEFACTORY_H
#define CHROMOSOMEFACTORY_H

#include "NNTrainData.h"
#include "Chromosome.h"
#include "OneMaxChromosome.h"
#include "NNChromosome.h"

/**
 * Chromosome factory. Provides abstract interface for chromosome creation.
 */
class ChromosomeFactory
{
  public:
    virtual Chromosome *getChromosome() = 0;
};


/**
 * OneMax Chromosome factory. Produces a binary chromosome with a defined
 * lenght.
 */
class OneMaxChromosomeFactory : public ChromosomeFactory
{
  public:

    /**
     * Constructor.
     *
     * @param cLength Chromosome length, in bits.
     */
    OneMaxChromosomeFactory(unsigned int cLength) : length(cLength) {};

    /**
     * Produce a chromosome, initialized randomly.
     *
     * @return Produced chromosome.
     */
    Chromosome *getChromosome()
    {
        return new OneMaxChromosome(length);
    }

    /**
     * Return the specified length of the chromosomes.
     */
    unsigned int getChromosomeLength()
    {
        return length;
    }

  private:
    /** Chromosome length, in bits **/
    unsigned int length;
};


/**
 * Neural Net Chromosome factory. Produces a neural network chromosome.
 */
//class NNChromosomeFactory : public ChromosomeFactory
//{
//  public:
//
//    /**
//     * Constructor.
//     *
//     * @param numLayers Number of layers.
//     * @param lSizes An array of size for each layer.
//     */
//    NNChromosomeFactory(unsigned numLayers, const unsigned* lSizes,
//            const NNTrainData& trainData) :
//        numLayers(numLayers), lSizes(lSizes), trainData(trainData) { }
//
//    /**
//     * Produce a chromosome, initialized randomly.
//     *
//     * @return Produced chromosome.
//     */
//    Chromosome *getChromosome()
//    {
//        return new NNChromosome(numLayers, lSizes, &trainData);
//    }
//
//  private:
//    /** Chromosome length, in bits **/
//    unsigned numLayers;
//    const unsigned* lSizes;
//    const NNTrainData& trainData;
//};


#endif /* CHROMOSOMEFACTORY_H */

/* End of file */
