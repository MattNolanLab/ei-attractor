/**
 * Name:        RecNeuralNet.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Characters recognition neural net.
 */

#ifndef RECNEURALNET_H
#define RECNEURALNET_H

#include <iostream>


#include "NeuralNetwork.h"
#include "debug.h"

class RecNeuralNet
{
  public:
    /**
     * Create neural network.
     *
     * @param numChars Number of charNets
     * @param charNets One character recognition NNs
     * @param codingNet Binary coding and filtering NN
     */
    RecNeuralNet(unsigned numChars, const NeuralNetwork** charNets,
            const NeuralNetwork* codingNet);

    /**
     * Load net from file.
     */
    RecNeuralNet(const char* fileName);

    /**
     * Create Net from file.
     *
     * @param s File input stream.
     */
    RecNeuralNet(const std::istream& s);

    ~RecNeuralNet();

    const double* eval(const double* inputs);

    void save(std::ostream& s) const;


  private:
    NeuralNetwork** charNets;
    NeuralNetwork*  codingNet;
    unsigned numChars;
    double* charNetOuts;
};

#endif /* RECNEURALNET_H */

/* End of file RecNeuralNet.h */
