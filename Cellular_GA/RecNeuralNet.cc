/**
 * Name:        RecNeuralNet.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Character recognition neural Net
 */

#include <sstream>
#include <iostream>
#include <fstream>

#include "RecNeuralNet.h"
#include "error.h"
#include "debug.h"

using namespace std;

/* ------------------------------------------------------------------------ */

RecNeuralNet::RecNeuralNet(unsigned numChars, const NeuralNetwork** cNets,
            const NeuralNetwork* codNet) :
    numChars(numChars)
{
    if (codNet->getNumInputs() != numChars)
    {
        ostringstream err;
        err << "An attempt to build RecNeuralNet, but coding layer input "
               "size doesn't match number of characters" << endl;
        throw NNException(err);
    }

    for (unsigned i = 0; i < numChars; i++)
    {
        if (cNets[i]->getNumOutputs() != 1)
        {
            ostringstream err;
            err << "An attempt to build RecNeuralNet, but char layer[" << i << "] "
                   "is not a one neuron output" << endl;
            throw NNException(err);
        }
    }

    charNets = new NeuralNetwork*[numChars];
    for (unsigned i = 0; i < numChars; i++)
        charNets[i] = cNets[i]->copy();

    codingNet = codNet->copy();
    charNetOuts = new double[numChars];
}

/* ------------------------------------------------------------------------ */

//RecNeuralNet::RecNeuralNet(const char* fileName)
//{
//    ifstream f(fileName);
//
//    if (!f.is_open())
//    {
//        ostringstream err;
//        err << "Could not open RecNeuralNetFile: " << fileName <<
//               " for reading" << endl;
//        throw NNException(err);
//    }
//
//    f >> numChars;
//
//    if (!f.good())
//    {
//        ostringstream err;
//        err << "Corrupted recNeuralNet file: " << fileName << endl;
//        throw NNException(err);
//    }
//
//    charNets = new NeuralNetwork*[numChars];
//
//    for (unsigned i = 0; i < numChars; i++)
//        charNets[i] = new NeuralNetwork(f);
//
//    codingNet = new NeuralNetwork(f);
//    charNetOuts = new double[numChars];
//}

/* ------------------------------------------------------------------------ */

RecNeuralNet::~RecNeuralNet()
{
    for (unsigned i = 0; i < numChars; i++)
        delete charNets[i];
    delete [] charNets;

    delete codingNet;
    delete []charNetOuts;
}

/* ------------------------------------------------------------------------ */

const double* RecNeuralNet::eval(const double* inputs)
{
    for (unsigned i = 0; i < numChars; i++)
        charNetOuts[i] = (charNets[i]->eval(inputs))[0];

    return codingNet->eval(charNetOuts);
}

/* ------------------------------------------------------------------------ */

void RecNeuralNet::save(ostream& s) const
{
    s << numChars << endl;

    for (unsigned i = 0; i < numChars; i++)
        charNets[i]->print(s);

    codingNet->print(s);
}

/* ------------------------------------------------------------------------ */


/* End of file RecNeuralNet.cc */
