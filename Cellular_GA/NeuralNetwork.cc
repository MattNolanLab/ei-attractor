/**
 * Name:        NeuralNetwork.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural Network implementation.
 */

#include <cmath>
#include <cassert>

#include "debug.h"

#include "NeuralNetwork.h"
#include "global.h"

using namespace std;

const double NNFullLayer::lambda = 1.25;

/* ------------------------------------------------------------------------ */

NNLayer::~NNLayer() {}

/* ------------------------------------------------------------------------ */

NNFullLayer::NNFullLayer(unsigned items, NNLayer* prev, bool randomInit) :
    prevLayer(prev), items(items) 
{
    assert(prev != NULL);

    outputs = new double[items];

    /*
     * Initialize neuron weights.
     * Each neuron has incoming weights from all neurons from the previous
     * layer.
     */
    weights = new double*[items];
    unsigned prevItems = prevLayer->getNumItems();

    /*
     * Allocate weights. +1, the last weight is the bias
     */
    for (unsigned i = 0; i < items; i++)
        weights[i] = new double[prevItems + 1];

    if (randomInit)
        this->randomInit();
}

/* ------------------------------------------------------------------------ */

NNFullLayer::~NNFullLayer()
{
    delete []outputs;

    for (unsigned i = 0; i < items; i++)
        delete []weights[i];

    delete []weights;
}

/* ------------------------------------------------------------------------ */

void NNFullLayer::randomInit()
{
    unsigned numWeights = prevLayer->getNumItems() + 1;

    for (unsigned n = 0; n < items; n++)
    {
        for (unsigned i = 0; i < numWeights; i++)
        {
            double init = 2.0 * ::uniformRandom() - 1.0;
            weights[n][i] = init;
        }
    }
}

/* ------------------------------------------------------------------------ */

NNLayer* NNFullLayer::copy(NNLayer* prev) const
{
    assert(prev != NULL);

    if (prev->getNumItems() != prevLayer->getNumItems())
    {
        ostringstream err;
        err << "An attemt to copy a layer, but the new previous layer "
               "has different number of nodes." << endl;
        throw NNException(err);
    }

    NNFullLayer* aCopy = new NNFullLayer(items, prev, false);
    unsigned numWeights = prev->getNumItems() + 1;
    for (unsigned int node = 0; node < items; node++)
    {
        for (unsigned wi = 0; wi < numWeights; wi++)
        {
            aCopy->weights[node][wi] = weights[node][wi];
        }
    }

    return aCopy;
}

/* ------------------------------------------------------------------------ */

double NNFullLayer::computeNeuronPotential(unsigned n, const double* in)
{
    unsigned prevItems = prevLayer->getNumItems();

    double sum = 0.0;

    for (unsigned i = 0; i < prevItems; i++)
        sum += weights[n][i] * in[i];

    /* Apply bias */
    sum += weights[n][prevItems];

    return sum;
}

/* ------------------------------------------------------------------------ */

void NNFullLayer::biasWeights(unsigned n, double mRate)
{
    assert(prevLayer != NULL);
    assert(mRate >= 0 && mRate <= 1);

    static const double expLambda = 128.0;

    if (!(n < items))
    {
        ostringstream err;
        err << "Attempt to bias weights in non-existent neuron" << endl;
        throw NNException(err);
    }

    unsigned numWeights = prevLayer->getNumItems() + 1;
    double nMutW = mRate*numWeights;
    double bias;

    if (nMutW >= 10)
    {
        for (unsigned i = 0; i < nMutW; i++)
        {
            int wIndex = (int) (numWeights * ::uniformRandom());
            assert(wIndex >= 0 && (unsigned)wIndex < numWeights);

            bias = ::expRandom(expLambda);
            //cerr << bias << endl;
            if (::uniformRandom() < 0.5)
                bias = -bias;
            weights[n][wIndex] += bias;

            if ((fabs(weights[n][i])) > 1.5)
                multiplyWeights(n, 0.90);
        }
    }
    else
    {
        for (unsigned i = 0; i < numWeights; i++)
        {
            if (::uniformRandom() < mRate)
            {
                bias = ::expRandom(expLambda);
                if (::uniformRandom() < 0.5)
                    bias = -bias;
                weights[n][i] += bias;

                if ((fabs(weights[n][i])) > 1.5)
                    multiplyWeights(n, 0.9);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */

void NNFullLayer::multiplyWeights(unsigned n, double m)
{
    assert(n < items);
    assert(m >= 0);

    unsigned nWeights = prevLayer->getNumItems() + 1;
    for (unsigned i = 0; i < nWeights; i++)
        weights[n][i] *= m;
}

/* ------------------------------------------------------------------------ */

//void NNFullLayer::changeWeight(unsigned n, unsigned w, double val)
//{
//    assert(prevLayer != NULL);
//    assert(n < items);
//    assert(w < prevLayer->getNumItems() + 1);
//
//    weights[n][w] = val;
//}

/* ------------------------------------------------------------------------ */

void NNFullLayer::combineUniform(const NNLayer* l)
{
    const NNFullLayer* layer = dynamic_cast<const NNFullLayer* >(l);

    assert(items == layer->items);
    assert(prevLayer->getNumItems() == layer->prevLayer->getNumItems());

    unsigned prevItems = prevLayer->getNumItems() + 1;
    unsigned changed = 0;
   
    for (unsigned i = 0; i < items; i++)
    {
        for (unsigned p = 0; p < prevItems; p++)
        {
            double rnd = ::uniformRandom();
            
            if (rnd < 0.5)
            {
                weights[i][p] = layer->getWeight(i, p);
                changed++;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */

void NNFullLayer::averageWeights(const NNLayer* l)
{
    const NNFullLayer* layer = dynamic_cast<const NNFullLayer* >(l);

    assert(items == layer->items);
    assert(prevLayer->getNumItems() == layer->prevLayer->getNumItems());

    unsigned prevItems = prevLayer->getNumItems() + 1;
   
    for (unsigned i = 0; i < items; i++)
    {
        for (unsigned p = 0; p < prevItems; p++)
        {
            weights[i][p] = (weights[i][p] + layer->getWeight(i, p)) / 2;
        }
    }
}

/* ------------------------------------------------------------------------ */

double NNFullLayer::computeActivationFunc(double ksi)
{
    /* Sigmoidal activation function */
    double f = 1.0 / (1.0 + exp(-lambda * ksi));
    assert(f >= 0);
    //cerr << f << endl;
    return f;
}

/* ------------------------------------------------------------------------ */

const double* NNFullLayer::eval()
{
    const double* prevOutputs = prevLayer->eval();

    /* Compute each neuron's output */
    for (unsigned n = 0; n < items; n++)
    {
        double ksi = computeNeuronPotential(n, prevOutputs);
        outputs[n] = computeActivationFunc(ksi);
        if (outputs[n] > 1.0)
            cerr << "!!!!!Warning: neuron output > 1.0: " << outputs[n] << endl;
    }

    return outputs;
}

/* ------------------------------------------------------------------------ */

void NNFullLayer::print(std::ostream& s) const
{
    unsigned wNum = prevLayer->getNumItems() + 1;

    s << wNum << endl;

    /* Neuron per line */
    for (unsigned i = 0; i < items; i++)
    {
        for (unsigned w = 0; w < wNum; w++)
        {
            s << weights[i][w] << " ";
        }
        s << endl;
    }
}

/* ------------------------------------------------------------------------ */


NeuralNetwork::NeuralNetwork(unsigned numLayers, const unsigned* lSizes,
        bool random) :
    numLayers(numLayers)
{
    if (!checkNumLayers(numLayers))
    {
        ostringstream err;
        err << "Attempt to create a perceptron feedforward neural"
               " network with less than 2 layers." << endl;
        throw NNException(err);
    }

    layers = new NNLayer*[numLayers];

    /** Create the first - input layer **/
    layers[0] = new NNInputLayer(lSizes[0]);

    /* Create network layers - the first one is the bottom one */
    for (unsigned i = 1; i < numLayers; i++)
        layers[i] =  new NNFullLayer(lSizes[i], layers[i - 1], random);

}

/* ------------------------------------------------------------------------ */

NeuralNetwork::NeuralNetwork(unsigned numLayers, const NNLayer** l) :
    numLayers(numLayers)
{
    assert(l != NULL);
    assert(dynamic_cast<const NNInputLayer* >(l[0]) != NULL);

    layers = new NNLayer*[numLayers];

    /* Make a copy of the first (input) layer */
    layers[0] = l[0]->copy(NULL); 

    for (unsigned i = 1; i < numLayers; i++)
        layers[i] = l[i]->copy(layers[i - 1]);
}

/* ------------------------------------------------------------------------ */

NeuralNetwork::NeuralNetwork(const NeuralNetwork& n) :
    numLayers(n.numLayers)
{
    layers = new NNLayer*[numLayers];

    /* Make a copy of the input layer */
    layers[0] = n.layers[0]->copy(NULL);

    /* Make copues of Full layers */
    for (unsigned i = 1; i < numLayers; i++)
        layers[i] = n.layers[i]->copy(layers[i - 1]);
}

/* ------------------------------------------------------------------------ */

NeuralNetwork::~NeuralNetwork()
{
    for (unsigned i = 0; i < numLayers; i++)
        delete layers[i];

    delete []layers;
}


/* ------------------------------------------------------------------------ */

const double* NeuralNetwork::eval(const double* nnInputs)
{
    /* Feed the first layer */
    (static_cast<NNInputLayer* >(layers[0]))->feed(nnInputs);

    /* Eval the last layer. Each layer will call eval to its prevLayer */
    return layers[numLayers - 1]->eval();
}

/* ------------------------------------------------------------------------ */

NeuralNetwork* NeuralNetwork::copy() const 
{
    return new NeuralNetwork(*this);
}

/* ------------------------------------------------------------------------ */

void NeuralNetwork::combineUniform(const NeuralNetwork* n)
{
    assert(n != NULL);

    for (unsigned l = 1; l < numLayers; l++)
        layers[l]->combineUniform(n->layers[l]);
}

/* ------------------------------------------------------------------------ */

void NeuralNetwork::averageWeights(const NeuralNetwork* n)
{
    assert(n != NULL);

    for (unsigned l = 1; l < numLayers; l++)
        layers[l]->averageWeights(n->layers[l]);
}

/* ------------------------------------------------------------------------ */

void NeuralNetwork::print(ostream& s)
{
    /* Print summary */
    s << numLayers << endl;

    /* Print weights -- do not print input layer */
    for (unsigned l = 0; l < numLayers; l++)
        layers[l]->print(s);
}

/* ------------------------------------------------------------------------ */

void NeuralNetwork::randomBias(double mRate)
{
    for (unsigned l = 1; l < numLayers; l++)
    {
        for (unsigned n = 0; n < layers[l]->getNumItems(); n++)
        {
            layers[l]->biasWeights(n, mRate);
        }
    }
}

/* ------------------------------------------------------------------------ */

/* End of file NeuralNetwork.cc */
