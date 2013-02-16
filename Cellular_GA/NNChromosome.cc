/**
 * Name:        NNChromosome.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural Network chromosome implementation
 */


#include <cmath>

#include <map>

#include "NNChromosome.h"
#include "debug.h"
#include "global.h"

using namespace std;

const double NNChromosome::errSumTreshold = 0.80;
const double NNChromosome::errBitTreshold = 0.20;


/* ------------------------------------------------------------------------ */
int NNChromosome::instances = 0;
const double NNChromosome::lambda = 1.0;

NNChromosome::NNChromosome(unsigned numLayers, const unsigned* lSizes,
        const NNTrainData* td) :
    Chromosome(), trainData(td)
{
    instances++;

    fitnessValid = false;
    fitness = NAN;
    nn = new NeuralNetwork(numLayers, lSizes, true);
}

/* ------------------------------------------------------------------------ */

NNChromosome::NNChromosome(unsigned numLayers, const NNLayer** layers,
        const NNTrainData*td) :
    Chromosome(), trainData(td)
{
    instances++;

    fitnessValid = false;
    fitness = NAN;
    nn = new NeuralNetwork(numLayers, layers);
}

/* ------------------------------------------------------------------------ */

NNChromosome::NNChromosome(const NNChromosome& c) :
    Chromosome(), trainData(c.trainData)
{
    instances++;

    fitnessValid = false;
    fitness = NAN;
    nn = c.nn->copy();
}

/* ------------------------------------------------------------------------ */

NNChromosome::~NNChromosome()
{
    instances--;
    delete nn;
}

/* ------------------------------------------------------------------------ */

Chromosome* NNChromosome::copy() const
{
    return new NNChromosome(*this);
}

/* ------------------------------------------------------------------------ */

double NNChromosome::getFitness()
{
    if (fitnessValid) return fitness;

    const double* nnOutputs = NULL;
    const double* wantOutputs = NULL;
    unsigned numData = trainData->getNumUsedItems();
    double errSum = 0.0;
    double oneCharErr = 0.0;
    unsigned badBits = 0;

    /* For each train item, evaluate the sum of squared errors */
    for (unsigned i = 0; i < numData; i++)
    {
        nnOutputs = nn->eval(trainData->getNNInputs(i));
        wantOutputs = trainData->getNNOutputs(i);
        oneCharErr = 0.0;
        
        for (unsigned o = 0; o < nn->getNumOutputs(); o++)
        {
            double err = nnOutputs[o] - wantOutputs[o];
            oneCharErr += err*err;

            if (fabs(err) > errBitTreshold)
                badBits++;
        }
        //if (oneCharErr > 0.01)
            errSum += oneCharErr;
    }

    double maxSum = sqrt(trainData->getNumUsedItems()*nn->getNumOutputs());

    fitnessValid = true;
    fitness = (maxSum - sqrt(errSum)) / maxSum;

    //if (badBits == 0 && fitness >= errSumTreshold)
    //    fitness = 1.0;

    return fitness;
}

/* ------------------------------------------------------------------------ */

Chromosome *NNChromosome::crossover(const Chromosome *c) const
{
    const NNChromosome* chr = dynamic_cast<const NNChromosome *>(c);
    assert(chr != NULL);
    assert(checkTopologyMatch(this, chr));

    /* Create new NN layers */
    //const NNLayer** newLayers = combineLayers(chr);

    //Chromosome* offspring = new NNChromosome(nn->getNumLayers(), newLayers, this->trainData);
    //delete []newLayers;
    
    NNChromosome* offspring = static_cast<NNChromosome* >(this->copy());
    //offspring->nn->averageWeights(chr->nn);
    offspring->nn->combineUniform(chr->nn);

    return offspring;
}

/* ------------------------------------------------------------------------ */

const NNLayer** NNChromosome::combineLayers(const NNChromosome* chr) const
{
    const NNLayer** newLayers = 
        (const NNLayer** ) new NNLayer*[chr->nn->getNumLayers()];

    for (unsigned i = 0; i < chr->nn->getNumLayers(); i++)
    {
        double rnd = rand() / (RAND_MAX + 1.0);
        if (rnd < 0.5)
            newLayers[i] = this->nn->getLayer(i);
        else
            newLayers[i] = chr->nn->getLayer(i);
    }

    return newLayers;
}

/* ------------------------------------------------------------------------ */

void NNChromosome::mutate(double mRate)
{
    assert(mRate >= 0 && mRate <= 1.0);

    fitnessValid = false;
    nn->randomBias(mRate);
}

/* ------------------------------------------------------------------------ */

void NNChromosome::print(std::ostream &s, const NNTrainData* td) const
{
    s << endl;
    //nn->print(s);

    bool charOK;

    map<unsigned char, int> numInst;
    /* Number characters instances recognized OK */
    map<unsigned char, int> instOK;
    map<unsigned char, double> instSum;

    map<unsigned char, int>::iterator it;


    /* Print eval info */
    for (unsigned i = 0; i < td->getNumUsedItems(); i++)
    {
        unsigned char recChar = td->getChar(i);

        const double* outs = nn->eval(td->getNNInputs(i));
        const double* wantOuts = td->getNNOutputs(i);
        charOK = true;

        for (int o = NN_OUTPUT_WIDTH - 1; o >= 0; o--)
        {
            if (nearbyint(outs[o]) != wantOuts[o])
                charOK = false;
        }

        numInst[recChar]++;
        if (charOK)
            instOK[recChar]++;
    }

    s << "------------------------------------------------------" << endl;
    /* Print out all characters */
    for (it = numInst.begin(); it != numInst.end(); it++)
    {
        s << "Char: " << (*it).first  << ": " << instOK[(*it).first] <<
            " OK out of: " << (*it).second << endl;
    }

}

/* ------------------------------------------------------------------------ */

void NNChromosome::print(std::ostream &s) const
{
    this->print(s, trainData);
}

/* ------------------------------------------------------------------------ */

bool NNChromosome::
checkTopologyMatch(const NNChromosome* c1, const NNChromosome* c2) const
{
    if (c1->nn->getNumLayers() != c2->nn->getNumLayers())
        return false;

    for (unsigned l = 0; l < c1->nn->getNumLayers(); l++)
    {
        if (c1->nn->getLayer(l)->getNumItems() !=
            c2->nn->getLayer(l)->getNumItems())
        {
            return false;
        }
    }

    return true;
}

/* ------------------------------------------------------------------------ */

/* End of file NNChromosme.cc */
