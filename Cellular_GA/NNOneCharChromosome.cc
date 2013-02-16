/**
 * Name:        NNOneCharChromosome.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural Network chromosome implementation - one character
 *              recognition.
 */

#include <map>

#include "NNOneCharChromosome.h"
#include "debug.h"
#include "global.h"

using namespace std;

/* ------------------------------------------------------------------------ */

NNOneCharChromosome::
NNOneCharChromosome(unsigned inputSize, unsigned NNNodes, const NNTrainData* td,
            unsigned char symbol) :
    NNChromosome(3, getLSizes(inputSize, NNNodes), td), symbol(symbol)
{
}

/* ------------------------------------------------------------------------ */

NNOneCharChromosome::NNOneCharChromosome(const NNOneCharChromosome& c) :
    NNChromosome(c), symbol(c.symbol)
{
}

/* ------------------------------------------------------------------------ */

NNOneCharChromosome::NNOneCharChromosome(const NNLayer** l, const NNTrainData* td,
        unsigned char s) :
    NNChromosome(3, l, td), symbol(s)
{
}

/* ------------------------------------------------------------------------ */

Chromosome* NNOneCharChromosome::copy() const
{
    return new NNOneCharChromosome(*this);
}

/* ------------------------------------------------------------------------ */

double NNOneCharChromosome::getFitness()
{
    if (fitnessValid) return fitness;

    const double* nnOutputs = NULL;
    double wantOutput;
    unsigned numData = trainData->getNumUsedItems();
    unsigned badBits = 0;
    double err = 0.0;
    double errSum = 0.0;

    lastBadBit = trainData->getNumUsedItems();

    /*
     * Evaluate each char. If the char is our symbol, NN should return 1.
     * Otherwise it should return 0
     */
    for (unsigned i = 0; i < numData; i++)
    {
        if ((unsigned) trainData->getChar(i) == symbol)
            wantOutput = 1.0;
        else
            wantOutput = 0;

        nnOutputs = nn->eval(trainData->getNNInputs(i));
        err = nnOutputs[0] - wantOutput;

        //if (trainData->getChar(i) != symbol)
        //    if (fabs(err) < 0.01) err = 0;

        errSum += err*err;

        if (fabs(nnOutputs[0] - wantOutput) > errBitTreshold)
        {
            badBits++;
            lastBadBit = i;
        }
    }

    double maxSum = sqrt(trainData->getNumUsedItems()); /* One neuron output !! */

    fitnessValid = true;
    fitness = (maxSum - sqrt(errSum)) / maxSum;

    if (badBits == 0)
    {
        if(fitness >= errSumTreshold)
            fitness = 1.0;
    }
    //else
    //{
    //    /* Penalize */
    //    fitness -= 2.5 * badBits/trainData->getNumUsedItems();
    //}

    return fitness;
}

/* ------------------------------------------------------------------------ */

Chromosome *NNOneCharChromosome::crossover(const Chromosome *c) const
{
    const NNOneCharChromosome* chr = dynamic_cast<const NNOneCharChromosome *>(c);
    assert(chr != NULL);
    assert(checkTopologyMatch(this, chr));

    /* Create new NN layers */
    const NNLayer** newLayers = combineLayers(chr);

    Chromosome* offspring = new NNOneCharChromosome(newLayers, this->trainData, symbol);
    delete []newLayers;
    
    //NNChromosome* offspring = static_cast<NNChromosome* >(this->copy());
    //offspring->nn->combineUniform(chr->nn);

    return offspring;
}

/* ------------------------------------------------------------------------ */

void NNOneCharChromosome::print(std::ostream &s, const NNTrainData* td) const
{
    s << endl;
    nn->print(s);


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

        numInst[recChar]++;

        if (recChar == symbol)
        {
            if ((int)nearbyint(outs[0]) == 1)
                instOK[recChar]++;
        }
        else
        {
            if ((int)nearbyint(outs[0]) == 0)
                instOK[recChar]++;
        }

        instSum[recChar] += outs[0];
    }

    s << "------------------------------------------------------" << endl;
    /* Print out all characters */
    for (it = numInst.begin(); it != numInst.end(); it++)
    {
        s << "Char: " << (*it).first  << " " << instOK[(*it).first] <<
            " OK out of: " << (*it).second <<
            ", average: " << instSum[(*it).first]/ (*it).second << endl;
    }

    s << "Trained Char: " << symbol << endl;
    s << "------------------------------------------------------" << endl;
}

/* ------------------------------------------------------------------------ */

void NNOneCharChromosome::print(std::ostream &s) const
{
    this->print(s, trainData);
}

/* ------------------------------------------------------------------------ */

/* End of file NNOneCharChromosme.cc */
