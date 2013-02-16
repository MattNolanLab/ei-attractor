/**
 * Name:        testNN.cc
 * Project: 	ARC-2008-L Diffuse GA + Neural network
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Train neural network using GA.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <cstdlib>
#include <ctime>

#include "args.h"
#include "global.h"
#include "CellularGA.h"
#include "Chromosome.h"
#include "NNChromosome.h"
#include "ChromosomeFactory.h"
#include "NNTrainData.h"

using namespace std;


#ifdef ARC_MPI_ENABLED
#   include <mpi.h>
#endif


/* Function prototypes */
void printHelp();
void testNet(NeuralNetwork& nn, const NNTrainData& td, ostream& s,
        bool printStat, const char* testFile);


const char trainBitmapFile[] = "../fonts/trainData.vystrel0.05.gauss0.5.bmp";
const char trainTxtFile[]    = "../fonts/trainData.vystrel0.05.gauss0.5.txt";
//const char trainBitmapFile[] = "../fonts/trainData.allchars.bmp";
//const char trainTxtFile[]    = "../fonts/trainData.allchars.txt";
//const char trainBitmapFile[] = "../fonts/trainData.bmp";
//const char trainTxtFile[]    = "../fonts/trainData.txt";

const unsigned numTestFiles = 10;
const char* (testBitmapFiles[]) = {
    "../fonts/testData.vystrel.0.02.bmp",
    "../fonts/testData.vystrel.0.04.bmp",
    "../fonts/testData.vystrel.0.06.bmp",
    "../fonts/testData.vystrel.0.08.bmp",
    "../fonts/testData.vystrel.0.1.bmp",
    "../fonts/testData.vystrel.0.12.bmp",
    "../fonts/testData.vystrel.0.14.bmp",
    "../fonts/testData.vystrel.0.16.bmp",
    "../fonts/testData.vystrel.0.18.bmp",
    "../fonts/testData.vystrel.0.20.bmp"
};

const char* testTxtFiles[] = {
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt",
    "../fonts/testData.txt"
};


/* Chromosome specification */
const unsigned numLayers = 3;
const unsigned numHiddenNeurons1 = 95;
const unsigned numHiddenNeurons2 = 10;
unsigned lSizes[numLayers];


class NNGA : public CellularGA
{
  public:
    NNGA(InitParams init) : CellularGA(init) {}

    bool stopCondition() 
    {
        return globalBestIndividualFitness() >= 0.99;
    }
};

int main(int argc, char **argv)
{

    /* ------------------------------------------------------------------- */
    /* Start of GA code */
    /* ------------------------------------------------------------------- */
    for (unsigned i = 0; i < numTestFiles; i++)
    {
        ifstream test(testBitmapFiles[i]);
        if (!test.is_open())
        {
            cerr << "Could not open test file: " << testBitmapFiles[i] << endl;
            exit(EXIT_FAILURE);
        }
        test.close();

        test.open(testTxtFiles[i]);
        if (!test.is_open())
        {
            cerr << "Could not open test file: " << testTxtFiles[i] << endl;
            exit(EXIT_FAILURE);
        }
        test.close();
    }


    try 
    {
        srand(time(NULL));

        InitParams params = parseArgs(argc, argv);

        /* Initialize */
        NNTrainData trainData(trainBitmapFile, trainTxtFile);
        //trainData.setNumUsedItems(usedItems);

        lSizes[0] = trainData.getCharWidth() * trainData.getCharHeight();
        lSizes[1] = numHiddenNeurons1;
        //lSizes[2] = numHiddenNeurons2;
        lSizes[2] = NN_OUTPUT_WIDTH; /* defined in NNTrainData.h */

        NNChromosomeFactory factory(numLayers, lSizes, trainData);

        params.chrFactory = &factory;

        for (unsigned i = 0; i < trainData.getNumUsedItems(); i++)
            trainData.printNNInput(i, cout);

        NNGA* ga = new NNGA(params);
        ga->run();

        ostream& s = cout;
        s << "<gastat> maxGen: " << params.maxGen << " xsize: " << params.gridXSize <<
              " ysize: " << params.gridYSize <<
              " gennum: " << ga->getGenNum() << " bestFitness: " <<
            ga->globalBestIndividualFitness() << endl;

        const Chromosome* c = ga->getGlobalBestIndividual();
        NNChromosome* nnc = static_cast<NNChromosome* >(c->copy());

        for (unsigned i = 0; i < numTestFiles; i++)
        {
            /* Evaluate best individual for generalization */
            if (testBitmapFiles != NULL)
            {
                NNTrainData testData(testBitmapFiles[i], testTxtFiles[i]);
                testNet(*nnc->getNN(), testData, s, true, testBitmapFiles[i]);
            }
        }


        delete ga;
        delete nnc;
    }
    catch (GenericException &e)
    {
        ostringstream s;
        s << e.what() << endl;
        printErr(true, s);
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        finalize();
        exit(EXIT_FAILURE);
    }


    /* ------------------------------------------------------------------- */
    /* End of GA code */
    /* ------------------------------------------------------------------- */
}

void testNet(NeuralNetwork& nn, const NNTrainData& td, ostream& s,
        bool printStat, const char* testFile)
{
    map<unsigned char, int> numInst;
    /* Number characters instances recognized OK */
    map<unsigned char, int> instOK;
    map<unsigned char, double> instSum;

    map<unsigned char, int>::iterator it;


    /* Print eval info */
    for (unsigned i = 0; i < td.getNumUsedItems(); i++)
    {
        unsigned char recChar = td.getChar(i);

        const double* outs = nn.eval(td.getNNInputs(i));
        const double* wantOuts = td.getNNOutputs(i);
        bool charOK = true;

        for (int o = NN_OUTPUT_WIDTH - 1; o >= 0; o--)
        {
            if (nearbyint(outs[o]) != wantOuts[o])
                charOK = false;
        }

        numInst[recChar]++;
        instOK[recChar] += (int) charOK;
    }

    unsigned numItems = 0;
    unsigned numItemsOK = 0;

    s << "------------------------------------------------------" << endl;
    /* Print out all characters */
    for (it = numInst.begin(); it != numInst.end(); it++)
    {
        s << "Char: " << (*it).first  << " " << instOK[(*it).first] <<
            " OK out of: " << (*it).second << endl;

        numItems += (*it).second;
        numItemsOK += instOK[(*it).first];
    }

    if (printStat)
    {
        /* Print overall statistics */
        s << "<" << testFile << "><generalizationstat> numChars: " << numItems << " numCharsOK: " <<
            numItemsOK << endl;
    }
}



/** Print help end exit **/
void printHelp()
{
    cout <<
        "testNN [-h] [-g gen] [-c cRate] [-m mRate] [-xsize xs] [-ysize ys] "
        " [-r rLength]" << endl
        << endl <<
        "Parameters description" << endl << endl <<

        "  -g     Maximal number of generations. If unspecified or zero, then" << endl <<
        "         unlimited." << endl <<
        "  -ct    Crossover type: " << C_RANDOMWALK << " - randomWalk, " << C_AREA << " - best individual in nearby area." << endl <<
        "  -c     Crossover rate. Default value is 0.7." << endl <<
        "  -m     Mutation rate. Default value is 0.0001." << endl <<
        "  -xsize Number of columns in the population grid. Defaul value is 10." << endl <<
        "  -ysize Number of rows in the populatin grid. Default value is 10." << endl <<
        "  -r     Random walk length/area radius - used in crossover selection" <<
        "  -h     Print this help." << endl;
}


/* End of file testNN.cc */
