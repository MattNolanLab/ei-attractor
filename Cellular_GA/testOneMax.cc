/**
 * Name:        testGA.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Test simple OneMax GA implementation.
 */

#include <sstream>

#include <cstdlib>
#include <ctime>

#include "args.h"
#include "global.h"
#include "CellularGA.h"
#include "Chromosome.h"
#include "OneMaxChromosome.h"
#include "ChromosomeFactory.h"


#ifdef ARC_MPI_ENABLED
#   include <mpi.h>
#endif


/* Function prototypes */
void printHelp();


/* Chromosome specification */
#ifndef CHR_LENGTH
#   define CHR_LENGTH 10000
#endif
const int chrLength = CHR_LENGTH;
OneMaxChromosomeFactory factory(chrLength);


/**
 * OneMax GA
 */
class OneMaxGA : public CellularGA
{
  public:

    OneMaxGA(InitParams init) : CellularGA(init) {};

    bool stopCondition()
    {
        return globalBestIndividualFitness() >= chrLength;
    };
};

int main(int argc, char **argv)
{

#ifdef ARC_MPI_ENABLED
    double start;
    double stop;
    int myRank;
    int procNum;

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        cerr << "MPI_Init() failed!" << endl;
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    start = MPI_Wtime();
#endif

    /* ------------------------------------------------------------------- */
    /* Start of GA code */
    /* ------------------------------------------------------------------- */
    int genNum = 0;
    OneMaxGA* ga = NULL;
    InitParams init = parseArgs(argc, argv);

    try 
    {
        unsigned int i;

        init.chrFactory = &factory;

        srand(time(NULL) + i);
        ga = new OneMaxGA(init);
        ga->run();

        genNum = ga->getGenNum();
    }
    catch (GAException &e)
    {
        delete ga;
        ostringstream s;
        s << e.what() << endl;
        printErr(true, s);
    }
    catch (exception &e)
    {
        delete ga;

        cerr << e.what() << endl;
        finalize();
        exit(EXIT_FAILURE);
    }




    /* ------------------------------------------------------------------- */
    /* End of GA code */
    /* ------------------------------------------------------------------- */

#ifdef ARC_MPI_ENABLED
    stop = MPI_Wtime();

    //cerr << "Average chromosome transmission time: " << ga->getAvgChrTime()
    //        << endl;

    //OneMaxChromosome::printNumberOfInstances(cerr);
    
    if (myRank == ROOT_NODE)
    {
        cerr << "<gastat> procnum: " << procNum << " gennum: " << genNum <<
            " time: " <<
        cerr.width(5);
        cerr << right << stop - start << endl;
    }

    MPI_Finalize();
#else
    cerr << "<gastat> maxGen: " << init.maxGen << " xsize: " << init.gridXSize <<
          " ysize: " << init.gridYSize <<
          " chrLength: " << chrLength << " gennum: " << genNum << " bestFitness: " <<
        ga->globalBestIndividualFitness() << endl;
#endif

    delete ga;

}




/** Print help end exit **/
void printHelp()
{
    cout <<
        "testGA [-h] [-g gen] [-c cRate] [-m mRate] [-xsize xs] [-ysize ys] "
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
        "         each individual can communicate only with its neighbours." << endl << 
        "  -h     Print this help." << endl;
}


/* End of file testGA.cc */
