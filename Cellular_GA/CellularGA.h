/**
 * Name:        CellularGA.h
 * Project: 	ARC-2008-L Cellular GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Cellular genetic algorithm header file
 */

#ifndef CELLULARGA_H
#define CELLULARGA_H

#include "global.h"
#include "Chromosome.h"
#include "ChromosomeFactory.h"
#include "error.h"


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
#   include <mpi.h>
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */


/**
 * GA statistics structure.
 */
struct GAStat
{
    /** Best Chromosome **/
    Chromosome *bestChromosome;


    /** Best chromosome fitness **/
    double bestFitness;


    /** Population Fitness sum **/
    double fitnessSum;

    GAStat()
    {
        bestChromosome = NULL;
    }

    ~GAStat()
    {
        if (bestChromosome != NULL && !bestChromosome->isLocal())
            delete bestChromosome;
    }

    /** Replace the best with a copy of c **/
    void replaceBest(const Chromosome *c)
    {
        assert(!(bestChromosome != NULL && bestChromosome->isLocal()));

        delete bestChromosome;
        bestChromosome = c->copy();
        bestChromosome->setLocal(false);
    }
};

/**
 * CellularGA init params.
 */
struct InitParams
{
    /** Generation number limit **/
    int maxGen; 
    /** Number of created offspring ratio in one generation. **/
    float crossoverRate;  
    /** Generic mutation ratio **/
    float mutationRate;
    /** Population grid number of columns **/
    int gridXSize;
    /** Population grid number of rows **/
    int gridYSize;
    /** Random walk length **/
    unsigned int randomWalkLength;
    /** Crossover type. See types declarations **/
    unsigned int crossoverType;
    /** Chromosome factory **/
    ChromosomeFactory *chrFactory;
};


/**
 * 2D grid/grid size.
 * Represents two values packed in a structure.
 */
struct Grid
{
    /** X value/size **/
    int x;
    /** Y Value/size **/
    int y;

    Grid() : x(0), y(0) {};
    Grid(unsigned int x, unsigned int y) : x(x), y(y) {};
};


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
/**
 * Global MPI variables
 */
class GlobalMPI
{
    /** Size of the (square) matrix in MPI, ie. Process matrix **/
    int MatrixSize;

    /**
     * Sizes of the matrix each processor has. Processors in the last row or
     * columnt may have smaller matrix than this size.
     */
    Grid ProcGridSize;


  public:
    GlobalMPI(const InitParams& p);

    /**
     * Get size of the square matrix of processes.
     */
    int getMatrixSize() const { return MatrixSize; }

    /**
     * Get processor grid size.
     */
    Grid getProcGridSize() const { return ProcGridSize; }

    /**
     * Get my position in the MPIMatrixGrid.
     *
     * @return My position.
     */
    Grid getMyMatrixPos() const
    {
        int myRank = getMyRank();
        return
            Grid(myRank % MatrixSize, myRank / MatrixSize);
    }


    /**
     * Get my process rank
     */
    int getMyRank() const
    {
        int err;
        int myRank = -1;
        err = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        CHECK_MPI_ERROR(err);
    
        return myRank;
    }


    /** 
     * Get number of processes.
     */
    int getProcNum() const
    {
        int procNum = -1;
        int err;

        err = MPI_Comm_size(MPI_COMM_WORLD, &procNum);
        CHECK_MPI_ERROR(err);

        return procNum;
    }


    /**
     * Get process rank for the given Chromosome global position.
     */
    int translateGlobalPosToRank(Grid p) const
    {
        int procX, procY;

        if (p.x >= MatrixSize * ProcGridSize.x)
            procX = MatrixSize - 1;
        else
            procX = p.x / (ProcGridSize.x);

        if (p.y >= MatrixSize * ProcGridSize.y)
            procY = MatrixSize - 1;
        else
            procY = p.y / (ProcGridSize.y);

        return procY * MatrixSize + procX;
    }


    /**
     * Translate Global position to my local position.
     */
    Grid translateGlobalPosToLocal(Grid gp) const
    {
        assert(translateGlobalPosToRank(gp) == getMyRank());

        Grid mPos;  /* My matrix position */
        mPos = getMyMatrixPos();

        return Grid(gp.x - mPos.x * ProcGridSize.x,
                    gp.y - mPos.y * ProcGridSize.y);
    }

    /**
     * Get the processor rank at the given position (at MPI processors grid).
     * This function has only
     * meaning when ARC_MPI_ENABLED is defined.
     *
     * @param p Specified global postition.
     * @return Process rank
     */
    int getProcessRank(Grid p) const
    {
        /* For now, only grid(block matrix) mapping is supported */
        int procX = p.x / ProcGridSize.x;
        int procY = p.y / ProcGridSize.y;
    
        return procY * MatrixSize + procX;
    }
};
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */


/**
 * Cellular genetic algorithm.
 */
class CellularGA
{
  public:

    /**
     * Parse command line arguments and return initialization parameters.
     */
    static InitParams parseArgs(int argc, char **argv);
    
    /**
     * Constructor
     *
     * @param mutationRate 
     * @paratm popSize Population size
     */
    CellularGA(const InitParams initParams);
    ~CellularGA();

    /**
     * Start the genetic algorithm
     */
    void run();

    /**
     * Step one generation further.
     */
    void step();

    /** Get mutation rate **/
    float getMutRate() const;

    /** Get crossover rate **/
    float getCrossRate() const;

    /**
     * Get average chromosome transmision time, in seconds.
     */
    double getAvgChrTime() { return chrTime / chrReq; }

    /** Get number of generations elapsed so far **/
    int getGenNum() { return generation; }

    /**
     * Get the best chromosome on the global grid. If MPI is used, this function
     * may only be called from the node 0.
     */
    const Chromosome *getGlobalBestIndividual();

    /**
     * Get the best chromosome's fitness.
     */
    double globalBestIndividualFitness();

    /**
     * Return local average population fitness.
     */
    double globalAvgFitness();


  protected:

    /**
     * Is the stop condition met?
     * This function must be defined in the specific subclass.
     *
     * @return True if the stop condition is met, false otherwise.
     */
    virtual bool stopCondition() = 0;

  private:

    /** Stop flag/command from root node */
    bool stopFlag;

    /** Local population **/
    Chromosome ***pop;

    /** Current generation number **/
    int generation;

    /** Individuals swapped in the local generation **/
    int individualsSwapped;

    /** X mapping for random walk **/
    static const int randomWalkXMap[];
    static const int randomWalkYMap[];

    /** Best individual reference **/
    Chromosome *bestIndividual;

    /** 
     * Summed fitness of the whole local population
     */
    double fitnessSum;


    /** Global statistics **/
    GAStat globalStat;


    /** Debug Statistics **/
    double chrTime;
    int chrReq;



    /** Get local grid parameters **/
    Grid getLocalGrid() const;


    /**
     * Recalculate all possible statistics over the local population.
     */
    void recalculatePopulation();


    /**
     * Do a random walk initiating from the point p. Pick up chromosome with the
     * best fitness along the walk and return it.
     *
     * @return Chromosome with the best fitness along the walk.
     */
    Chromosome *randomWalk(Grid p);

    /**
     * Get best individual in the square neighbourhood.
     *
     * @param p Center of the area, local grid.
     * @param r Area radius.
     */
    Chromosome* bestInArea(Grid p, unsigned r);


    /**
     * Get new random position with step 1 (manhattan metric).
     */
    Grid getRandomWalkNewPosition(Grid p);


    /**
     * Return a chromosome, determined by the _global_ position p.
     *
     * @param p Position of the chromosome in the GLOBAL grid. This is used in
     *   when the ARC_MPI_ENABLED flag is defined. This function may communicate
     *   an individual from different demes, using MPI. In case MPI is disabled,
     *   the global grid == local grid.
     * @return Chromosome at the given position. If MPI is enabled, this
     * chromosome will have isLocal flag set to false and must be properly
     * handled and deleted by the application.
     */
    Chromosome *getGlobalChromosome(Grid p);


    /** 
     * Return the _local_ best individual
     * @return The best individual.
     */
    Chromosome *getBestIndividual() const { return bestIndividual; }


    /**
     * Return local average population fitness.
     */
    double avgFitness()
    {
        return fitnessSum / (lGrid.x * lGrid.y);
    }

    

    /**
     * Get Best individual's fitness.
     *
     * @return Its fitness
     */
    double bestIndividualFitness() {
        return bestIndividual->getFitness();
    }


    /**
     * Initialize the population with random chromosomes.
     */
    void initPopulation();


    /** Get random value in the range <0, max) **/
    unsigned long getRandomInt(unsigned long max) const
    {
        return (int)(max * (rand() / (RAND_MAX + 1.0)));
    }

    /** Get a random point in the local grid **/
    Grid getLocalRandomPoint() const
    {
        return Grid (getRandomInt(lGrid.x), getRandomInt(lGrid.y));
    }


    /**
     * Translate position at local grid to global grid coordinates.
     *
     * @param p Local grid position.
     * @return Global grid position.
     */
    Grid translateLocalToGlobalPos(Grid p) const;


    /**
     * Delete and swap an individual in the local grid at specified position,
     * with specified new individual.
     *
     * @param p Position in the local grid
     * @param c New Chromosome that will replace an old one.
     */
    void swapIndividual(Grid p, Chromosome *c);



    /* -------------------------------------------------------------------
     * MPI communication functions and members
     * ------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
#   define MPINEIGHBOUR_MAX 4

    static const int GAMPI_TAG_CHRSEND;
    static const int GAMPI_TAG_CHRREQ;
    static const int GAMPI_TAG_FSUM;
    static const int GAMPI_TAG_WHATNEXT;
    static const int GAMPI_TAG_BESTCHR;

    /* GAMPI_WHATNEXT tag message context */
#   define STOP             0 /* Finish */
#   define CONTINUE         1 /* Continue into the next generation */
#   define GET_CHROMOSOME   2 /* Give me your best chromosome */

    int remoteChromosomesRequested;
    int chromosomesRequested;

    /** MPI parameters **/
    GlobalMPI gMPI;
    
    /**
     * MPI neighbours list
     */
    struct MPINeighbours
    {
        /**
         * An array of MPI requests for Chromosomes from the neighbours.
         */
        MPI_Request requests[MPINEIGHBOUR_MAX];

        /**
         * Array of ranks for the MPI Chromosome requests.
         */
        int ranks[MPINEIGHBOUR_MAX];

        /**
         * Buffers, two int values, X and Y global grid positions.
         */
        int bufs[MPINEIGHBOUR_MAX][2];

        MPINeighbours(const GlobalMPI &g);


        /**
         * Check if there are any requests and if so, indicate that in the
         * returning parameters. If the function returns true, the user _MUST_
         * send a Chromosome to the specified processor. If the chromosome is
         * not sent, deadlock will happen.
         *
         * @param p Global requested chromosome position.
         * @param r Destination rank.
         * @return True if a chromosome is requested.
         */
        bool checkRequests(Grid &p, int &r);


        /**
         * Check if there is any request present at position given by index.
         *
         * @param p Global requested chromosome position.
         * @param r Destination rank.
         * @param index Neighbours index.
         * @return True if chromosome is requested.
         */
        bool checkRequests(Grid& p, int& r, int index);


        /**
         * Get the size of the neighbours array.
         */
        int getSize() { return size; }

      private:
        const GlobalMPI &gMPI;

        /**
         * List size
         */
        int size;

    };


    /** 
     * MPI statistics collector.
     */
    struct MPIStats;
    friend struct MPIStats;
    struct MPIStats
    {
      public:
        /**
         * Constructor.
         *
         * @param g  Global MPI functions.
         * @param ch Temporary chromosome.
         */
        MPIStats(CellularGA &ga, const GlobalMPI &g, Chromosome *ch);

        ~MPIStats();


        /** 
         * Gets statistics.
         * 
         * @param stat This process's local GA statistics.
         * @return True if all statistics are collected. Then the stat parameter
         *         will hold proper statistics (and rewrite the values in the
         *         stat parameter).
         */
        bool getStat(GAStat &stat);


        /**
         * Begin new non-blocking statistics receives.
         */
        void beginStat();

      private:

        /**
         * [0] --> Fitness sum
         * [1] --> Fitness of the best chromosome
         */
        double (*fitness)[2];

        /** MPI Fitness requests array **/
        MPI_Request *fitnessReqs;

        /** CellularGA handler **/
        CellularGA &cga;

        /** Global MPI handler reference **/
        const GlobalMPI &gMPI;

        /** Temporary chromosome **/
        Chromosome *tempChr;


        /** Number of received statistics **/
        int numRecv;
    };


    /**
     * My MPI neighbours
     */
    MPINeighbours neighbours;


    /**
     * MPI Statistics management.
     * Only ROOT_NODE will allocate this.
     */
    MPIStats *mpiStats;


    /**
     * Print statistics on chromosome remote communication.
     */
    void printRemoteChromosomeRatio()
    {
        cout << "Ratio of remote chromosomes: " <<
            ((float) remoteChromosomesRequested) / chromosomesRequested << endl;

        chromosomesRequested = 0;
        remoteChromosomesRequested = 0;
    }


    /**
     * Request a chromosome at the local position p in processor at specified
     * rank
     *
     * @param p global position
     * @param rank Process rank
     */
    Chromosome *getMPIChromosome(Grid p, int rank);


    /**
     * Send statistics to the ROOT processor (usually rank 0). Blocks the user,
     * but doesn't block MPI communication to avoid deadlock.
     *
     * @return True is the GA should continue.
     */
    bool sendMPIStatistics();


    /**
     * Begin statistics collection. You should use this function only at
     * processor with ROOT_NODE rank number.
     *
     */
    void recvMPIStatistics();


    /**
     * Send chromosome to processor p. Blocks the user and block MPI
     * communication --> user must run non blocking recv at the same time at
     * processor p.
     *
     * @param p    Chromosome position in the local grid
     * @param rank Rank of the receiving process
     */
    void sendMPIChromosome(Grid p, int rank, int tag = GAMPI_TAG_CHRSEND);


    /**
     * Send chromosome to processor p. Blocks the user and block MPI
     * communication --> user must run non blocking recv at the same time at
     * processor p.
     *
     * @param c    Chromosome pointer.
     * @param rank Rank of the receiving process
     */
    void
    sendMPIChromosome(const Chromosome *c, int rank, int tag = GAMPI_TAG_CHRSEND);

    /**
     * Poll MPI action. Use this function to check if an MPI action needs be
     * performed. Doesn't block.
     */
    void pollMPIAction();


    /**
     * Set needed MPI parameters.
     */
    void setMPIParams();


#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */


    /**
     * GA initialization parameters.
     * This is valid for whole grid. The grid may be divided into several
     * subgrids, in order to partition the problem into several processors.
     * Therefore, not the popX/YSize parameters are used, but the local ones.
     */
    const InitParams params;


    /** Local grid parameters **/
    const Grid lGrid;

};


#endif /* CELLULARGA_H */

/* End of file CellularGA.h */
