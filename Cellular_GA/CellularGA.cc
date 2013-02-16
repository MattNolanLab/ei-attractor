/**
 * Name:        CellularGA.cc
 * Project: 	ARC-2008-L Cellular GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Cellular genetic algorithm.
 */

#include <iostream>
#include <exception>
#include <sstream>
#include <cmath>

#include "CellularGA.h"
#include "debug.h"
#include "error.h"

using namespace std;



/**
 * This is a cellular GA implementation which uses an MPI interface. The main
 * idea is to divide the whole grid into several processes (processors).
 * Each process takes a part of the whole grid of individuals and runs the
 * diffuse genetic algorithm (DGA) on its own (the operation of the DGA is
 * described later). If any of the processes needs an individual from the other
 * process, the individual is communicated via the MPI interace functions.
 *
 * This communication is transparent to the whole genetic algorithm and is
 * compile-time selectable. If the user wants to use just the pseudo-parallel
 * variant, it is possible and no MPI usage/linking will be done.
 *
 * The operation of a Diffuse GA is as follows:
 *
 * 1) Create the grid (population). Create only local grid, based on the
 *   processors number, grid parameters and processors mapping
 *   (block/cyclic-block).
 * 2) Choose a random point P in the grid. Make a random walk from P of a
 *   defined length (less than the length of the grid) and pick up the best
 *   individual. Repeat the walk to pick up another individual.
 * 3) Using the two selected individuals, make crossover on them to produce
 *   an offspring.
 * 4) Mutate the offspring
 * 5) If the offspring's fitness is better than the fitness of the individual
 *   at position P, replace that individual at P.
 * 6) Repeat steps 2) - 5) until the specified number of offsprings has been
 *   created. At this point, this is considerd a new generation
 * 7) Repeat steps 2) -6) until the stop condition is met.
 */

/* ------------------------------------------------------------------------ */

const int CellularGA::randomWalkXMap[] = { 0, 1, 0, -1 };
const int CellularGA::randomWalkYMap[] = {-1, 0, 1,  0 };

/* ------------------------------------------------------------------------ */

CellularGA::CellularGA(const InitParams initParams) :
/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    gMPI(initParams), 
    neighbours(gMPI),
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    params(initParams), lGrid(getLocalGrid())
{

    bestIndividual = NULL;
    fitnessSum = 0.0;
    generation = 1;
    stopFlag = false;
    individualsSwapped = 0;

    globalStat.bestChromosome = NULL;
    globalStat.fitnessSum = -1.0;


/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    remoteChromosomesRequested = 0;
    chromosomesRequested = 0;

    if (gMPI.getMyRank() == ROOT_NODE)
    {
        mpiStats = new MPIStats(*this, gMPI, params.chrFactory->getChromosome());
    }
    else
    {
        mpiStats = NULL;
    }

    if (params.gridXSize < gMPI.getMatrixSize() ||
        params.gridYSize < gMPI.getMatrixSize())
    {
        ostringstream msg;
        msg << "Too small grid size for this number of processors! " <<
               "Increase the grid size to be at least (" <<
               gMPI.getMatrixSize() << ", " << gMPI.getMatrixSize() <<
               ")";
        throw GAException(msg);
    }

#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    initPopulation();

    if (params.gridXSize == 1 && params.gridYSize == 1)
    {
        this->~CellularGA();
        ostringstream msg;
        msg << "Grid of size (1, 1) is useless!";
        throw GAException(msg);
    }

    dbg_msg("lGrid(x,y): (" << lGrid.x << ", " << lGrid.y << ")" << endl);

    //MPI_Finalize();
    //exit(1);
}

/* ------------------------------------------------------------------------ */

CellularGA::~CellularGA()
{

    /* Delete population */
    for (int x = 0; x < lGrid.x; x++)
        for (int y = 0; y < lGrid.y; y++)
            delete pop[x][y];

    for (int i = 0; i < lGrid.x; i++)
        delete []pop[i];
    delete []pop;

/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    delete mpiStats;
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

}

/* ------------------------------------------------------------------------ */

void CellularGA::initPopulation()
{
    pop = new Chromosome**[lGrid.x];
    for (int i = 0; i < lGrid.x; i++)
    {
        pop[i] = new Chromosome*[lGrid.y];
    }

    /* Create grid, composed of random individuals */
    for (int x = 0; x < lGrid.x; x++)
    {
        for (int y = 0; y < lGrid.y; y++)
        {
            pop[x][y] = params.chrFactory->getChromosome();
        }
    }

    recalculatePopulation();
}

/* ------------------------------------------------------------------------ */

Grid CellularGA::getLocalGrid() const
{
#ifdef ARC_MPI_ENABLED
    Grid ret;

    /*
     * Check if the global grid size is divisible by the number of processors.
     * If so, the local grid size == MPIProcGridSize. Otherwise the processors
     * in the last row/column will have shorter grids.
     */
    if (params.gridXSize % gMPI.getMatrixSize() == 0)
    {
        ret.x = params.gridXSize / gMPI.getMatrixSize();
    }
    else
    {
        int localGridX = floor((float)params.gridXSize / gMPI.getMatrixSize());
        ret.x = localGridX;

        //dbg_msg("My matrix pos X: " << gMPI.getMyMatrixPos().x << endl);
        
        if (gMPI.getMyMatrixPos().x == gMPI.getMatrixSize() - 1)
            ret.x = params.gridXSize - (gMPI.getMatrixSize() - 1) * localGridX;

    }

    if (params.gridYSize % gMPI.getMatrixSize() == 0)
    {
        ret.y = params.gridYSize / gMPI.getMatrixSize();
    }
    else
    {
        int localGridY = floor((float)params.gridYSize / gMPI.getMatrixSize());
        ret.y = localGridY;

        //dbg_msg("My matrix pos Y: " << gMPI.getMyMatrixPos().y << endl);

        if(gMPI.getMyMatrixPos().y == gMPI.getMatrixSize() - 1)
            ret.y = params.gridYSize - (gMPI.getMatrixSize() - 1) * localGridY;
    }

    return ret;
#else
    return Grid(params.gridXSize, params.gridYSize);
#endif

}


/* ------------------------------------------------------------------------ */

void CellularGA::run()
{
    //cout << "Initial best individual: fitness: " << globalBestIndividualFitness()
    //     << ", value: "; bestIndividual->print(cout);

    while (!stopCondition() &&
           !stopFlag &&
           (generation <= params.maxGen || params.maxGen == 0))
    {
        step();
    }


/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    if (gMPI.getMyRank() == ROOT_NODE)
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    {
        cout << "Best individual: ";
        bestIndividual->print(cout);
        cout << "fitness: " << globalBestIndividualFitness() << endl;
    }
}

/* ------------------------------------------------------------------------ */

void CellularGA::step()
{
    /* Determine the number of offsprings produced in one generation */
    float offNum = params.crossoverRate * lGrid.x * lGrid.y;
    int genOffsprings = (int) offNum;
    static double lastBestFitness = bestIndividualFitness();

    /* If the number of chromosome to reproduce is < 1.0, take it randomly */
    if (offNum < 1.0)
    {
        float rnd = rand() / (RAND_MAX + 1.0);
        if (rnd <= params.crossoverRate)
            genOffsprings = 1;
        else
            genOffsprings = 0;
    }

    for (int cnt = 0; cnt < genOffsprings; cnt++)
    {
        /* Determine a random point in the local grid */
        Grid p = getLocalRandomPoint();

#ifdef ARC_MPI_ENABLED
        pollMPIAction();
#endif
    
        /* Pick up chromosome for crossover and mutation */
        const Chromosome* c1;
        if (params.crossoverType == C_RANDOMWALK)
            c1 = randomWalk(p);
        else if (params.crossoverType == C_AREA)
            c1 = bestInArea(p, params.randomWalkLength);
        else
            assert(false);

#ifdef ARC_MPI_ENABLED
        pollMPIAction();
#endif
    
        Chromosome *offspring = pop[p.x][p.y]->crossover(c1);

        /**
         * If any of the c1 or c2 are non-local delete them.
         */
        if (!c1->isLocal()) delete c1;

        offspring->mutate(params.mutationRate);
    
        //if (offspring->getFitness() > pop[p.x][p.y]->getFitness())
        if (pop[p.x][p.y] != bestIndividual ||
            offspring->getFitness() > pop[p.x][p.y]->getFitness())
        {
            swapIndividual(p, offspring);
            individualsSwapped++;
        }
        else
        {
            delete offspring;
        }

    }

    /* Pick some inidividuals and mutate them */
    //double nMut = params.mutationRate * lGrid.x * lGrid.y;
    //for (unsigned i = 0; i < nMut; i++)
    //{
    //    Grid p = getLocalRandomPoint();
    //    if (pop[p.x][p.y] != bestIndividual)
    //        pop[p.x][p.y]->mutate(params.mutationRate);
    //}
    //recalculatePopulation();



/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    if (gMPI.getMyRank() == ROOT_NODE)
    {
        recvMPIStatistics();
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */


        if (generation % 1 == 0)//bestIndividualFitness() != lastBestFitness)
        {
            lastBestFitness = bestIndividualFitness();

            cout << "Best individual: ";
            getGlobalBestIndividual()->print(cout);
            cout << endl;

            cout << "G[";
            cout.width(5);
            cout << right << generation << "]: Avg. fitness: ";
            cout.setf(ios::fixed); cout.precision(5);
            cout << avgFitness() << ", Best fitness: ";
            cout.width(5);
            cout << globalBestIndividualFitness() << endl;

            //cerr << "Individuals swapped: " << individualsSwapped << endl;
        }

/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    }
    else
    {
        if (!sendMPIStatistics())
            stopFlag = true;    /* Finish command (got from ROOT_NODE) */
    }

    dbg_func(printRemoteChromosomeRatio());
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    generation++;
}

/* ------------------------------------------------------------------------ */

Chromosome *CellularGA::randomWalk(Grid p)
{
    /* Get the chromosome with the highest, originating from position p */
    Grid gPos = translateLocalToGlobalPos(p); 
    Chromosome *maxChr = getGlobalChromosome(gPos);
    Chromosome *tempChr = NULL;

    for (unsigned int i = 0; i < params.randomWalkLength; i++)
    {
        gPos = getRandomWalkNewPosition(gPos);
        tempChr = getGlobalChromosome(gPos);
        if (tempChr->getFitness() > maxChr->getFitness())
        {
            if (!maxChr->isLocal())
                delete maxChr;
            maxChr = tempChr;
        }
        else
        {
            if (!tempChr->isLocal())
                delete tempChr;
        }
    }

    return maxChr;
}

/* ------------------------------------------------------------------------ */

Chromosome* CellularGA::bestInArea(Grid p, unsigned r)
{
    Grid gPos = translateLocalToGlobalPos(p);

    int startX = gPos.x - r;
    int startY = gPos.y - r;

    if (startX < 0) startX = 0;
    if (startY < 0) startY = 0;

    int endX = gPos.x + r;
    int endY = gPos.y + r;

    if (endX >= params.gridXSize) endX = params.gridXSize - 1;
    if (endY >= params.gridYSize) endY = params.gridYSize - 1;

    Chromosome* maxChr = getGlobalChromosome(Grid(startX, startY));
    Grid pos;

    for (pos.x = startX; pos.x <= endX; pos.x++)
    {
        for (pos.y = startY; pos.y <= endY; pos.y++)
        {
            if (pos.x != gPos.x && pos.y != gPos.y)
            {
                Chromosome* c = getGlobalChromosome(pos);
                if (c->getFitness() > maxChr->getFitness())
                {
                    if (!maxChr->isLocal()) delete maxChr;
                    maxChr = c;
                }
            }
        }
    }

    return maxChr;
}

/* ------------------------------------------------------------------------ */

Grid CellularGA::getRandomWalkNewPosition(Grid p)
{
    int rnd = (int) (4.0 * (rand() / (RAND_MAX + 1.0)));

    while (p.x + randomWalkXMap[rnd] < 0 ||
           p.x + randomWalkXMap[rnd] >= params.gridXSize ||
           p.y + randomWalkYMap[rnd] < 0 ||
           p.y + randomWalkYMap[rnd] >= params.gridYSize)
    {
        rnd++;
        rnd %= 4;
    }

    return Grid(p.x + randomWalkXMap[rnd], p.y + randomWalkYMap[rnd]);
}

/* ------------------------------------------------------------------------ */

Chromosome *CellularGA::getGlobalChromosome(Grid p)
{
    assert(p.x < params.gridXSize);
    assert(p.y < params.gridYSize);

    Chromosome *ret = NULL;

#ifdef ARC_MPI_ENABLED
    double start = MPI_Wtime();
    int myRank = gMPI.getMyRank();
    int chrProcessRank = gMPI.translateGlobalPosToRank(p);

    chromosomesRequested++; /* Debug */

    /* 
     * Check if the chromosome is on my local grid or it is global. Take
     * appropriate action then.
     */
    if (chrProcessRank == myRank)
    {
        Grid local = gMPI.translateGlobalPosToLocal(p);
        ret = pop[local.x][local.y];
    }
    else
    {
        remoteChromosomesRequested++;
        ret =  getMPIChromosome(p, chrProcessRank);
    }

    double stop = MPI_Wtime();
    chrTime += stop - start;
    chrReq++;
#else
    ret = pop[p.x][p.y];
#endif /* ARC_MPI_ENABLED */

    return ret;
}

/* ------------------------------------------------------------------------ */

Grid CellularGA::translateLocalToGlobalPos(Grid p) const
{
#ifdef ARC_MPI_ENABLED
    Grid myPos = gMPI.getMyMatrixPos();
    return Grid(myPos.x * gMPI.getProcGridSize().x + p.x,
                myPos.y * gMPI.getProcGridSize().y + p.y);
#else
    return p;
#endif
}

/* ------------------------------------------------------------------------ */

void CellularGA::recalculatePopulation()
{
    Chromosome *bestChr = pop[0][0];
    fitnessSum = 0.0;

    /** Take the local grid and sum up the fitness **/
    for (int r = 0; r < lGrid.x; r++)
    {
        for (int c = 0; c < lGrid.y; c++)
        {
            fitnessSum += pop[r][c]->getFitness();
            if (pop[r][c]->getFitness() > bestChr->getFitness())
                bestChr = pop[r][c];
        }
    }

    bestIndividual = bestChr;
}

/* ------------------------------------------------------------------------ */

void CellularGA::swapIndividual(Grid p, Chromosome *c)
{
    fitnessSum -= pop[p.x][p.y]->getFitness();
    fitnessSum += c->getFitness();

    if (c->getFitness() > bestIndividual->getFitness())
        bestIndividual = c;

    /* If we want to remove the best individual */
    if (bestIndividual == pop[p.x][p.y])
        bestIndividual = c;

    delete pop[p.x][p.y];
    pop[p.x][p.y] = c;

}


/* ------------------------------------------------------------------------ */

double CellularGA::globalBestIndividualFitness()
{

/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    return globalStat.bestFitness;
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    return bestIndividualFitness();
}

/* ------------------------------------------------------------------------ */

const Chromosome *CellularGA::getGlobalBestIndividual()
{
/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    return globalStat.bestChromosome;
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    return getBestIndividual();
}

/* ------------------------------------------------------------------------ */

double CellularGA::globalAvgFitness()
{
/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    return globalStat.fitnessSum / (params.gridXSize * params.gridYSize);
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */

    return avgFitness();
}

/* ------------------------------------------------------------------------ */

/* -------------------------------------------------------------------
 * MPI communication functions
 * ------------------------------------------------------------------- */
#ifdef ARC_MPI_ENABLED
const int CellularGA::GAMPI_TAG_CHRSEND     = 0;
const int CellularGA::GAMPI_TAG_CHRREQ      = 1;
const int CellularGA::GAMPI_TAG_FSUM        = 2;
const int CellularGA::GAMPI_TAG_WHATNEXT    = 3;
const int CellularGA::GAMPI_TAG_BESTCHR     = 4;

/* ------------------------------------------------------------------------ */

Chromosome *CellularGA::getMPIChromosome(Grid p, int rank)
{
    //dbg_msg(endl);

    MPI_Request req;
    int reqCount;
    MPI_Datatype reqType;
    void *reqBuf = NULL;

    /* Create new chromosome and start MPI_Irecv() */
    Chromosome *ch = params.chrFactory->getChromosome();
    ch->getMPIProps(&reqBuf, &reqCount, &reqType);

    assert(rank < gMPI.getProcNum());
    int err;
    err = MPI_Irecv(reqBuf, reqCount, reqType, rank, GAMPI_TAG_CHRSEND,
            MPI_COMM_WORLD, &req);
    CHECK_MPI_ERROR(err);

    /*
     * Send request - can bo blocking - the other side MUST use non-blocking
     * Irecv for chromosome requests
     * */
    int buf[2];
    buf[0] = p.x;
    buf[1] = p.y;
    err = MPI_Send((void *) buf, 2, MPI_INT, rank, GAMPI_TAG_CHRREQ,
            MPI_COMM_WORLD);
    CHECK_MPI_ERROR(err);


    /* Complete non-blocking recv, and poll for requests in the loop */
    int flag = 0;
    while (!flag)
    {
        MPI_Status reqStatus;
        err = MPI_Test(&req, &flag, &reqStatus);
        CHECK_MPI_ERROR(err);
        pollMPIAction();
    }

    
    /* Chromosome received to reqBuf */

    ch->replaceWithMPIChromosome(reqBuf, reqCount);
    ch->setLocal(false);

    return ch;
}

/* ------------------------------------------------------------------------ */

void CellularGA::sendMPIChromosome(Grid p, int rank, int tag)
{
    sendMPIChromosome(pop[p.x][p.y], rank, tag);
}

/* ------------------------------------------------------------------------ */

void CellularGA::sendMPIChromosome(const Chromosome *c, int rank, int tag)
{
    const void *buf = NULL;
    int count;
    MPI_Datatype type;

    c->getMPISendData(&buf, &count, &type);

    int err;
    err = MPI_Send(const_cast<void *>(buf), count, type, rank, tag,
            MPI_COMM_WORLD);
    CHECK_MPI_ERROR(err);
}

/* ------------------------------------------------------------------------ */

void CellularGA::pollMPIAction()
{
    Grid gPos;
    int rank = -1;
    int neighSize = neighbours.getSize();
    
    for (int i = 0; i < neighSize; i++)
    {
        if (neighbours.checkRequests(gPos, rank, i))
        {
            Grid lPos = gMPI.translateGlobalPosToLocal(gPos);
            sendMPIChromosome(lPos, rank);
        }
    }
}

/* ------------------------------------------------------------------------ */

CellularGA::MPINeighbours::MPINeighbours(const GlobalMPI &g) :
    gMPI(g)
{
    int myRank = gMPI.getMyRank();

    /* Initialize Neighbour requests */
    Grid myPos = gMPI.getMyMatrixPos();

    //dbg_msg("myPos: " << myPos.y << ", " << myPos.x << ", MatrixSize: " <<
    //        gMPI.getMatrixSize() << endl);

    size = 0;
    if (myPos.y > 0)
    {
        ranks[size] = myRank - gMPI.getMatrixSize();
        size++;
    }
    if (myPos.x < gMPI.getMatrixSize() - 1)
    {
        ranks[size] = myRank + 1;
        size++;
    }
    if (myPos.y < gMPI.getMatrixSize() - 1)
    {
        ranks[size] = myRank + gMPI.getMatrixSize();
        size++;
    }
    if (myPos.x > 0)
    {
        ranks[size] = myRank - 1;
        size++;
    }

    for (int i = 0; i < size; i++)
    {
        assert(ranks[i] < gMPI.getProcNum());
        int err;
        err = MPI_Irecv(bufs[i], 2, MPI_INT, ranks[i], GAMPI_TAG_CHRREQ,
                MPI_COMM_WORLD, &requests[i]);
        CHECK_MPI_ERROR(err);
    }

}

/* ------------------------------------------------------------------------ */

bool CellularGA::MPINeighbours::checkRequests(Grid &p, int &r)
{
    int index = -1;
    int flag = 0;
    MPI_Status status;
    int err;

    err = MPI_Testany(size, requests, &index, &flag, &status);
    CHECK_MPI_ERROR(err);

    if(flag && size != 0)
    {
        CHECK_MPI_ERROR(err);

        p.x = bufs[index][0];
        p.y = bufs[index][1];
        r = ranks[index];

        assert(ranks[index] < gMPI.getProcNum());
        /* Renew request */
        err = MPI_Irecv(bufs[index], 2, MPI_INT, ranks[index],
                GAMPI_TAG_CHRREQ, MPI_COMM_WORLD, &requests[index]);
        CHECK_MPI_ERROR(err);

        return true;
    }

    return false;
}

/* ------------------------------------------------------------------------ */

bool CellularGA::MPINeighbours::checkRequests(Grid& p, int& r, int index)
{
    assert(index >= 0 && index < size);

    int flag = 0;
    MPI_Status status;
    int err;

    err = MPI_Test(&requests[index], &flag, &status);
    CHECK_MPI_ERROR(err);

    if(flag)
    {
        p.x = bufs[index][0];
        p.y = bufs[index][1];
        r = ranks[index];

        assert(ranks[index] < gMPI.getProcNum());
        /* Renew request */
        err = MPI_Irecv(bufs[index], 2, MPI_INT, ranks[index],
                GAMPI_TAG_CHRREQ, MPI_COMM_WORLD, &requests[index]);
        CHECK_MPI_ERROR(err);

        return true;
    }

    return false;
}

/* ------------------------------------------------------------------------ */

GlobalMPI::GlobalMPI(const InitParams& p) :
    MatrixSize((int) sqrt(getProcNum())),
    ProcGridSize(floor((float) p.gridXSize / MatrixSize),
                 floor((float) p.gridYSize / MatrixSize))
{
    int procNum = getProcNum();

    if (MatrixSize * MatrixSize != procNum)
    {
        ostringstream msg;
        msg << procNum << " processors cannot make square matrix.";
        throw GAException(msg);
    }
}

/* ------------------------------------------------------------------------ */

bool CellularGA::sendMPIStatistics()
{
    assert(gMPI.getMyRank() != ROOT_NODE);

    recalculatePopulation();

    /**
     * Send fitness sum and best chromosome to the root node (process 0).
     */
    int err;
    double sendBuf[2];  /** @see MPIStats::fitness description **/
    sendBuf[0] = fitnessSum;
    sendBuf[1] = bestIndividualFitness(); /* Local */


    err = MPI_Send((void *) sendBuf, 2, MPI_DOUBLE, ROOT_NODE, GAMPI_TAG_FSUM,
            MPI_COMM_WORLD);
    CHECK_MPI_ERROR(err);

    dbg_msg("Send GA stats" << endl);


    /* Wait for the continue/end command: 0 - continue, 1 - finish */
    MPI_Request request;
    MPI_Status status;
    int whatNext;       /* GAMPI_TAG_WHATNEXT context */
    int mpiFlag = 0;    /* MPI Irecv flag !! */

    for (int i = 0; i < 2; i++)
    {
        dbg_msg("Begin WHATNEXT recv" << endl);

        assert(ROOT_NODE < gMPI.getProcNum());
        err = MPI_Irecv((void *) &whatNext, 1, MPI_INT, ROOT_NODE, GAMPI_TAG_WHATNEXT,
                MPI_COMM_WORLD, &request);
        CHECK_MPI_ERROR(err);

        while (!mpiFlag)
        {
            pollMPIAction();
            err = MPI_Test(&request, &mpiFlag, &status);
            CHECK_MPI_ERROR(err);
        }

        mpiFlag = 0;

        dbg_msg("Received WHATNEXT, value: " << whatNext << endl);

        /* Determine what action to take */
        switch (whatNext)
        {
            case CONTINUE:
                return true;
                break;

            case STOP:
                return false;
                break;

            case GET_CHROMOSOME:
                assert(i != 1); /* Only once !! */

                dbg_msg("Sending chromosome to ROOT_NODE" << endl);

                /* Send my best chromosome to ROOT_NODE */
                sendMPIChromosome(getBestIndividual(), ROOT_NODE, GAMPI_TAG_BESTCHR);
                break;

            default:
                assert(0);
                break;
        }
    }

    assert(0);
    return false;
}

/* ------------------------------------------------------------------------ */

void CellularGA::recvMPIStatistics()
{
    dbg_msg("Waiting for stat" << endl);

    GAStat stat;
    stat.bestChromosome = getBestIndividual();
    stat.bestFitness = bestIndividualFitness();
    while(!mpiStats->getStat(stat)) pollMPIAction();

}

/* ------------------------------------------------------------------------ */

CellularGA::MPIStats::MPIStats(CellularGA &ga, const GlobalMPI &g, Chromosome *ch) :
    cga(ga), gMPI(g), tempChr(ch)
{
    assert(ch != NULL);
    assert(gMPI.getMyRank() == ROOT_NODE);

    int procNum = gMPI.getProcNum();

    /* Allocate receive buffers */
    fitness = new double[procNum][2];
    fitnessReqs = new MPI_Request[procNum];

    numRecv = 0;
    beginStat();
}

/* ------------------------------------------------------------------------ */

CellularGA::MPIStats::~MPIStats()
{
    delete []fitness;
    delete []fitnessReqs;

    delete tempChr;
}

/* ------------------------------------------------------------------------ */

void CellularGA::MPIStats::beginStat()
{
    assert(numRecv == gMPI.getProcNum() - 1 || numRecv == 0);

    int myRank = gMPI.getMyRank();
    int procNum = gMPI.getProcNum();
    int err;

    numRecv = 0;

    /* For each process, except me, initialize the non-blocking receive **/
    for (int i = 0; i < procNum; i++)
    {
        if (i != myRank)
        {
            assert(i < gMPI.getProcNum());
            err = MPI_Irecv(&fitness[i], 2, MPI_DOUBLE, i, GAMPI_TAG_FSUM,
                    MPI_COMM_WORLD, &fitnessReqs[i]);
            CHECK_MPI_ERROR(err);
        }
        else
        {
            fitnessReqs[i] = NULL;
        }
    }
}

/* ------------------------------------------------------------------------ */

bool CellularGA::MPIStats::getStat(GAStat &stat)
{
    int err;
    int index;
    int flag;
    MPI_Status status;
    int procNum = gMPI.getProcNum();
    int myRank = gMPI.getMyRank();


    err = MPI_Testany(procNum, fitnessReqs, &index, &flag, &status);
    CHECK_MPI_ERROR(err);

    /* Don't do any check if there is only one process */
    if (flag && procNum != 1)
        numRecv++;

    if (numRecv == procNum - 1)
    {
        dbg_msg("Finished gathering stats" << endl);

        /* Insert my best fitness in the fitness field */
        fitness[myRank][0] = stat.fitnessSum;
        fitness[myRank][1] = stat.bestFitness;
        cga.globalStat.fitnessSum = 0.0;

        //double fitnessSum = 0.0;
        double maxFitness = fitness[0][1]; /** See fitness member array **/
        double maxFIndex = 0;

        /* Get maximal fitness - including myself */
        for (int i = 1; i < procNum; i++)
        {
            cga.globalStat.fitnessSum += fitness[i][0];
            if (fitness[i][1] > maxFitness)
            {
                maxFitness = fitness[i][1];
                maxFIndex = i;
            }
        }

        if (maxFIndex == myRank)
        {
            dbg_msg("Local chromosome is the best" << endl);

            /* Set global stat */
            cga.globalStat.replaceBest(stat.bestChromosome);
            cga.globalStat.bestFitness = stat.bestFitness;
            /* fitnessSum processed in the cycle above */
            /* * * * * * * * * */

            int msg;
            if (cga.stopCondition())
                msg = STOP;
            else
                msg = CONTINUE;

            for (int rank = 0; rank < procNum; rank++)
            {
                if (rank != myRank)
                {
                    int err;
                    err = MPI_Send((void *) &msg, 1, MPI_INT, rank, GAMPI_TAG_WHATNEXT,
                            MPI_COMM_WORLD);
                    CHECK_MPI_ERROR(err);
                }
            }
        }
        else
        {
            dbg_msg("Remote chromosome is the best, retreiving..." << endl);

            /* The best chromosome is not local - get it */
            int err;
            int buf = GET_CHROMOSOME;
            err = MPI_Send((void *) &buf, 1, MPI_INT, maxFIndex,
                    GAMPI_TAG_WHATNEXT, MPI_COMM_WORLD);
            CHECK_MPI_ERROR(err);

            void *recvBuf = NULL;
            int recvCount;
            MPI_Datatype recvType;
            MPI_Status recvStatus;

            tempChr->getMPIProps(&recvBuf, &recvCount, &recvType);
            MPI_Recv(recvBuf, recvCount, recvType, maxFIndex, GAMPI_TAG_BESTCHR,
                    MPI_COMM_WORLD, &recvStatus);

            /* Set global stat */
            delete cga.globalStat.bestChromosome;
            cga.globalStat.bestChromosome = tempChr->getNewFromMPI(recvBuf,
                    recvCount);
            cga.globalStat.bestChromosome->setLocal(false);
            cga.globalStat.bestFitness = stat.bestFitness;
            /* fitnessSum processed in the cycle above */
            /* * * * * * * * * */

            dbg_msg("Chromosome retreived from process " << maxFIndex << endl);

            /* Check if the stop condition is met */
            /* TODO: stopCondtion() */
            int m;
            if (cga.stopCondition())
                m = STOP;
            else
                m = CONTINUE;

            dbg_msg("Sending 2nd WHATNEXT to all: " << m << endl);

            /* Send final command */
            for (int rank = 0; rank < procNum; rank++)
            {
                if (rank != myRank)
                {
                    int err;
                    err = MPI_Send((void *) &m, 1, MPI_INT, rank,
                            GAMPI_TAG_WHATNEXT, MPI_COMM_WORLD);
                    CHECK_MPI_ERROR(err);
                }
            }
        }

        beginStat();

        return true;
    }
    else
    {
        return false;
    }
}

/* ------------------------------------------------------------------------ */

#endif /* ARC_MPI_ENABLED */


/* End of file CellularGA.cc */
