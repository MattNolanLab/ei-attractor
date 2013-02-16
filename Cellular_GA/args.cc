/**
 * Name:        args.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Arguments and shared stuff for the CellularGA.
 */


#include <sstream>

#include "args.h"


/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
#   include "mpi.h"
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */


/* ------------------------------------------------------------------------ */

InitParams parseArgs(int argc, char **argv)
{
    // Initialize argument counter
    int cnt = 1; 
    int intArg = -1;
    float floatArg = -1.0;

    ostringstream errMsg;

    InitParams p;

    /* Defaul paremeters */
    p.maxGen            = 0;
    p.crossoverRate     = 0.7;
    p.crossoverType     = C_RANDOMWALK;
    p.mutationRate      = 0.001;
    p.gridXSize         = 10;
    p.gridYSize         = 10;
    p.randomWalkLength  = 1;
    p.chrFactory        = NULL;  /* user defined */
    
    while (cnt < argc)
    {
        /* Parameter -h */
        if (parseArg(argv, "-h", &cnt))
        {
            printHelp();
            finalize();
            exit(EXIT_SUCCESS);
        }
    
    
        /* Parameter -g - maxGen */
        else if (parseIntArg(argv, "-g", &cnt, argc, &intArg))
        {
            p.maxGen = intArg;

            if (p.maxGen < 0)
            {
                errMsg << "Maximal number of generation must be >= 0" << endl;
                printErr(true, errMsg);
            }
        }
    
        /* Parameter -c -- crossoverRate */
        else if (parseFloatArg(argv, "-c", &cnt, argc, &floatArg))
        {
            p.crossoverRate = floatArg;

            if (p.crossoverRate <= 0 || p.crossoverRate > 1)
            {
                errMsg << "Crossover rate must be int the range (0, 1>" << endl;
                printErr(true, errMsg);
            }
        }
    
        /* Parameter -m -- mutationRate */
        else if (parseFloatArg(argv, "-m", &cnt, argc, &floatArg))
        {
            p.mutationRate = floatArg;

            if (p.mutationRate <= 0 || p.mutationRate > 1)
            {
                errMsg << "Mutation rate must be int the range (0, 1>" << endl;
                printErr(true, errMsg);
            }
        }
    
        /* Parameter -xsize -- gridXSize */
        else if (parseIntArg(argv, "-xsize", &cnt, argc, &intArg))
        {
            p.gridXSize = intArg;

            if (!(p.gridXSize > 0))
            {
                errMsg << "X size of the population grid must be > 0" << endl;
                printErr(true, errMsg);
            }
        }

        /* Parameter -ysize -- gridYSize */
        else if (parseIntArg(argv, "-ysize", &cnt, argc, &intArg))
        {
            p.gridYSize = intArg;

            if (!(p.gridYSize > 0))
            {
                errMsg << "Y size of the population grid must be > 0" << endl;
                printErr(true, errMsg);
            }
        }
    
        /* Parameter -r -- randomWalkLength/area radius */
        else if (parseIntArg(argv, "-r", &cnt, argc, &intArg))
        {
            p.randomWalkLength = intArg;

            if (!(p.randomWalkLength >= 1))
            {
                errMsg << "Random walk length/crossover area radius must be more than 0" << endl;
                printErr(true, errMsg);
            }
        }
    
    
        /* Parameter -ct -- Crossover type */
        else if (parseIntArg(argv, "-ct", &cnt, argc, &intArg))
        {
            p.crossoverType = intArg;
            if (!(intArg >= 0))
            {
                errMsg << "Invalid Crossover type" << endl;
                printErr(true, errMsg);
            }
        }
    
    
        /* Unknown parameter */
        else
        {
            ostringstream errMsg;
            errMsg << "Unknown parameter: " << argv[cnt] << endl;
            printErr(true, errMsg);
        }
    }

    /** Check crossover type **/
    if (p.crossoverType != C_RANDOMWALK &&
        p.crossoverType != C_AREA)
    {
            ostringstream errMsg;
            errMsg << "Invalid Crossover type" << endl;
            printErr(true, errMsg);
    }

    return p;
}

/* ------------------------------------------------------------------------ */

bool parseFloatArg(char **argv, const char *wantArg, int *cnt, int argc,
        float *val)
{
    bool match = false;
    char *endptr = NULL;
    ostringstream errMsg;

    if (strcmp(argv[*cnt], wantArg) == 0)
    {
        match = true;
        (*cnt)++;
    
        if (!check_counter(*cnt, argc))
        {
            errMsg << argv[*cnt - 1] << " parameter requires another argument" << endl;
            printErr(true, errMsg);
        }
    
        *val = strtof(argv[*cnt], &endptr);
    
        if (*endptr != '\0')
        {
            errMsg << "Could not parse number after " << argv[*cnt - 1] << ": " <<
                argv[*cnt] << endl;
            printErr(true, errMsg);
        }

        (*cnt)++;
    }

    return match;
}

/* ------------------------------------------------------------------------ */

bool parseIntArg(char **argv, const char *wantArg, int *cnt, int argc,
        int *val)
{
    bool match = false;
    char *endptr = NULL;
    ostringstream errMsg;

    if (strcmp(argv[*cnt], wantArg) == 0)
    {
        match = true;
        (*cnt)++;
    
        if (!check_counter(*cnt, argc))
        {
            errMsg << argv[*cnt - 1] << " parameter requires another argument" << endl;
            printErr(true, errMsg);
        }
    
        *val = strtol(argv[*cnt], &endptr, 10);
    
        if (*endptr != '\0')
        {
            errMsg << "Could not parse number after " << argv[*cnt - 1] << ": " <<
                argv[*cnt] << endl;
            printErr(true, errMsg);
        }

        (*cnt)++;
    }

    return match;
}

/* ------------------------------------------------------------------------ */

bool parseArg(char **argv, const char *wantArg, int *cnt)
{
    bool match = false;

    if (strcmp(argv[*cnt], wantArg) == 0)
    {
        match = true;
        (*cnt)++;
    }

    return match;
}

/* ------------------------------------------------------------------------ */

void printErr(int exit_flag, const ostringstream &errMsg)
{
/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    int myRank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == ROOT_NODE)
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */
    {
        cerr << errMsg.str();
    }

    if (exit_flag)
    {
        finalize();
        exit(EXIT_FAILURE);
    }
}

/* ------------------------------------------------------------------------ */

int check_counter(int counter, int max)
{
    return counter < max && counter > 0 ? 1 : 0;
}

/* ------------------------------------------------------------------------ */

/** Finalize all communication (used only in MPI version) **/
void finalize()
{
/* Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
    MPI_Finalize();
#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff */
}

/* ------------------------------------------------------------------------ */

/* End of file args.cc */
