/**
 * Name:        error.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Error functions definitions.
 */

#include <iostream>
#include <cstdlib>

#include "error.h"


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED

/**
 * Print MPI error message to cerr
 */
void printMPIError(int errNum)
{
#ifdef NO_MPI_ERROR_PRINT
    cerr << "MPI error occured (MPI_Error_string() disabled), MPI error code: "
        << errNum << endl;
#else
    static char errString[STRERR_LEN_MPI] = {'\0'};
    int errLen = 0;

    MPI_Error_string(errNum, errString, &errLen);
    std::cerr << "MPI error occured. MPI Error message: " << errString;

    exit(EXIT_FAILURE);
#endif /* NO_MPI_ERROR_PRINT */
}

#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */

/* End of file error.cc */
