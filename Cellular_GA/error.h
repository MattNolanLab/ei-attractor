/**
 * Name:        error.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Error handling.
 */


#ifndef ERROR_H
#define ERROR_H

#include <exception>
#include <iostream>
#include <string>
#include <sstream>

/* ------------------------------------------------------------------------ */
/* Start of ARC MPI block */
#ifdef ARC_MPI_ENABLED

#include "mpi.h"

#define STRERR_LEN_MPI 1000

/**
 * Print MPI error message to cerr
 */
void printMPIError(int errNum);


/* Check and print error message if error occured */
#define CHECK_MPI_ERROR(errNum)         \
{                                       \
    if (errNum != MPI_SUCCESS)          \
    {                                   \
        cerr << __func__ << std::endl;  \
        printMPIError(errNum);          \
        exit(EXIT_FAILURE);             \
    }                                   \
}

#endif /* ARC_MPI_ENABLED */
/* End of ARC MPI block */
/* ------------------------------------------------------------------------ */

/**
 * Generic error exception
 */
class GenericException : public std::exception
{
  protected:
    std::string msg;

  public:
    GenericException()
    {
        msg = std::string("An unexpected error occured");
    }

    /** Initialize with specific message **/
    GenericException(const std::ostringstream& s)
    {
        msg = s.str();
    }

    ~GenericException() throw() {}

    const char* what()
    {
        return msg.c_str();
    }
};

/* ------------------------------------------------------------------------ */

/**
 * Genetic algorithm exception
 */
class GAException : public GenericException
{
  public:
    GAException() : GenericException() {}

    GAException(const std::ostringstream& s) : GenericException(s) {}

    ~GAException() throw() {}
};


/**
 * Neural Network exception
 */
class NNException : public GenericException
{
  public:
    NNException() : GenericException() {}

    NNException(const std::ostringstream& s) : GenericException(s) {}

    ~NNException() throw() {}

    const char* what()
    {
        return (std::string("Neural network error: ") + msg).c_str();
    }
};

#endif /* ERROR_H */

/* End of file */
