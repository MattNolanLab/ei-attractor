/**
 * Name:        args.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Arguments header file.
 */

#ifndef ARGS_H
#define ARGS_H

#include "CellularGA.h"


/** Parse arguments and return initialization parameters **/
InitParams parseArgs(int argc, char **argv);


/**
 * Parse float value - in the form 'switch_name' val
 *
 * in return changes cnt value
 **/
bool parseFloatArg(char **argv, const char *wantArg, int *cnt, int argc,
        float *val);


/**
 * Parse int value - in the form 'switch_name' val
 **/
bool parseIntArg(char **argv, const char *wantArg, int *cnt, int argc,
        int *val);


/**
 * Parse argument- in the form 'switch_name'
 **/
bool parseArg(char **argv, const char *wantArg, int *cnt);


/** Print error message and exit program */
void printErr(int exit_flag, const ostringstream &errMsg);


/** Checks if the counter is in the range (0, max - 1) **/
int check_counter(int counter, int max);


/** Finalize all communication (used only in MPI version) **/
void finalize();


/**
 * Print help.
 * user defined.
 */
void printHelp();


#endif /* ARGS_H */

/* End of file args.h */
