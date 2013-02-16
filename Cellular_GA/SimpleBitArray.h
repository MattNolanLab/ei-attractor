/**
 * Name:        SimpleBitArray.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Simple Bit array manipulation class.
 */

#ifndef SIMPLEBITARRAY_H
#define SIMPLEBITARRAY_H

#include <ostream>

#include <cassert>
#include <cstdlib>

#include "debug.h"

using namespace std;

/**
 * Class for bit field manipulation.
 */
class SimpleBitArray
{

  public:
    /**
     * Constructor.
     * Creates a bitfield with specified length.
     *
     * @param size Size of the array
     * @param initVal Initialization value. 0, 1; -1 --> random
     * @param init If true, initialize with initVal
     */
    SimpleBitArray(int size, int initVal = -1, bool init = false) : size(size)
    {
        dat = new unsigned char[size];

        
        if (init)
        {
            if (initVal == -1)
            {
                for (int i = 0; i < size; i++)
                {
                    dat[i] = (int) (2.0 * rand() / (RAND_MAX + 1.0));
                }
            }
            else if (initVal == 0 || initVal == 1)
            {
                for (int i = 0; i < size; i++)
                {
                    dat[i] = initVal;
                }
            }
            else
            {
                assert(false);
            }
        }
    }

    /**
     * Construct array from a prepared data.
     */
    SimpleBitArray(int size, unsigned char *init) : dat(init), size(size) {}


    /** 
     * Copy constructor. Makes deep copy.
     */
    SimpleBitArray(const SimpleBitArray& a)
    {
        size = a.size;
        dat = new unsigned char[size];

        for (int i = 0; i < size; i++)
            dat[i] = a.dat[i];
    }

    ~SimpleBitArray() { delete []dat; };


    /**
     * Get bit at specified position
     *
     * @return Specified bits, 0 or 1.
     * @TODO Array boundaries check.
     */
    int getBit(int i) const
    {
        assert(dat[i] == 0 || dat[i] == 1);
        
        return dat[i];
    }


    /**
     * Set bit at specified position
     *
     * @param i Index
     * @param val value, must be 0 or 1 otherwise error
     * @TODO Array boundaries check.
     */
    void setBit(int i, int val) {
        assert(val == 0 || val == 1);

        dat[i] = val;
    }


    /**
     * Flib a bit at specified position
     */
    void flipBit(int pos)
    {
        assert(pos >= 0 && pos < size);

        if (dat[pos] == 1)
            dat[pos] = 0;
        else
            dat[pos] = 1;
    }

    /**
     * Return the length of the field.
     */
    int getSize() const
    {
        return size;
    }


    /**
     * Return raw data, read-only.
     */
    const unsigned char *getRawData()
    {
        return dat;
    }

    void printDebug(ostream &s)
    {
        for (int i = 0; i < size; i++)
            s << "." << (int) dat[i];
        s << std::endl;
    }

  private:
    unsigned char *dat;
    int size;

};

#endif /* SIMPLEBITARRAY_H */


/* End of file BitArray.h */
