/**
 * Name:        BitArray.h
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Bit array manipulation class.
 */

#ifndef BITARRAY_H
#define BITARRAY_H

#include <cstdlib>
#include <climits>


/*
 * Define the name of the type to implement the bit array.
 * Define the number of implemented bits to make macros more readable 
 */
#define ITYPE unsigned



/**
 * Class for bit field manipulation.
 */
class BitArray
{

  public:
    /**
     * Constructor.
     * Creates a bitfield with specified length.
     */
    BitArray(unsigned long size) : size(size, bool randomize = false)
    {
        const arraySize = ((size)/(IBITS) + ((size)%(IBITS) != 0));
        dat = new ITYPE[arraySize];
        
        if (randomize)
        {
            for (unsigned long i = 0; i < arraySize; i++)
            {
                data[i] = random();
            }
        }
        else
        {
            data[i] = 0;
        }
    }


    /**
     * Get bit at specified position
     *
     * @return Specified bits, 0 or 1.
     * @TODO Array boundaries check.
     */
    int getBit(unsigned long i)
    {
        //if( i >= size )
        //  Error("Index %i mimo rozsah 0..%lu\n",i ,bf->size - 1);
          
        return (dat[i/IBITS] >> i%IBITS) & 1;
    }


    /**
     * Set bit at specified position
     *
     * @param i Index
     * @param val value: 0 stands for 0, other value for 1
     * @TODO Array boundaries check.
     */
    void setBit(unsigned long i, int val) {
        //if( i >= bf->size )
        //  Error("Index %i mimo rozsah 0..%lu\n",i , bf->size - 1);
     
        if( val == 0  )
            dat[i/IBITS] &= ~(1 << i%IBITS);
        else
            dat[i/IBITS] |= (1 << i%IBITS);
    }

    /**
     * Return the length of the field.
     */
    unsigned long getSize() { return length; }

  private:
    unsigned long size;
    ITYPE *dat;

    /** Size of the type used as an array **/
    static const TYPE_SIZE;
    /** Number of bits in this array type used **/
    static const IBITS;

}

BitArray::TYPE_SIZE = sizeof(ITYPE);
BitArray::IBITS = (TYPE_SIZE)*CHAR_BIT;


#endif /* BITARRAY_H */


/* End of file BitArray.h */
