/**
 * Name:        OneMaxChromosome.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: OneMax chromosome implementation
 */


#include "OneMaxChromosome.h"
#include "debug.h"



/* ------------------------------------------------------------------------ */
int OneMaxChromosome::instances = 0;

OneMaxChromosome::OneMaxChromosome(unsigned int length) :
    Chromosome()
{
    instances++;

    fitnessValid = false;
    array = new SimpleBitArray(length, -1, true);
}

/* ------------------------------------------------------------------------ */

OneMaxChromosome::OneMaxChromosome(unsigned char *d, int s) :
    Chromosome()
{
    instances++;

    fitnessValid = false;
    array = new SimpleBitArray(s, d);
}

/* ------------------------------------------------------------------------ */

OneMaxChromosome::OneMaxChromosome(const OneMaxChromosome& c)
  : Chromosome()
{
    instances++;

    array = new SimpleBitArray(*c.array);
}

/* ------------------------------------------------------------------------ */

OneMaxChromosome::~OneMaxChromosome()
{
    instances--;
    delete array;
}

/* ------------------------------------------------------------------------ */

Chromosome* OneMaxChromosome::copy() const
{
    return new OneMaxChromosome(*this);
}

/* ------------------------------------------------------------------------ */

double OneMaxChromosome::getFitness()
{
    if (fitnessValid) return fitness;

    unsigned long sum = 0;

    int size = array->getSize();
    for (int i = 0; i < size; i++)
    {
        sum += array->getBit(i);
    }

    //fitness = ((double) sum) / size;
    fitness = sum;
    fitnessValid = true;

    return fitness;
}

/* ------------------------------------------------------------------------ */

Chromosome *OneMaxChromosome::crossover(const Chromosome *c) const
{
    const OneMaxChromosome *chr = dynamic_cast<const OneMaxChromosome *>(c);

    assert(array->getSize() == chr->array->getSize());

    OneMaxChromosome *offspring = new OneMaxChromosome(array->getSize());


    /* Perform uniform crossover on the chromosomes */
    for (int i = 0; i < array->getSize(); i++)
    {
        if ((random() % 2) == 0)
            offspring->array->setBit(i, array->getBit(i));
        else
            offspring->array->setBit(i, chr->array->getBit(i));
    }
    
    ///** One point crossover **/
    //int crossoverPoint = (int) array->getSize() * rand() / (RAND_MAX + 1.0);
    //for (int i = 0; i < array->getSize(); i++)
    //{
    //    if (i < crossoverPoint)
    //        offspring->array->setBit(i, array->getBit(i));
    //    else
    //        offspring->array->setBit(i, chr->array->getBit(i));
    //}

    return offspring;
}

/* ------------------------------------------------------------------------ */

void OneMaxChromosome::mutate(double mRate)
{
    assert(mRate >= 0 && mRate <= 1.0);

    int flipped = 0;

    /* Inspiried by GALib */

    /* If the mutation rate is too small, make a bit flip one by one */
    if (mRate * array->getSize() <= 1.0)
    {
        for (int i = 0; i < array->getSize(); i++)
        {
            if (getRandomDouble() < mRate)
            {
                array->flipBit(i);
                flipped++;
            }
        }
    }

    /* Else pick up the desired number of alleles to flip */
    else
    {
        int mutNumber = mRate * array->getSize();
        for (int i = 0; i < mutNumber; i++)
        {
            array->flipBit((getRandomDouble()*array->getSize()));
            flipped++;
        }
    }

    //dbg_msg("Flipped bits: " << flipped << endl);
}

/* ------------------------------------------------------------------------ */

void OneMaxChromosome::print(std::ostream &s) const
{
    if (array->getSize() > 1000)
        return;

    for (int i = 0; i < array->getSize(); i++)
    {
        s << array->getBit(i);
    }

    s << std::endl;
}


/* ------------------------------------------------------------------------
   Begin MPI Stuff */
#ifdef ARC_MPI_ENABLED
void OneMaxChromosome::
        getMPISendData(const void **buf, int *count, MPI_Datatype *type) const
{
    *count = array->getSize();
    *type = MPI_UNSIGNED_CHAR;
    *buf = (const void *) array->getRawData();
}

/* ------------------------------------------------------------------------ */

bool OneMaxChromosome::replaceWithMPIChromosome(void *buf, int size)
{
    assert(array->getSize() == size);
    delete array;
    array = new SimpleBitArray(size, (unsigned char *) buf);

    return true;
}

/* ------------------------------------------------------------------------ */

void OneMaxChromosome::getMPIProps(void **buf, int *count, MPI_Datatype *type)
{
    *buf = (void *) new unsigned char[array->getSize()];
    *count = array->getSize();
    *type = MPI_UNSIGNED_CHAR;
}

/* ------------------------------------------------------------------------ */

Chromosome *OneMaxChromosome::getNewFromMPI(void *buf, int siz) const
{
    return new OneMaxChromosome((unsigned char *) buf, siz);
}

#endif /* ARC_MPI_ENABLED */
/* End MPI Stuff
   ------------------------------------------------------------------------ */

/* End of file OneMaxChromosme.cc */
