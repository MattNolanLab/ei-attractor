/**
 * Name:        NNTrainData.h
 * Project: 	BIN-2008-L Neural network trained using GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Loads an manipulates image train data.
 */

#ifndef NNTRAINDATA_H
#define NNTRAINDATA_H

#include <iostream>

#include "NeuralNetwork.h"
#include "error.h"

#define NN_OUTPUT_WIDTH 7


/**
 * Neural network train data.
 */
class NNTrainData
{

  public:

    /**
     * Loads training data from two files, bitmap and corresponding characters
     * for each bitmap glyph.
     */
    NNTrainData(const char *bitmapFile, const char *charFile);

    /**
     * Create new training data, use current outputs, but instead of current
     * inputs, use given inputs.
     */
    NNTrainData(const NNTrainData* td, unsigned numInputs, double** inputs,
            unsigned inputSize);


    ~NNTrainData()
    {
        for (unsigned i = 0 ; i < numItems; i++)
        {
            delete []nnInputs[i];
            delete []nnOutputs[i];
        }
        delete []nnInputs;
        delete []nnOutputs;
    }


    /**
     * Get NN inputs.
     */
    const double* getNNInputs(unsigned i) const {return nnInputs[i];}

    /** Get NN output **/
    const double* getNNOutputs(unsigned i) const {return nnOutputs[i];}

    /** Get character width **/
    unsigned getCharWidth() const {return charWidth;}

    /** Get character height **/
    unsigned getCharHeight() const {return charHeight;}

    /** Get number of items **/
    unsigned getNumItems() const { return numItems; }

    /** Get training char **/
    char getChar(unsigned i) const { return chars[i]; }

    /** Print neural network input in ascii into the stream s**/
    void printNNInput(unsigned i, std::ostream& s) const;

    /** Set number of items taken for training **/
    void setNumUsedItems(unsigned i)
    {
        if (i > numItems)
        {
            std::ostringstream err;
            err << "Attempt to use " << i << " items for training, " 
                << "but only " << numItems << " available" << std::endl;
            throw NNException(err);
        }
        numUsedItems = i;
    }

    /** Get number of items used for training **/
    unsigned getNumUsedItems() const { return numUsedItems; }

  protected:

    /**
     * Loads data from bitmap file. Bitmap _must_ be grayscale.
     *
     * @param bitmapF BMP file path.
     * @param iw      Image width;
     * @param ih      Image height.
     * @return Image data
     */
    unsigned char* getImageData(const char* filePath, unsigned* iw,
            unsigned* ih);
    
    /**
     * Apply thresholding on image data and normalize.
     *
     * @param data Image data
     * @param size Size of the data in bytes.
     */
    void tresholdNormalize(unsigned char* data, unsigned size);
  private:

    /** Array of Neural network inputs **/
    double** nnInputs;
    /** Array of Neural network outputs **/
    double** nnOutputs;

    /** Number of training items **/
    unsigned numItems;

    /** Number of used training items **/
    unsigned numUsedItems;

    /** Character width (px) **/
    unsigned charWidth;

    /** Character height **/
    unsigned charHeight;

    /** Actual characters for training **/
    std::string chars;
};

#endif /* NNTRAINDATA_H */

/* End of file NNTrainData.h */
