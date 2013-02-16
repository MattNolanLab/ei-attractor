/**
 * Name:        NNTrainData.cc
 * Project: 	ARC-2008-L Diffuse GA
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural network train data class.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <cassert>

#include "NNTrainData.h"
#include "error.h"
#include "debug.h"


using namespace std;


NNTrainData::NNTrainData(const char *bitmapF, const char *charF) :
    chars("")
{

    /* Read character map */
    ifstream charFile(charF);

    if (!charFile.is_open())
    {
        ostringstream err;
        err << "Could not open character file: " << charF << endl;
        throw NNException(err);
    }

    char c;
    while (charFile.good())
    {
        charFile >> noskipws >> c;
        if (c != '\n')
            chars += c;
        else break;
    }

    if (chars.length() == 0)
    {
        ostringstream err;
        err << "Character file doesn't contain any characters on the first line" << endl;
        throw NNException(err);
    }

    /* Translate characters into NN outputs - 7 bits used */
    nnOutputs = new double*[chars.length()];

    /* Translate each char to the NN output */
    for (unsigned int i = 0; i < chars.length(); i++)
    {
        nnOutputs[i] = new double[NN_OUTPUT_WIDTH];
        unsigned char c = chars[i];

        for (int cbit = 0; cbit < NN_OUTPUT_WIDTH; cbit++)
        {
            nnOutputs[i][cbit] = (double)(c & 1);
            c >>= 1;
        }
    }

    /** Fill in the NN inputs */
    unsigned imageWidth;
    unsigned imageHeight;
    unsigned char* imageData = 
        getImageData(bitmapF, &imageWidth, &imageHeight);
    //tresholdNormalize(imageData, imageWidth*imageHeight);

    charWidth = imageWidth / chars.length();
    charHeight = imageHeight;

    if (charWidth*chars.length() != imageWidth)
    {
        cerr <<
            "Warning: Cannot divide image width " << imageWidth <<
            " into " << chars.length() << " characters" << endl;
    }

    nnInputs = new double*[chars.length()];
    for (unsigned i = 0; i < chars.length(); i++)
    {
        nnInputs[i] = new double[charWidth*charHeight];
    }

    for (unsigned r = 0; r < imageHeight; r++)
    {
        for (unsigned c = 0; c < imageWidth; c++)
        {
            //dbg_msg("item: " << c/charWidth <<
            //        ", item offset: " << (r*charWidth) + (c %charWidth) << endl);

            nnInputs[c/charWidth][(r*charWidth) + (c % charWidth)] = 
                (double)(imageData[r*imageWidth + c]) / 256.0f;
        }
    }

    dbg_msg("imageWidth: " << imageWidth <<
            ", imageHeight: " << imageHeight <<
            ", charWidth: " << charWidth <<
            ", charHeight: " << charHeight << endl);


    numItems = chars.length();
    numUsedItems = numItems;
    delete []imageData;
}

/* ------------------------------------------------------------------------ */

unsigned char* NNTrainData::getImageData(const char* filePath, unsigned* iw,
        unsigned* ih)
{
    ifstream is(filePath, ifstream::in | ifstream::binary);
    ostringstream err;

    if (!is.is_open())
    {
        err << "Could not open binary image file: " << filePath << " for reading"
            << endl;
        throw NNException(err);
    }


    /* Read BMP header */
    char header[54];
    is.read(header, 54);

    if (is.gcount() != 54)
    {
        is.close();
        err << "Binary image file: " << filePath << " doesn't contain valid " <<
               "BMP header" << endl;
        throw NNException(err);
    }

    /* Check magic number */
    if (header[0] != 0x42 || header[1] != 0x4d)
    {
        /* Not a BMP file */
        is.close();
        err << "Binary image file: " << filePath << " is not a BMP file" << endl;
        throw NNException(err);
    }

    /* Check bit depth - should be 8 */
    /* And check color pallete, should be 256 */
    uint16_t bitDepth     = *((uint16_t* ) &header[28]);
    uint32_t colorPalette = *((uint32_t* ) &header[46]);
    if (bitDepth != 8 || colorPalette != 256)
    {
        is.close();
        err << "Binary image file: " << filePath << " is not grayscale" << endl;
        throw NNException(err);
    }

    /* Skip color palette */
    is.seekg(256*4, ios_base::cur);

    /* Get data */
    uint32_t imageWidth  = *((uint32_t* ) &header[18]);
    uint32_t imageHeight = *((uint32_t* ) &header[22]);

    int size = imageWidth*imageHeight;
    char* imageData = new char[size];
    char* dataPointer = imageData;      /* help pointer */

    for (unsigned r = 0; r < imageHeight; r++)
    {
        is.read(dataPointer, imageWidth);

        if (is.gcount() != (int) imageWidth)
        {
            is.close();
            err << "Binary image file: " << filePath << " corrupted" << endl;
            throw NNException(err);
        }

        if (imageWidth % 4 != 0)
            is.seekg(4 - (imageWidth % 4), ios_base::cur);

        dataPointer += imageWidth;
    }

    is.close();

    *iw = imageWidth;
    *ih = imageHeight;

    //for (int r = imageHeight - 1; r >= 0; r--)
    //{
    //    for (unsigned c = 0; c < imageWidth; c++)
    //    {
    //        if ((unsigned char)imageData[r*imageWidth + c] < 128)
    //            cerr << "X";
    //        else
    //            cerr << " ";
    //    }
    //    cerr << endl;
    //}

    return (unsigned char* )imageData;
}

/* ------------------------------------------------------------------------ */

void NNTrainData::printNNInput(unsigned i, ostream& s) const
{
    assert(i < numItems);
    //static const int treshold = 128;

    s << "Character: " << chars[i] << endl;

    for (int r = charHeight - 1; r >= 0; r--)
    {
        for (unsigned c = 0; c < charWidth; c++)
        {
            if (nnInputs[i][r*charWidth + c] < 0.5)
                s << "X";
            else
                s << " ";
        }

        s << endl;
    }

    /* Print desired output */
    for (int o = NN_OUTPUT_WIDTH - 1; o >= 0; o--)
        s << nnOutputs[i][o] << " ";

    s << endl;
}

/* ------------------------------------------------------------------------ */

void NNTrainData::tresholdNormalize(unsigned char* data, unsigned size)
{
    for (unsigned i = 0; i < size; i++)
    {
        if (data[i] < 128)
            data[i] = 0;
        else
            data[i] = 1;
    }
}

/* ------------------------------------------------------------------------ */

NNTrainData::NNTrainData(const NNTrainData* td, unsigned numInputs,
        double** inputs, unsigned inputSize)
    : chars(td->chars)
{
    assert(td->chars.length() == numInputs);

    /** Fill nnInputs **/
    nnInputs = new double*[chars.length()];
    for (unsigned c = 0; c < chars.length(); c++)
    {
        nnInputs[c] = new double[inputSize];
        for (unsigned i = 0; i < inputSize; i++)
            nnInputs[c][i] = inputs[c][i];
    }

    /* Fake char dims - TODO */
    charWidth = 1;
    charHeight = inputSize;

    /* Fill nnOutputs */
    nnOutputs = new double*[chars.length()];
    for (unsigned c = 0; c < chars.length(); c++)
    {
        nnOutputs[c] = new double[NN_OUTPUT_WIDTH];
        for (unsigned o = 0; o < NN_OUTPUT_WIDTH; o++)
            nnOutputs[c][o] = td->nnOutputs[c][o];
    }

    numItems = td->numItems;
    numUsedItems = numItems;
}

/* ------------------------------------------------------------------------ */

/* End of file NNTrainData.cc */
