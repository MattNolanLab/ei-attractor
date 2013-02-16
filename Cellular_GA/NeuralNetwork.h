/**
 * Name:        NeuralNetwork.h
 * Project: 	BIN-2008-L Neural Network
 * Author:  	Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
 * Date:    	2008-04
 * Description: Neural network header file.
 */

#ifndef NEURALNETWORK_H
#define NEURALNETWORK_H

#include <sstream>
#include <cassert>

#include "debug.h"
#include "error.h"

/**
 * Neural Network layer interface.
 */
class NNLayer
{
  public:
    virtual ~NNLayer() = 0;

    /**
     * Make a copy of myself.
     */
    virtual NNLayer* copy(NNLayer* prev) const = 0;
    
    
    /** Get number of neurons in the layer **/
    virtual unsigned getNumItems() const = 0;

    /**
     * Bias all weights in the specified neuron.
     *
     * @param n Neuron index
     * @param bias Weights bias.
     */
    virtual void biasWeights(unsigned n, double mRate) = 0;


    /**
     * Get weight value
     *
     * @param n Neuron index.
     * @param w weight index.
     */
    virtual double getWeight(unsigned n, unsigned w) const = 0;
        
    
    /** Evaluate outputs **/
    virtual const double* eval() = 0;

    /**
     * Combine my weights with another layer
     *
     * @param l Another layer.
     */
    virtual void combineUniform(const NNLayer* l) = 0;

    /**
     * Make averages of my weigths an weights of l and save into my weigths
     *
     * @param l Another layer.
     */
    virtual void averageWeights(const NNLayer* l) = 0;

    virtual void print(std::ostream& s) const = 0;
};


/**
 * Neural network layer class
 */
class NNFullLayer : public NNLayer
{

  public:
    /**
     * Initializes the layer.
     *
     * @param items Number of neurons in the layer.
     * @param prev Pointer to the previous layer.
     * @param randomInit Initialize the weights randomly.
     */
    NNFullLayer(unsigned items, NNLayer* prev, bool randomInit);

    ~NNFullLayer();

    /**
     * Make a copy of myself.
     */
    NNLayer* copy(NNLayer* prev) const;
    

    /** Get number of neurons in the layer **/
    unsigned getNumItems() const {return items;}

    /** Evaluate outputs **/
    const double* eval();

    /** Initialize the weights randomly **/
    void randomInit();

    /**
     * Bias all weights in the specified neuron.
     *
     * @param n Neuron index
     * @param bias Weights bias.
     */
    void biasWeights(unsigned n, double mRate);

    /**
     * Multiply all weights in neuron by a constant
     * 
     * @param n Neuron index
     * @param m Multiplication constant
     */
    void multiplyWeights(unsigned n, double m);


    /**
     * Change weight of a neuron.
     *
     * @param n Neuron index.
     * @param w Weight index.
     * @param val New value;
     */
    //void changeWeight(unsigned n, unsigne w, double val);
    
    /**
     * Combine my weights with another layer
     *
     * @param l Another layer.
     */
    void combineUniform(const NNLayer* l);

    /**
     * Make averages of my weigths an weights of l and save into my weigths
     *
     * @param l Another layer.
     */
    void averageWeights(const NNLayer* l);

    /**
     * Get weight value
     *
     * @param n Neuron index.
     * @param w weight index.
     */
    double getWeight(unsigned n, unsigned w) const
    {
        assert(n < items);
        assert(w < prevLayer->getNumItems() + 1);

        return weights[n][w];
    }
        

    /**
     * Print weights into stream
     */
    void print(std::ostream& s) const;


  protected:
      /*
       * Computes activation function.
       *
       * @param ksi Neuron potential.
       * @return Neuron output.
       */
      double computeActivationFunc(double ksi);


      /**
       * Copute neuron potential (ksi).
       *
       * @param n Neuron index
       * @param in Neuron inputs.
       * @return Neuron potential.
       */
      double computeNeuronPotential(unsigned n, const double* in);
    

  private:
    /** Previous network layer reference **/
    NNLayer* const prevLayer;

    /** Number of neurons in the layer **/
    const unsigned items;

    /** Output values array **/
    double* outputs;

    /** Neurons weights **/
    double** weights;

    /** Sigmoid activation function steepness factor **/
    static const double lambda;
};



/**
 * Input layer class.
 */
class NNInputLayer : public NNLayer
{
  public:
    NNInputLayer(unsigned numItems) : numItems(numItems), inputs(NULL) {}

    /** Load from stream **/
    //NNInputLayer(const istream& s) :
    //    inputs(NULL);
    //{
    //    if (!s.good())
    //    {
    //        std::ostringstream err;
    //        err << "Corrupted file when loading Neural Network layer" << endl;
    //        throw NNException(err);
    //    }

    //    s >> numItems
    //}

    /**
     * Make a copy of myself.
     */
    NNLayer* copy(NNLayer* prev) const
    {
        assert(prev == NULL);
        return new NNInputLayer(numItems);
    }
    
    
    /** Get number of neurons in the layer **/
    unsigned getNumItems() const {return numItems;}

    /** Evaluate outputs **/
    const double* eval() { return inputs; }

    /** Set the inputs pointer **/
    void feed(const double* i) {inputs = i;};

    /**
     * Bias all weights in the specified neuron.
     * This function will cause exception. This layer has no weights.
     *
     * @param n Neuron index
     * @param bias Weights bias.
     */
    void biasWeights(unsigned n, double mRate)
    {
        n = 0; mRate = 0; /* Fake, not used */
        std::ostringstream err;
        err << "An attempt to bias weights in an Input layer." << std::endl;
        throw NNException(err);
    }

    /**
     * Get weight value. Causes and error here.
     *
     * @param n Neuron index.
     * @param w weight index.
     */
    double getWeight(unsigned n, unsigned w) const
    {
        n = 0; w = 0; /* Fake, not used */
        std::ostringstream err;
        err << "An attempt to get weight value in an Input layer." << std::endl;
        throw NNException(err);
    }
    
    void combineUniform(const NNLayer* l)
    {
        l = NULL; /* Fake, not used */
        std::ostringstream err;
        err << "An attempt to combine input layer." << std::endl;
        throw NNException(err);
    }

    void averageWeights(const NNLayer* l)
    {
        l = NULL; /* Fake, not used */
        std::ostringstream err;
        err << "An attempt to average input layer." << std::endl;
        throw NNException(err);
    }

    void print(std::ostream& s) const
    {
        s << numItems << std::endl;
    }
    
  private:
      /** Number of items in the layer **/
      unsigned numItems;

      /** Inputs pointer **/
      const double* inputs;
};



/**
 * Neural Network class
 */
class NeuralNetwork
{

  public:
      /**
       * Constructor.
       *
       * @param numLayers Number of layers (the first layer is the network's
       *          input only).
       * @param lSizes Layer sizes. An array, with number of neurons for each
       *          layer.
       * @param random true: initialize randomly.
       */
      NeuralNetwork(unsigned numLayers, const unsigned* lSizes, bool random);

      /**
       * Make a new instantion of the neural net, using existing layers.
       *
       * @param numLayers Number of layers.
       * @param l An array of layer pointers
       */
      NeuralNetwork(unsigned numLayers, const NNLayer** l);

      ~NeuralNetwork();

      /**
       * Mape a copy of myself.
       */
      NeuralNetwork* copy() const;

      /**
       * Evaluate the network for the given input.
       *
       * @param nnInputs An array of neural network inputs. Must be the size of
       *          the first layer.
       * @return Neural network outputs. Size of the array is the size of the
       *          last layer.
       */
      const double* eval(const double* nnInputs);


      /** Get number of neurons in the first layer (network inputs) **/
      unsigned getNumInputs() const
      {
          return layers[0]->getNumItems();
      }


      /** Get numbe of neurons in the last layer (network outputs) **/
      unsigned getNumOutputs() const
      {
          return layers[numLayers - 1]->getNumItems();
      }

      /**
       * Get a layer read-only.
       *
       * @param i Layer index - bottom layer == 0
       */
      const NNLayer* getLayer(unsigned i) const
      {
          assert(i < numLayers);
          return layers[i];
      }

      /** Get number of layers **/
      unsigned getNumLayers() const { return numLayers; }

      /**
       * Get layer readWrite.
       */
      NNLayer* getLayerForChange(unsigned i)
      {
          assert(i > 0 && i < numLayers);
          return layers[i];
      }

      /**
       * Combine uniformly my weights with another neural network.
       *
       * @param n Another neutal network.
       */
      void combineUniform(const NeuralNetwork* n);

      /**
       * Average weights of me and another net n
       *
       * @param n Another neutal network.
       */
      void averageWeights(const NeuralNetwork* n);

      /**
       * Print network data into stream
       */
      void print(std::ostream& s);


      /** Get number of neurons in the network -- input neurons excluded **/
      unsigned getNumItems()
      {
          unsigned sum = 0;
          for (unsigned i = 1; i < numLayers; i++)
              sum += layers[i]->getNumItems();
          return sum;
      }

      /**
       * Bias weights of all neurons with a probability.
       *
       * @param mRate (Mutation) bias rate.
       */
      void randomBias(double mRate);

  private:

      NeuralNetwork(const NeuralNetwork& n);

      /** Number of layers in the network **/
      const unsigned numLayers;

      /** Neural network layers - the first one is the bottom one **/
      NNLayer** layers;

      /** Check if the number of layers is sane **/
      bool checkNumLayers(unsigned num)
      {
          return num > 1;
      }

};

#endif /* NEURALNETWORK_H */


/* End of file */
