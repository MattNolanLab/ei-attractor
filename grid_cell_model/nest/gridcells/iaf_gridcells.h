/*
 *  iaf_gridcells.h
 *
 *  This file has been created from the NEST aeif_cond_exp native model.
 *  Modified by Lukas Solanka <l.solanka@sms.ed.ac.uk> to incorporate it in the grid cell simultaions.
 *
 *  Copyright (C) 2010 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef IAF_GRIDCELLS_H
#define IAF_GRIDCELLS_H

#include "config.h"

#ifdef HAVE_GSL_1_11

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

/* 
iaf_gridcells - Conductance based exponential integrate-and-fire neuron. Forward Euler integration method.

Description:
iaf_gridcells is a non-adaptive, exponential integrate and fire neuron. Refractory properties can be set by AHP currents. The model can receive AMPA, NMDA, and GABA_A synapses.

The membrane potential is given by the following differential equation:
C dV/dt= -g_L(V-E_L)+g_L*Delta_T*exp((V-V_T)/Delta_T) - g_XXXX(t)(V-E_XXXX) + I_e

where

XXXX is all of:
  * AMPA
  * NMDA
  * GABA_A

AMPA conductances are modeled as a single exponential. NMDA, GABA_A conductances are modeled as a difference of exponentials, peaked at the incoming weight values

Parameters: 
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m               double - Membrane potential in mV
  g_AMPA            double - AMPA synaptic conductance in nS.
  g_NMDA            double - NMDA synaptic conductance in nS.
  g_GABA_A          double - GABA_A synaptic conductance in nS.

Membrane Parameters:
  V_peak            double - Spike detection threshold in mV.
  V_reset           double - Reset value for V_m after a spike. In mV.
  t_ref             double - Duration of refractory period in ms. 
  g_L               double - Leak conductance in nS.
  C_m               double - Capacity of the membrane in pF
  E_L               double - Leak reversal potential in mV. 
  E_AHP             double - Reversal potential of afterhyperpolarisation in mV (AHP)
  tau_AHP           double - decay time of AHP in ms (exp function)
  I_e               double - Constant external input current in pA.

Spike initiation parameters:
  Delta_T           double - Slope factor in mV
  V_th              double - Spike initiation threshold in mV

Synaptic parameters
  E_AMPA            double - AMPA reversal potential in mV
  E_NMDA            double - NMDA reversal potential in mV
  E_GABA_A          double - GABA_A reversal potential in mV

  tau_AMPA_fall     double - decay time of AMPA synaptic conductance in ms
  tau_NMDA_fall     double - decay time of NMDA synaptic conductance in ms
  tau_NMDA_rise     double - rise time of NMDA synaptic conductance in ms
  tau_GABA_A_fall   double - decay time of GABA_A synaptic conductance in ms
  tau_GABA_A_rise   double - rise time of GABA_A synaptic conductance in ms


Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

SeeAlso: iaf_cond_exp, aeif_cond_alpha
*/

namespace nest
{
  namespace names
  {
    /**
     * Add state and parameter names to the nest::names namespace
     */
    const Name E_AMPA("E_AMPA");
    const Name E_NMDA("E_NMDA");
    const Name E_GABA_A("E_GABA_A");
    const Name tau_AMPA_fall("tau_AMPA_fall");
    //const Name tau_NMDA_rise("tau_NMDA_rise");
    const Name tau_NMDA_fall("tau_NMDA_fall");
    const Name tau_GABA_A_rise("tau_GABA_A_rise");
    const Name tau_GABA_A_fall("tau_GABA_A_fall");
    const Name tau_AHP("tau_AHP");
    const Name g_AMPA("g_AMPA");
    const Name g_NMDA("g_NMDA");
    const Name g_GABA_A("g_GABA_A");
  }


  /**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  extern "C"
  int iaf_gridcells_dynamics (double, const double*, double*, void*);

  class iaf_gridcells:
    public Archiving_Node
  {
    
  public:        
    
    iaf_gridcells();
    iaf_gridcells(const iaf_gridcells&);
    ~iaf_gridcells();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */
    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection &, port);
    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &);
    
    port connect_sender(SpikeEvent &, port);
    port connect_sender(CurrentEvent &, port);
    port connect_sender(DataLoggingRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    /**
     * Synapse types to connect to
     * @note Excluded upper and lower bounds are defined as INF_, SUP_.
     *       Excluding port 0 avoids accidental connections.
     */
    enum SynapseTypes {
        AMPA = 0,
        NMDA,
        GABA_A,
        SYNAPSE_TYPES_SIZE };
    
    void init_node_(const Node &proto);
    void init_state_(const Node &proto);
    void init_buffers_();
    void calibrate();
    void update(const Time &, const long_t, const long_t);

    // make dynamics function quasi-member
    friend int iaf_gridcells_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<iaf_gridcells>;
    friend class UniversalDataLogger<iaf_gridcells>;


  private:
    // ---------------------------------------------------------------- 

    //! Independent parameters
    struct Parameters
    {
      double_t V_peak;     //!< Spike detection threshold in mV
      double_t V_reset;    //!< Reset Potential in mV
      double_t t_ref;      //!< Refractory period in ms

      double_t g_L;         //!< Leak Conductance in nS
      double_t C_m;         //!< Membrane Capacitance in pF
      double_t E_L;         //!< Leak reversal Potential (aka resting potential) in mV
      double_t E_AMPA;
      double_t E_NMDA;
      double_t E_GABA_A;
      double_t Delta_T;     //!< Slope faktor in ms.
      double_t V_th;        //!< Spike threshold in mV.
      double_t I_e;         //!< Intrinsic current in pA.

      double_t tau_AMPA_fall;
      double_t tau_NMDA_rise;
      double_t tau_NMDA_fall;
      double_t tau_GABA_A_rise;
      double_t tau_GABA_A_fall;
  
      double_t tau_AHP;
  
      Parameters();  //!< Sets default parameter values

      void get(DictionaryDatum &) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum &);  //!< Set values from dicitonary
    };

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_
    {
      /**
       * Enumeration identifying elements in state array State_::y_.
       * The state vector must be passed to GSL as a C array. This enum
       * identifies the elements of the vector. It must be public to be
       * accessible from the iteration function.
       */  
      enum StateVecElems
      {
        V_M   = 0,
        // All synaptic conductances that receive spikes must be the last in the list
        G_AMPA,
        G_NMDA,
        G_GABA_A,
        STATE_VEC_SIZE
      };


      double_t       y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
      double_t    dydt_[STATE_VEC_SIZE];  //!<  derivatives
      int_t       r_;                  //!< number of refractory steps remaining

      State_(const Parameters &);      //!< Default initialization
      State_(const State_ &);
      State_& operator = (const State_ &);

      void get(DictionaryDatum &) const;
      void set(const DictionaryDatum &, const Parameters &);
    };

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_
    {
      Buffers_(iaf_gridcells &);                    //!<Sets buffer pointers to 0
      Buffers_(const Buffers_ &, iaf_gridcells &);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      UniversalDataLogger<iaf_gridcells> logger_;

      /** buffers and sums up incoming spikes/currents */
      std::vector<RingBuffer> spike_inputs_;
      RingBuffer currents_;

      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      double_t step_;             //!< step size in ms
      double   IntegrationStep_;  //!< current integration time step, updated by GSL

      /** 
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      double_t I_stim_;
    };

    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_
    {
      int_t RefractoryCounts_;
    };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    double_t get_y_elem_() const {
        //if (elem == State_::G_GABA_A) assert(S_.y_[elem] >= 0);
        return S_.y_[elem];
    }

    // ---------------------------------------------------------------- 

    Parameters P;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static RecordablesMap<iaf_gridcells> recordablesMap_;
  };

  inline  
  port iaf_gridcells::check_connection(Connection &c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port iaf_gridcells::connect_sender(SpikeEvent &, port receptor_type)
  {
    if (receptor_type < 0 || receptor_type >= SYNAPSE_TYPES_SIZE)
      throw UnknownReceptorType(receptor_type, get_name());
    return receptor_type;
  }

  inline
  port iaf_gridcells::connect_sender(DataLoggingRequest& dlr, 
				     port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  port iaf_gridcells::connect_sender(CurrentEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  void iaf_gridcells::get_status(DictionaryDatum &d) const
  {
    P.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    DictionaryDatum receptor_type = new Dictionary();
  
    (*receptor_type)["AMPA"] = AMPA;
    (*receptor_type)["NMDA"] = NMDA;
    (*receptor_type)["GABA_A"] = GABA_A;
  
    (*d)["receptor_types"] = receptor_type;
    (*d)[names::recordables] = recordablesMap_.get_list();
  }

  inline
  void iaf_gridcells::set_status(const DictionaryDatum &d)
  {
    Parameters ptmp = P;  // temporary copy in case of errors
    ptmp.set(d);            // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);      // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    Archiving_Node::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P = ptmp;
    S_ = stmp;
  }
  
} // namespace

#endif // HAVE_GSL_1_11
#endif // AEIF_COND_EXP_H
