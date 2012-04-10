/*
 *  aeif_cond_exp_custom.h
 *
 *  This file is part of NEST
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

#ifndef AEIF_COND_EXP_H
#define AEIF_COND_EXP_H

#include "config.h"

#ifdef HAVE_GSL_1_11

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

/* BeginDocumentation
Name: aeif_cond_exp_custom - Conductance based exponential integrate-and-fire neuron model according to Brette and Gerstner (2005).

Description:
aeif_cond_exp_custom is the adaptive exponential integrate and fire neuron according to Brette and Gerstner (2005).

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg solver with adaptive stepsize to integrate
the differential equation.

The membrane potential is given by the following differential equation:
C dV/dt= -g_L(V-E_L)+g_L*Delta_T*exp((V-V_T)/Delta_T)-g_e(t)(V-E_e) -g_i(t)(V-E_i)-w +I_e

and

tau_w * dw/dt= a(V-E_L) -W

Parameters: 
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m        double - Membrane potential in mV
  g_ex       double - Excitatory synaptic conductance in nS.
  g_in       double - Inhibitory synaptic conductance in nS.
  w          double - Spike-adaptation current in pA.

Membrane Parameters:
  C_m        double - Capacity of the membrane in pF
  t_ref      double - Duration of refractory period in ms. 
  V_peak     double - Spike detection threshold in mV.
  V_reset    double - Reset value for V_m after a spike. In mV.
  E_L        double - Leak reversal potential in mV. 
  g_L        double - Leak conductance in nS.
  I_e        double - Constant external input current in pA.

Spike adaptation parameters:
  a          double - Subthreshold adaptation in nS.
  b          double - Spike-triggered adaptation in pA.
  Delta_T    double - Slope factor in mV
  tau_w      double - Adaptation time constant in ms
  V_t        double - Spike initiation threshold in mV (V_th can also be used for compatibility).

Synaptic parameters
  E_ex       double - Excitatory reversal potential in mV.
  tau_syn_ex double - Rise time of excitatory synaptic conductance in ms (exp function).
  E_in       double - Inhibitory reversal potential in mV.
  tau_syn_in double - Rise time of the inhibitory synaptic conductance in ms (exp function).


Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: Brette R and Gerstner W (2005) Adaptive Exponential Integrate-and-Fire Model as 
            an Effective Description of Neuronal Activity. J Neurophysiol 94:3637-3642

SeeAlso: iaf_cond_exp, aeif_cond_alpha
*/

namespace nest
{
  namespace names
  {
    /**
     * Add state and parameter names to the nest::names namespace
     */
    const Name E_AMPA;
    const Name E_NMDA;
    const Name E_GABA_A;
    const Name tau_AMPA_fall;
    const Name tau_NMDA_rise;
    const Name tau_NMDA_fall;
    const Name tau_GABA_A_rise;
    const Name tau_GABA_A_fall;
    const Name tau_AHP;

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
  int aeif_cond_exp_custom_dynamics (double, const double*, double*, void*);

  class aeif_cond_exp_custom:
    public Archiving_Node
  {
    
  public:        
    
    aeif_cond_exp_custom();
    aeif_cond_exp_custom(const aeif_cond_exp_custom&);
    ~aeif_cond_exp_custom();

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
    enum SynapseTypes { INF_SPIKE_RECEPTOR = 0,
        AMPA, NMDA, AMPA_NMDA, GABA_A,
        SUP_SPIKE_RECEPTOR };
    
    void init_node_(const Node &proto);
    void init_state_(const Node &proto);
    void init_buffers_();
    void calibrate();
    void update(const Time &, const long_t, const long_t);

    // make dynamics function quasi-member
    friend int aeif_cond_exp_custom_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<aeif_cond_exp_custom>;
    friend class UniversalDataLogger<aeif_cond_exp_custom>;


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
        G_AMPA,
        G_NMDA,
        G_GABA_A,
        STATE_VEC_SIZE
      };


      double_t    y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
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
      Buffers_(aeif_cond_exp_custom &);                    //!<Sets buffer pointers to 0
      Buffers_(const Buffers_ &, aeif_cond_exp_custom &);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      UniversalDataLogger<aeif_cond_exp_custom> logger_;

      /** buffers and sums up incoming spikes/currents */
      RingBuffer spike_exc_;
      RingBuffer spike_inh_;
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
    double_t get_y_elem_() const { return S_.y_[elem]; }

    // ---------------------------------------------------------------- 

    Parameters P;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static RecordablesMap<aeif_cond_exp_custom> recordablesMap_;
  };

  inline  
  port aeif_cond_exp_custom::check_connection(Connection &c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port aeif_cond_exp_custom::connect_sender(SpikeEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }

  inline
  port aeif_cond_exp_custom::connect_sender(DataLoggingRequest& dlr, 
				     port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  port aeif_cond_exp_custom::connect_sender(CurrentEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  void aeif_cond_exp_custom::get_status(DictionaryDatum &d) const
  {
    P.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    (*d)[names::recordables] = recordablesMap_.get_list();
  }

  inline
  void aeif_cond_exp_custom::set_status(const DictionaryDatum &d)
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
