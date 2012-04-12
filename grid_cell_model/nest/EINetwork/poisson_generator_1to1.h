/*
 *  poisson_generator_1to1.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef POISSON_GENERATOR_1TO1_H
#define POISSON_GENERATOR_1TO1_H

#include "nest.h"
#include "event.h"
#include "node.h"
#include "stimulating_device.h"
#include "connection.h"
#include "poisson_randomdev.h"

namespace nest
{
/*! Class poisson_generator simulates a large population
    of randomly (Poisson) firing neurons. It replaces the old
    neuron-intrinsic shot-noise generator
*/


/*BeginDocumentation
Name: poisson_generator - simulate neuron firing with Poisson processes statistics.
Description:
  The poisson_generator simulates a neuron that is firing with Poisson statistics,
  i.e. exponentially distributed interspike intervals. It will generate a _unique_
  spike train for each of it's targets. If you do not want this behavior and need
  the same spike train for all targets, you have to use a parrot neuron inbetween
  the poisson generator and the targets.

Parameters:
   The following parameters appear in the element's status dictionary:

   rate     double - mean firing rate in Hz
   origin   double - Time origin for device timer in ms
   start    double - begin of device application with resp. to origin in ms
   stop     double - end of device application with resp. to origin in ms

Sends: SpikeEvent

Remarks:
   A Poisson generator may, especially at high rates, emit more than one
   spike during a single time step. If this happens, the generator does
   not actually send out n spikes. Instead, it emits a single spike with
   n-fold synaptic weight for the sake of efficiency.

   The design decision to implement the Poisson generator as a device
   which sends spikes to all connected nodes on every time step and then
   discards the spikes that should not have happened generating random
   numbers at the recipient side via an event hook is twofold.

   On one hand, it leads to the saturation of the messaging network with
   an enormous amount of spikes, most of which will never get delivered
   and should not have been generated in the first place.

   On the other hand, a proper implementation of the Poisson generator
   needs to provide two basic features: (a) generated spike trains
   should be IID processes w.r.t. target neurons to which the generator
   is connected and (b) as long as virtual_num_proc is constant, each
   neuron should receive an identical Poisson spike train in order to
   guarantee reproducibility of the simulations across varying machine
   numbers.

   Therefore, first, as network()->send sends spikes to all the
   recipients, differentiation has to happen in the hook, second, the
   hook can use the RNG from the thread where the recipient neuron sits,
   which explains the current design of the generator. For details,
   refer to:

   http://ken.brainworks.uni-freiburg.de/cgi-bin/mailman/private/nest_developer/2011-January/002977.html

SeeAlso: poisson_generator_ps, Device, parrot_neuron
*/

  class poisson_generator_1to1 : public Node
  {

  public:

    /**
     * The generator is threaded, so the RNG to use is determined
     * at run-time, depending on thread.
     */
    poisson_generator_1to1();
    poisson_generator_1to1(poisson_generator_1to1 const&);

    bool has_proxies() const {return false;}


    port check_connection(Connection&, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  private:

    void init_node_(const Node&);
    void init_state_(const Node&);
    void init_buffers_();
    void calibrate();

    void update(Time const &, const long_t, const long_t);

    // ------------------------------------------------------------

    /**
     * Store independent parameters of the model.
     */
    struct Parameters_ {
      double_t rate_;   //!< process rate in Hz

      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };

    // ------------------------------------------------------------

    struct Variables_ {
      librandom::PoissonRandomDev poisson_dev_;  //!< Random deviate generator
    };

    // ------------------------------------------------------------

    StimulatingDevice<SpikeEvent> device_;
    Parameters_ P_;
    Variables_  V_;

  };

  inline
  port poisson_generator_1to1::check_connection(Connection& c, port receptor_type)
  {
    DSSpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  void poisson_generator_1to1::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    device_.get_status(d);
  }

  inline
  void poisson_generator_1to1::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty

    // We now know that ptmp is consistent. We do not write it back
    // to P_ before we are also sure that the properties to be set
    // in the parent class are internally consistent.
    device_.set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
  }

} // namespace nest

#endif
