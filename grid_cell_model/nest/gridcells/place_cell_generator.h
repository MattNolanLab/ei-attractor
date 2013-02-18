/*
 *  place_cell_generator.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef PLACE_CELL_GENERATOR_H
#define PLACE_CELL_GENERATOR_H

#include "nest.h"
#include "event.h"
#include "node.h"
#include "stimulating_device.h"
#include "connection.h"
#include "poisson_randomdev.h"

namespace nest
{

  class place_cell_generator : public nest::Node
  {

  public:

    /**
     * The generator is threaded, so the RNG to use is determined
     * at run-time, depending on thread.
     */
    place_cell_generator();
    place_cell_generator(place_cell_generator const&);


    using Node::event_hook;

    nest::port check_connection(nest::Connection&, nest::port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  private:

    void init_state_(const Node&);
    void init_buffers_();
    void calibrate();

    void update(nest::Time const &, const nest::long_t, const nest::long_t);
//    void event_hook(nest::DSSpikeEvent&);

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

    nest::StimulatingDevice<nest::SpikeEvent> device_;
    Parameters_ P_;
    Variables_  V_;

  };

  inline
  nest::port place_cell_generator::check_connection(nest::Connection& c, nest::port receptor_type)
  {
    nest::SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  void place_cell_generator::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    device_.get_status(d);
  }

  inline
  void place_cell_generator::set_status(const DictionaryDatum &d)
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
