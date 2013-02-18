/*
 *  place_cell_generator.cpp
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

#include "place_cell_generator.h"
#include "network.h"
#include "dict.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "exceptions.h"

using namespace nest;

/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

mynest::place_cell_generator::Parameters_::Parameters_()
  : rate_(0.0    )  // pA
{}


/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::place_cell_generator::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, nest::names::rate, rate_);
}

void mynest::place_cell_generator::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, nest::names::rate, rate_);
  if ( rate_ < 0 )
    throw BadProperty("The rate cannot be negative.");
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

mynest::place_cell_generator::place_cell_generator()
  : Node(),
    device_(),
    P_()
{}

mynest::place_cell_generator::place_cell_generator(const place_cell_generator& n)
  : Node(n),
    device_(n.device_),
    P_(n.P_)
{}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::place_cell_generator::init_state_(const Node& proto)
{
  const place_cell_generator& pr = downcast<place_cell_generator>(proto);

  device_.init_state(pr.device_);
}

void mynest::place_cell_generator::init_buffers_()
{
  device_.init_buffers();
}

void mynest::place_cell_generator::calibrate()
{
  device_.calibrate();

  // rate_ is in Hz, dt in ms, so we have to convert from s to ms
  V_.poisson_dev_.set_lambda(Time::get_resolution().get_ms() * P_.rate_ * 1e-3);
}


/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */

void mynest::place_cell_generator::update(Time const & T, const long_t from, const long_t to)
{
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  if ( P_.rate_ <= 0 )
    return;

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    if ( !device_.is_active( T + Time::step(lag) ) )
      continue;  // no spike at this lag


    librandom::RngPtr rng = net_->get_rng(get_thread());
    ulong_t n_spikes = V_.poisson_dev_.uldev(rng);

    if ( n_spikes > 0 ) // we must not send events with multiplicity 0
    {
      nest::SpikeEvent se;
      se.set_multiplicity(n_spikes);
      network()->send(*this, se, lag);
    }
  }
}

//void mynest::place_cell_generator::event_hook(DSSpikeEvent& e)
//{
//  librandom::RngPtr rng = net_->get_rng(get_thread());
//  ulong_t n_spikes = V_.poisson_dev_.uldev(rng);
//
//  if ( n_spikes > 0 ) // we must not send events with multiplicity 0
//  {
//    e.set_multiplicity(n_spikes);
//    e.get_receiver().handle(e);
//  }
//}
