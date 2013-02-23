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


namespace mynest 
{

vecType place_cell_generator::Parameters_::rat_pos_x(1, .0);
vecType place_cell_generator::Parameters_::rat_pos_y(1, .0);
double place_cell_generator::Parameters_::rat_pos_dt = 20.; // msec


/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

place_cell_generator::Parameters_::Parameters_()
    : rate(      0.0),  // pA
      ctr_x(      .0),  // cm
      ctr_y(      .0),  // cm
      field_size(10.)   // cm 
{}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void place_cell_generator::Parameters_::get(DictionaryDatum &d) const
{
    def<double>(d,  nest::names::rate,       rate);
    def<double>(d,  nest::names::ctr_x,      ctr_x);
    def<double>(d,  nest::names::ctr_y,      ctr_y);
    def<double>(d,  nest::names::field_size, field_size);
    def<vecType>(d, nest::names::rat_pos_x,  rat_pos_x);
    def<vecType>(d, nest::names::rat_pos_y,  rat_pos_y);
    def<double>(d,  nest::names::rat_pos_dt, rat_pos_dt);
}


void place_cell_generator::Parameters_::set(const DictionaryDatum& d)
{
    updateValue<double>(d,  nest::names::rate,       rate);
    updateValue<double>(d,  nest::names::ctr_x,      ctr_x);
    updateValue<double>(d,  nest::names::ctr_y,      ctr_y);
    updateValue<double>(d,  nest::names::field_size, field_size);
    updateValue<vecType>(d, nest::names::rat_pos_x,  rat_pos_x);
    updateValue<vecType>(d, nest::names::rat_pos_y,  rat_pos_y);
    updateValue<double>(d,  nest::names::rat_pos_dt, rat_pos_dt);

    // rat_pos must be at least one element long, if not add one
    if (rat_pos_x.size() == 0)
        rat_pos_x.push_back(.0);
    if (rat_pos_y.size() == 0)
        rat_pos_y.push_back(.0);


    if ( rate < 0 )
        throw BadProperty("Place cell model rate cannot be negative.");

    if (field_size < 0)
        throw BadProperty("Place cell field size must be >= 0.");

    if (rat_pos_dt < 0)
        throw BadProperty("Place cell positional data time step must be >= 0.");

    // Here check if rat_pos_dt is a multiple of simulation time step
    double sim_dt = Time::get_resolution().get_ms();
    if (((int) (rat_pos_dt / sim_dt)) != (rat_pos_dt / sim_dt))
        throw BadProperty("rat_pos_dt must be an integer multiple of the simulation time step.");
}



/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

place_cell_generator::place_cell_generator()
    : Node(), device_(), P_()
{
    sim_dt_per_pos_dt = 0;
    next_pos_step = 0;
    pos_it = 0;
}

place_cell_generator::place_cell_generator(const place_cell_generator& n)
    : Node(n), device_(n.device_), P_(n.P_)
{
    sim_dt_per_pos_dt = n.sim_dt_per_pos_dt;
    next_pos_step = n.next_pos_step;
    pos_it = n.pos_it;
}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void place_cell_generator::init_state_(const Node& proto)
{
    const place_cell_generator& pr = downcast<place_cell_generator>(proto);

    device_.init_state(pr.device_);
}


void place_cell_generator::init_buffers_()
{
    device_.init_buffers();
}


void place_cell_generator::calibrate()
{
    device_.calibrate();

    double sim_dt = Time::get_resolution().get_ms();

    // Calibrate how often to advance the positional index
    sim_dt_per_pos_dt = (int) (P_.rat_pos_dt / sim_dt);
    next_pos_step = sim_dt_per_pos_dt;
    pos_it = 0; 

    setFiringRate();
}


/* ----------------------------------------------------------------
 * Update function
 * ---------------------------------------------------------------- */

void place_cell_generator::update(Time const & T, const long_t from, const long_t to)
{
    assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
    assert(from < to);

    for (long_t lag = from; lag < to; ++lag)
    {
        long_t now = T.get_steps() + lag;

        if (device_.is_active(T + Time::step(lag)))
        {
            librandom::RngPtr rng = net_->get_rng(get_thread());
            ulong_t n_spikes = poisson_dev_.uldev(rng);

            if (n_spikes > 0) // we must not send events with multiplicity 0
            {
                nest::SpikeEvent se;
                se.set_multiplicity(n_spikes);
                network()->send(*this, se, lag);
            }

            // Advance the positional information if there is space in the
            // positional vector
            if (now >= next_pos_step) {
                next_pos_step = now + sim_dt_per_pos_dt;
                if (pos_it < P_.rat_pos_x.size() - 1)
                    pos_it++;
                setFiringRate();
            }
        }
    }
}


inline
double place_cell_generator::gaussianFunction()
{
    double dist_x = P_.rat_pos_x[pos_it] - P_.ctr_x;
    double dist_y = P_.rat_pos_y[pos_it] - P_.ctr_y;
    return exp(-(dist_x*dist_x + dist_y*dist_y) / 2. / (P_.field_size*P_.field_size));
}


inline
void place_cell_generator::setFiringRate()
{
    // rate is in Hz, dt in ms, so we have to convert from s to ms
    double sim_dt = Time::get_resolution().get_ms();
    poisson_dev_.set_lambda(sim_dt * P_.rate * 1e-3 * gaussianFunction());
}

} // namespace mynest
