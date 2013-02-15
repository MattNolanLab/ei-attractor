/*
 *  ramp_current_generator.cpp
 *
 *       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
 *       
 *       This program is free software: you can redistribute it and/or modify
 *       it under the terms of the GNU General Public License as published by
 *       the Free Software Foundation, either version 3 of the License, or
 *       (at your option) any later version.
 *       
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *       GNU General Public License for more details.
 *       
 *       You should have received a copy of the GNU General Public License
 *       along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "ramp_current_generator.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"


/* ---------------------------------------------------------------- 
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */
    
nest::ramp_current_generator::Parameters_::Parameters_()
  : a(0.0), b(0.0), // pA
    c(0.0) // dimensionless
{}


/* ---------------------------------------------------------------- 
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void nest::ramp_current_generator::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, "a", a);
  def<double>(d, "b", b);
  def<double>(d, "c", c);
}  

void nest::ramp_current_generator::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, "a", a);
  updateValue<double>(d, "b", b);
  updateValue<double>(d, "c", c);
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::ramp_current_generator::ramp_current_generator()
  : Node(),
    device_(), 
    P_()
{}

nest::ramp_current_generator::ramp_current_generator(const ramp_current_generator& n)
  : Node(n), 
    device_(n.device_),
    P_(n.P_)
{}


/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::ramp_current_generator::init_node_(const Node& proto)
{
  const ramp_current_generator& pr = downcast<ramp_current_generator>(proto);

  device_.init_parameters(pr.device_);
  
  P_ = pr.P_;
}

void nest::ramp_current_generator::init_state_(const Node& proto)
{ 
  const ramp_current_generator& pr = downcast<ramp_current_generator>(proto);

  device_.init_state(pr.device_);
}

void nest::ramp_current_generator::init_buffers_()
{ 
  device_.init_buffers();
}

void nest::ramp_current_generator::calibrate()
{
  device_.calibrate();
}


/* ---------------------------------------------------------------- 
 * Update function
 * ---------------------------------------------------------------- */

void nest::ramp_current_generator::update(Time const &origin, 
                                const long_t from, const long_t to)
{
  long_t start = origin.get_steps();
  long_t t = origin.get_ms();

  CurrentEvent ce;
  ce.set_current(std::max(0.0, (P_.a * t + P_.b) * P_.c));

  for ( long_t offs = from ; offs < to ; ++offs )
    if( device_.is_active( Time::step(start+offs) ) )
      network()->send(*this, ce, offs);
}
