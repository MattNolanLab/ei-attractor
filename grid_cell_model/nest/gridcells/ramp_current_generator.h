/*
 *  ramp_current_generator.h
 *  
 *   Generator that injects ramp-like current into its destinations
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


/* 
 * Generates ramp current according to the following equation:
 *  I = max(0, (a*t + b) * c),
 *
 *  where a, b, and c are constant parameters, t is simulation time. In 
 *  other words, the generator creates a ramp-like current, the steepness of
 *  which can further be modulated by parameter c.
 */ 

#ifndef RAMP_CURRENT_GENERATOR_H
#define RAMP_CURRENT_GENERATOR_H

#include <vector>
#include "nest.h"
#include "event.h"
#include "node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "stimulating_device.h"

namespace nest
{
/**
 * Ramp current generator.
 */
class ramp_current_generator : public Node
{

    public:        

        ramp_current_generator();
        ramp_current_generator(const ramp_current_generator&);

        bool has_proxies() const {return false;} 


        port check_connection(Connection&, port);

        void get_status(DictionaryDatum &) const;
        void set_status(const DictionaryDatum &);

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
            double_t    a;   // slope of the current (pA/ms)
            double_t    b;   // constant linear part (pA)
            double_t    c;   // multiplication factor (dimensionless)

            Parameters_();  //!< Sets default parameter values

            void get(DictionaryDatum&) const;  //!< Store current values in dictionary
            void set(const DictionaryDatum&);  //!< Set values from dicitonary
        };

        // ------------------------------------------------------------

        StimulatingDevice<CurrentEvent> device_;
        Parameters_ P_;
};

inline  
port ramp_current_generator::check_connection(Connection& c, port receptor_type)
{
    CurrentEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
}

inline
void ramp_current_generator::get_status(DictionaryDatum &d) const
{
    P_.get(d);
    device_.get_status(d);
}

inline
void ramp_current_generator::set_status(const DictionaryDatum &d)
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


} // namespace

#endif /* #ifndef DC_GENERATOR_H */

