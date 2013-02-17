/*
 *  gridcellsmodule.h
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

#ifndef GRIDCELLSMODULE_H
#define GRIDCELLSMODULE_H

#include "dynmodule.h"
#include "slifunction.h"

namespace nest
{
  class Network;
}

// Put your stuff into your own namespace.
namespace gridcells {
  
/**
 * Class defining your model.
 * @note For each model, you must define one such class, with a unique name.
 */
class GridCellsModule : public DynModule
{
public:

  // Interface functions ------------------------------------------
  
  /**
   * @note The constructor registers the module with the dynamic loader. 
   *       Initialization proper is performed by the init() method.
   */
  GridCellsModule();
  
  /**
   * @note The destructor does not do much in modules. Proper "downrigging"
   *       is the responsibility of the unregister() method.
   */
  ~GridCellsModule();

  /**
   * Initialize module by registering models with the network.
   * @param SLIInterpreter* SLI interpreter
   * @param nest::Network*  Network with which to register models
   * @note  Parameter Network is needed for historical compatibility
   *        only.
   */
  void init(SLIInterpreter*, nest::Network*);

  /**
   * Return the name of your model.
   */
  const std::string name(void) const;
  
  /**
   * Return the name of a sli file to execute when GridCellsModule is loaded.
   * This mechanism can be used to define SLI commands associated with your
   * module, in particular, set up type tries for functions you have defined.
   */
  const std::string commandstring(void) const;
     
public:
  
  };
} // namespace gridcells

#endif
