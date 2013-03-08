/*
 *  gridcellsmodule.cpp
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

// include necessary NEST headers
#include "config.h"
#include "network.h"
#include "model.h"
#include "dynamicloader.h"
#include "genericmodel.h"
#include "generic_connector.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "nestmodule.h"

// include headers with your own stuff
#include "gridcellsmodule.h"
#include "place_cell_generator.h"
#include "iaf_gridcells.h"

// -- Interface to dynamic module loader ---------------------------------------

/*
 * The dynamic module loader must be able to find your module. 
 * You make the module known to the loader by defining an instance of your 
 * module class in global scope. This instance must have the name
 *
 * <modulename>_LTX_mod
 *
 * The dynamicloader can then load modulename and search for symbol "mod" in it.
 */

mynest::GridCellsModule gridcellsmodule_LTX_mod;


// -- DynModule functions ------------------------------------------------------

mynest::GridCellsModule::GridCellsModule()
{ 
#   ifdef LINKED_MODULE
    // register this module at the dynamic loader
    // this is needed to allow for linking in this module at compile time
    // all registered modules will be initialized by the main app's dynamic loader
    nest::DynamicLoaderModule::registerLinkedModule(this);
#   endif     
}

mynest::GridCellsModule::~GridCellsModule()
{
}


const std::string mynest::GridCellsModule::name(void) const
{
    return std::string("Grid cells NEST module"); // Return name of the module
}

const std::string mynest::GridCellsModule::commandstring(void) const
{
    /* 1. Tell interpreter that we provide the C++ part of GridCellsModule with the
       current revision number. 
       2. Instruct the interpreter to check that gridcellsmodule-init.sli exists, 
       provides at least version 1.0 of the SLI interface to GridCellsModule, and
       to load it.
       */
    return std::string(
            "/gridcellsmodule /C++ ($Revision: 9902 $) provide-component "
            "/gridcellsmodule /SLI (7165) require-component"
            );
}

//-------------------------------------------------------------------------------------

void mynest::GridCellsModule::init(SLIInterpreter *i, nest::Network*)
{
    nest::register_model<place_cell_generator>(nest::NestModule::get_network(), 
            "place_cell_generator");

    nest::register_model<nest::iaf_gridcells>(nest::NestModule::get_network(), 
            "iaf_gridcells");
}


