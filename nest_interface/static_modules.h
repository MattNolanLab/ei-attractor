#ifndef STATIC_MODULES_H
#define STATIC_MODULES_H

#include "../models/modelsmodule.h"
#include "../precise/precisemodule.h"
#include "../topology/topologymodule.h"

#include "interpret.h"
#include "network.h"

void add_static_modules(SLIInterpreter& engine, nest::Network& net)
{
  engine.addmodule(new nest::ModelsModule(net));
  engine.addmodule(new nest::PreciseModule(net));
  engine.addmodule(new nest::TopologyModule(net));

}
#endif
