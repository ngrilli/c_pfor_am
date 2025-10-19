//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "c_pfor_amTestApp.h"
#include "MooseMain.h"

// Create a performance log
PerfLog Moose::perf_log("c_pfor_am");

// Begin the main program.
int
main(int argc, char * argv[])
{
  return Moose::main<c_pfor_amTestApp>(argc, argv);
}
